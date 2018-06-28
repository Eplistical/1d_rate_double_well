#include <cassert>
#include <cmath>
#include <memory>
#include <algorithm>
#include <string>

#include "misc/vector.hpp"
#include "misc/ioer.hpp"
#include "misc/timer.hpp"
#include "misc/randomer.hpp"
#include "misc/MPIer.hpp"

#include "config.hpp"
#include "loadinit.hpp"
#include "argparse.hpp"
#include "particle_1d/bfp_particle_1d.hpp"

using ptcl_t = BFP_Particle_1D<potential_t, inttable_mgr_t>;

using std::string;
using std::vector;

int main(int argc, char** argv) 
{
    MPIer::setup();

    // parse args
    if (argparse(argc, argv, para.workdir) == false) {
        MPIer::finalize(); 
        return 0;
    }

    // setup
    string START_TIME;
    if (MPIer::master) {
        START_TIME = timer::now();
        timer::tic();
    }
    randomer::seed(MPIer::assign_random_seed(para.random_seed));
    potential_t potential(para.mass, para.omega, para.g, para.dG, para.gamma0);
    inttable_mgr_t inttable_mgr("inttable.hdf5", potential);
    for (int r(0); r < MPIer::size; ++r) {
        if (MPIer::rank == r) {
            inttable_mgr.load_inttable();
            assert(inttable_mgr.kT == para.kT);
            assert(inttable_mgr.gamma0 == para.gamma0);
        }
        MPIer::barrier();
    }

    // prep init particles
    vector<double> x, v;
    if (MPIer::master) {
        loadinit(x, v);
    }
    MPIer::bcast(0, x, v);

    // assign jobs
    vector<uint32_t> mybatch = MPIer::assign_job(para.Ntraj);
    const int Njob(mybatch.size());

    // recorders
    const int Nrecord(static_cast<int>(para.Nstep / para.Anastep + 0.5));
    vector<double> sumN(Nrecord, 0.0);
    vector<double> sumEk(Nrecord, 0.0);

    // evolution
    for (int ijob(0); ijob < Njob; ++ijob) {
        int j(mybatch[ijob]);
        ptcl_t ptcl(x[j], v[j], para.mass, para.kT, para.nuclear_fric, potential, inttable_mgr);
        for (size_t istep(0); istep < para.Nstep; ++istep) {
            if (istep % para.Anastep == 0) {
                const int irecord(istep / para.Anastep);
                sumN[irecord] += ptcl.get_N();
                sumEk[irecord] += ptcl.get_Ek();
            }
            ptcl.evolve(para.dt);
        }
    }

    // gather data
    MPIer::barrier();
    for (int r(1); r < MPIer::size; ++r) {
        if (MPIer::master) {
            vector<double> buf;

            MPIer::recv(r, buf);
            sumN = sumN + buf;
            MPIer::recv(r, buf);
            sumEk = sumEk + buf;
        }
        else if (MPIer::rank == r) {
            MPIer::send(0, sumN);
            MPIer::send(0, sumEk);
        }
        MPIer::barrier();
    }

    // output
    if (MPIer::master) {
        string outfile_fullpath = para.workdir + "/bfp_mpi.out";
        ioer::info("outfile: ", outfile_fullpath);
        ioer::output_t out(outfile_fullpath);
        out.set_precision(10);
        out.set_width(20);
        out.info("# ", START_TIME);

        out.tabout("# t", "N", "Ek/kT");
        for (int i(0); i < Nrecord; ++i) {
            out.tabout( i * para.Anastep * para.dt, 
                    1.0 - sumN[i] / para.Ntraj,
                    sumEk[i] / para.Ntraj / para.kT
                    );
        }
        out.info("# MPI_size: ", MPIer::size);
        out.info("# ", timer::toc());
        out.info("# ", timer::now());
        out.close();
    }

    MPIer::finalize();
    return 0;
}
