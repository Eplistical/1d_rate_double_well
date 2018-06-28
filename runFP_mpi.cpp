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
#include "particle_1d/fp_particle_1d.hpp"

using ptcl_t = FP_Particle_1D<potential_t>;
using std::string;
using std::vector;

int main(int argc, char** argv) 
{
    MPIer::setup();

    // parse args
    if (argparse(MPIer::master, argc, argv, para.workdir) == false) {
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

    // load init trajectories
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
        ptcl_t ptcl(x[j], v[j], para.mass, para.kT, para.nuclear_fric, potential);
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
        vector<double> tarr(Nrecord, 0.0);
        for (int irecord(0); irecord < Nrecord; ++irecord) {
            tarr[irecord] = irecord * para.Anastep * para.dt; 
        }

        string outfile_fullpath = para.workdir + "/fp_mpi.h5";
        ioer::info("outfile: ", outfile_fullpath);
        ioer::h5file_t h5f(outfile_fullpath, std::ios::out);
        saveParatoh5(h5f);
        h5f.create_dataset(
                "t", tarr,
                "N0", 1.0 - sumN / para.Ntraj,
                "Ek", sumEk / para.Ntraj / para.kT
                );
        h5f.create_attr("para", 
                "MPI_size", MPIer::size,
                "start_time", START_TIME,
                "end_time", timer::now(),
                "elapsed_time", timer::toc()
                );
        h5f.close();
    }

    MPIer::finalize();
    return 0;
}
