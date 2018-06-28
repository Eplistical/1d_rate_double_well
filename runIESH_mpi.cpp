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
#include "hamiltonian_1d/impurity_near_uniform_band_hamiltonian.hpp"
#include "particle_1d/iesh_particle_1d.hpp"

using hamiltonian_t = ImpurityNearUniformBandHamiltonian_1D<potential_t>;
using ptcl_t = IESH_Particle_1D<hamiltonian_t>;

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
    const int Nele(static_cast<int>(para.Nbath / 2));
    const int Nhole(Nele + 1);
    potential_t potential(para.mass, para.omega, para.g, para.dG, para.gamma0);
    hamiltonian_t hamiltonian(potential, para.bandwidth, Nele + Nhole, Nele + Nhole - 1, para.gamma0);

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
    vector<double> sumx(Nrecord, 0.0);
    vector<double> sumv(Nrecord, 0.0);
    vector<double> sumEk(Nrecord, 0.0);
    vector<double> sumnu_Ep(Nrecord, 0.0);
    vector<double> sumel_Ep(Nrecord, 0.0);

    // evolution
    uint64_t hop_count(0);
    for (int ijob(0); ijob < Njob; ++ijob) {
        int j(mybatch[ijob]);
        ptcl_t ptcl(x[j], v[j], para.mass, para.kT, para.nuclear_fric, Nele, Nhole, para.Ndtq, para.thermal_tau, hamiltonian);
        for (size_t istep(0); istep < para.Nstep; ++istep) {
            if (istep % para.Anastep == 0) {
                const int irecord(istep / para.Anastep);
                sumN[irecord] += ptcl.get_N();
                sumx[irecord] += ptcl.x;
                sumv[irecord] += ptcl.v;
                sumEk[irecord] += ptcl.get_Ek();
                sumnu_Ep[irecord] += ptcl.get_nu_Ep();
                sumel_Ep[irecord] += ptcl.get_el_Ep();
            }
            ptcl.evolve(para.dt);
        }
        hop_count += ptcl.hop_count;
    }
    MPIer::barrier();

    // gather data
    vector< vector<double>* > data_ptr 
    {
        &sumN, &sumx, &sumv, &sumEk, &sumnu_Ep, &sumel_Ep
    };

    for (int r(1); r < MPIer::size; ++r) {
        if (MPIer::master) {
            vector<double> buf;
            for (auto& it : data_ptr) {
                MPIer::recv(r, buf);
                *it = *it + buf;
            }

            int ibuf;
            MPIer::recv(r, ibuf);
            hop_count += ibuf;
        }
        else if (MPIer::rank == r) {
            for (auto& it : data_ptr) {
                MPIer::send(0, *it);
            }
            MPIer::send(0, hop_count);
        }
        MPIer::barrier();
    }

    // output
    if (MPIer::master) {
        string outfile_fullpath = para.workdir + "/iesh_mpi.out";
        ioer::info("outfile: ", outfile_fullpath);
        ioer::output_t out(outfile_fullpath);
        out.set_precision(10);
        out.set_width(20);
        out.info("# ", START_TIME);

        out.tabout("# t", "1 - N0", "x", "v", "Ek/kT", "nu_Ep/kT", "el_Ep/kT", "Ep/kT", "Etot/kT");
        for (int i(0); i < Nrecord; ++i) {
            out.tabout( i * para.Anastep * para.dt, 
                    1.0 - sumN[i] / para.Ntraj,
                    sumx[i] / para.Ntraj,
                    sumv[i] / para.Ntraj,
                    sumEk[i] / para.Ntraj / para.kT,
                    sumnu_Ep[i] / para.Ntraj / para.kT,
                    sumel_Ep[i] / para.Ntraj / para.kT,
                    (sumnu_Ep[i] + sumel_Ep[i]) / para.Ntraj / para.kT,
                    (sumEk[i] + sumnu_Ep[i] + sumel_Ep[i]) / para.Ntraj / para.kT
                    );
        }
        out.info("# hop_count: ", hop_count);
        out.info("# MPI_size: ", MPIer::size);
        out.info("# ", timer::toc());
        out.info("# ", timer::now());
        out.close();
    }

    MPIer::finalize();
    return 0;
}
