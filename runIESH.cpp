#include <cassert>
#include <cmath>
#include <memory>
#include <algorithm>
#include <string>

#include "misc/vector.hpp"
#include "misc/ioer.hpp"
#include "misc/timer.hpp"
#include "misc/randomer.hpp"

#include "config.hpp"
#include "loadinit.hpp"
#include "hamiltonian_1d/impurity_near_uniform_band_hamiltonian.hpp"
#include "particle_1d/iesh_particle_1d.hpp"

using hamiltonian_t = ImpurityNearUniformBandHamiltonian_1D<potential_t>;
using ptcl_t = IESH_Particle_1D<hamiltonian_t>;

using std::string;
using std::vector;

int main(int argc, char** argv) 
{
    string START_TIME(timer::now());
    timer::tic();

    // setup
    ioer::output_t out("iesh.out");
    out.set_precision(10);
    out.set_width(20);
    out.info("# ", START_TIME);
    randomer::seed(para.random_seed);
    const int Nele(static_cast<int>(para.Nbath / 2));
    const int Nhole(Nele + 1);
    potential_t potential(para.mass, para.omega, para.g, para.dG, para.gamma0);
    hamiltonian_t hamiltonian(potential, para.bandwidth, Nele + Nhole, Nele + Nhole - 1, para.gamma0);

    // load init trajectories
    vector<double> x, v;
    loadinit(x, v);

    vector<ptcl_t> swarm;
    swarm.reserve(para.Ntraj);
    for (int j(0); j < para.Ntraj; ++j) {
        swarm.push_back(ptcl_t(x[j], v[j], para.mass, para.kT, para.nuclear_fric, Nele, Nhole, para.Ndtq, para.thermal_tau, hamiltonian));
    }

    // evolution
    for (size_t istep(0); istep < para.Nstep; ++istep) {
        double sumN(0.0);
        double sumEk(0.0);

        for (ptcl_t& ptcl : swarm) {
            ptcl.evolve(para.dt);
            sumN += ptcl.get_N();
            sumEk += ptcl.get_Ek();
        }

        if (istep % para.Anastep == 0) {
            out.tabout(istep * para.dt, 
                        1.0 - sumN / para.Ntraj, 
                        sumEk / para.Ntraj / para.kT);
        }
    }

    // done
    out.info("# ", timer::toc());
    out.info("# ", timer::now());
    out.close();

    return 0;
}
