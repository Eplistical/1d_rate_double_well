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
#include "inttable_mgr.hpp"
#include "particle_1d/bcme_with_ele_fric_diff_particle_1d.hpp"

using ptcl_t = BCME_Particle_1D<potential_t, inttable_mgr_t>;

int main(int argc, char** argv) 
{
    string START_TIME(timer::now());
    timer::tic();

    // setup
    randomer::seed(para.random_seed);
    ioer::output_t out("bcme.out");
    out.set_precision(10);
    out.set_width(20);
    out.info("# ", START_TIME);
    potential_t potential(para.mass, para.omega, para.g, para.dG, para.gamma0);
    inttable_mgr_t inttable_mgr("inttable.hdf5", potential);
    inttable_mgr.load_inttable();
    assert(inttable_mgr.kT == para.kT);
    assert(inttable_mgr.gamma0 == para.gamma0);

    // load init trajectories
    vector<double> x, v;
    loadinit(x, v);

    vector<ptcl_t> swarm;
    swarm.reserve(para.Ntraj);
    for (int j(0); j < para.Ntraj; ++j) {
        swarm.push_back(ptcl_t(x[j], v[j], para.mass, para.kT, para.surf0, para.nuclear_fric, potential, inttable_mgr));
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
            out.tabout(istep * para.dt, 1.0 - sumN / para.Ntraj, sumEk / para.Ntraj / para.kT);
        }
    }

    // done
    out.info("# ", timer::toc());
    out.info("# ", timer::now());
    out.close();

    return 0;
}
