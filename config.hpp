#ifndef _CONFIG_HPP
#define _CONFIG_HPP

#include "potential_1d/double_harmonic_potential_const_gamma_1d.hpp"
#include "inttable_mgr.hpp"
#include "misc/ioer.hpp"
#include <vector>
#include <string>

struct Para {
    const double gamma0 = 6.4e-3;
    const double omega = 2.0e-4;
    const double mass = 2000.0;
    const double g = 20.6097;
    const double kT = 9.500432590749929e-4;
    const double dG = -0.0038;
    const double nuclear_fric = 0.0;

    const int Nbath = 50; 
    const double bandwidth = 6.4e-2;
    const double thermal_tau = -1.0;

    const double dt = 1.0;
    const int Ndtq = 1;
    const size_t Nstep = 1e4;
    const size_t Anastep = 200;
    const int Ntraj = 100;

    const double mw2 = mass * omega * omega;
    const int random_seed = 22906779;

    const double kT0 = 5.0 * kT;
    const double avgx0 = 0;
    const int surf0 = 0;

    std::string workdir;
} para;

inline void saveParatoh5(ioer::h5file_t& h5f) {
    h5f.create_dataset("para", std::vector<double>(1, 0.0));
    h5f.create_attr("para", 
            "gamma0", para.gamma0,
            "omega", para.omega,
            "mass", para.mass,
            "g", para.g,
            "kT", para.kT,
            "dG", para.dG,
            "nuclear_fric", para.nuclear_fric,
            "Nbath", para.Nbath,
            "bandwidth", para.bandwidth,
            "thermal_tau", para.thermal_tau,
            "dt", para.dt,
            "Ndtq", para.Ndtq,
            "Nstep", para.Nstep,
            "Anastep", para.Anastep,
            "Ntraj", para.Ntraj,
            "random_seed", para.random_seed,
            "kT0", para.kT0,
            "avgx0", para.avgx0,
            "surf0", para.surf0
            );
}

using potential_t = DoubleHarmonicPotentialConstGamma_1D;
using inttable_mgr_t = InttableMgr<potential_t>;
using para_t = Para;


#endif // _CONFIG_HPP
