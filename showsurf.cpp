#include <cassert>
#include <cmath>
#include <cassert>
#include <memory>
#include <algorithm>
#include <string>

#include "misc/vector.hpp"
#include "misc/ioer.hpp"
#include "misc/timer.hpp"
#include "misc/randomer.hpp"

#include "config.hpp"
#include "loadinit.hpp"
#include "argparse.hpp"
#include "hamiltonian_1d/impurity_near_uniform_band_hamiltonian.hpp"
#include "particle_1d/iesh_particle_1d.hpp"

using hamiltonian_t = ImpurityNearUniformBandHamiltonian_1D<potential_t>;
using ptcl_t = IESH_Particle_1D<hamiltonian_t>;

std::string outfile = "showsurf.out";

int main(int argc, char** argv) 
{
    // setup
    randomer::seed(para.random_seed);
    ioer::output_t out(outfile);
    out.set_precision(6);
    const int Nele(static_cast<int>(para.Nbath / 2));
    const int Nhole(Nele + 1);
    potential_t potential(para.mass, para.omega, para.g, para.dG, para.gamma0);
    hamiltonian_t hamiltonian(potential, para.bandwidth, Nele + Nhole, Nele + Nhole - 1, para.gamma0);


    double xmin(-10.0);
    double xmax(30.0);
    const int N(2000);
    const double dx((xmax - xmin) / N);

    for (double x = xmin; x <= xmax; x += dx) {
        hamiltonian.update_H(x);
        hamiltonian.update_dc(x);
        out.tabout( x, 
                hamiltonian.cal_potential(x, 0), 
                hamiltonian.cal_potential(x, 1), 
                hamiltonian.cal_force(x, 0), 
                hamiltonian.cal_force(x, 1),
                hamiltonian.eva.at(1),
                hamiltonian.eva.at(4),
                hamiltonian.eva.at(9),
                hamiltonian.dc.at(1 + 2 * hamiltonian.dim),
                hamiltonian.dc.at(3 + 5 * hamiltonian.dim),
                hamiltonian.dc.at(4 + 7 * hamiltonian.dim)
              );
    }
    out.close();

    return 0;
}
