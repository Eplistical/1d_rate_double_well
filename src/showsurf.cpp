#include <cassert>
#include <cmath>
#include <cassert>
#include <memory>
#include <algorithm>
#include <string>

#include "misc/vector.hpp"
#include "misc/ioer.hpp"
#include "misc/timer.hpp"
#include "misc/fermi.hpp"
#include "misc/randomer.hpp"

#include "config.hpp"
#include "loadinit.hpp"
#include "argparse.hpp"

int main(int argc, char** argv) 
{
    // setup
    if (argparse(true, argc, argv) == false) {
        return 0;
    }
    assert(not para.outfile.empty());

    const int Nele(static_cast<int>(para.Nbath / 2));
    const int Nhole(Nele + 1);
    potential_t potential(para.mass, para.omega, para.g, para.dG, para.gamma0);

    assert(not para.inttablefile.empty());
    inttable_mgr_t inttable_mgr(para.inttablefile, potential);
    inttable_mgr.load_inttable();
    assert(inttable_mgr.kT == para.kT);
    assert(inttable_mgr.gamma0 == para.gamma0);

    const double xmin(-10.0);
    const double xmax(30.0);
    const int N(2000);
    const double dx((xmax - xmin) / N);

    vector<double> xarr(N, 0.0);
    vector<double> U0(N, 0.0), U1(N, 0.0);
    vector<double> Ub0(N, 0.0), Ub1(N, 0.0);
    xarr[0] = xmin;
    U0[0] = potential.cal_potential(xmin, 0);
    U1[0] = potential.cal_potential(xmin, 1);
    Ub0[0] = U0[0];
    Ub1[0] = U1[0];
    for (int i(1); i < N; ++i) {
        double x = xmin + i * dx;
        xarr[i] = x;
        U0[i] = U0[i - 1] - dx * potential.cal_force(x, 0);
        U1[i] = U1[i - 1] - dx * potential.cal_force(x, 1);

        double f(misc::fermi(potential.cal_h(x) / para.kT));
        double dhdx(potential.cal_dhdx(x));
        double fBCME(dhdx * f + inttable_mgr.retrieve("force", x));
        Ub0[i] = Ub0[i - 1] - dx * (potential.cal_force(x, 0) + fBCME);
        Ub1[i] = Ub1[i - 1] - dx * (potential.cal_force(x, 1) + fBCME);
    }
    ioer::h5file_t h5f(para.outfile, std::ios::out);
    saveParatoh5(h5f);
    h5f.create_dataset(
            "x", xarr,
            "U0", U0,
            "U1", U1,
            "Ub0", Ub0,
            "Ub1", Ub1
            );
    h5f.close();

    return 0;
}
