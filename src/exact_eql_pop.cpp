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

    const int Nele(static_cast<int>(para.Nbath / 2));
    const int Nhole(Nele + 1);
    potential_t potential(para.mass, para.omega, para.g, para.dG, para.gamma0);

    assert(not para.inttablefile.empty());
    inttable_mgr_t inttable_mgr(para.inttablefile, potential);
    inttable_mgr.load_inttable();
    assert(inttable_mgr.kT == para.kT);
    assert(inttable_mgr.gamma0 == para.gamma0);

    const double xmin(-100.0);
    const double xmax(100.0);
    const int N(100000);
    const double dx((xmax - xmin) / N);

    double Uadiab(potential.cal_potential(xmin, 0));

    double int_nx(0.0);
    double numer(0.0), denum(0.0);

    for (int i(0); i < N; ++i) {
        const double x(xmin + i * dx);
        double dhdx(potential.cal_dhdx(x));
        double f(misc::fermi(potential.cal_h(x) / para.kT));
        double n(dhdx * f + inttable_mgr.retrieve("n", x));
        int_nx += n * dx;
        Uadiab = potential.cal_potential(x, 0) + dhdx * int_nx;
        numer += n * exp(-Uadiab / para.kT) * dx;
        denum += exp(-Uadiab / para.kT) * dx;
    }
    ioer::tabout("# gamma0", "kT");
    ioer::tabout(para.gamma0, para.kT);
    ioer::tabout("# N0", "1 - N0");
    ioer::tabout(numer / denum, 1.0 - numer / denum);

    return 0;
}
