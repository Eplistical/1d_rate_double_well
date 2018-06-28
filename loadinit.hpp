#ifndef _LOADINIT_HPP
#define _LOADINIT_HPP

#include <vector>
#include <string>
#include <iostream>
#include "misc/randomer.hpp"
#include "misc/ioer.hpp"
#include "config.hpp"

inline void loadinit(std::vector<double>& x, std::vector<double>& v) 
{
    /*
    x = randomer::vnormal(para.Ntraj, para.avgx0, sqrt(para.kT0 / para.mw2));
    v = randomer::vnormal(para.Ntraj, 0.0, sqrt(para.kT0 / para.mass));
    */

    std::string initfile("./5kT100traj.dat");
    std::vector<double> xv(para.Ntraj * 2);
    x.resize(para.Ntraj);
    v.resize(para.Ntraj);
    ioer::input_t init(initfile, std::ios::in | std::ios::binary);
    int Ntraj_read;
    init.read(Ntraj_read);
    assert(Ntraj_read == para.Ntraj);
    init.read(xv);
    init.close();
    for (int i(0); i < para.Ntraj; ++i) {
        x[i] = xv[0 + i * 2];
        v[i] = xv[1 + i * 2];
    }
}

#endif // _LOADINIT_HPP
