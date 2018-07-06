#ifndef _ARGPARSE_HPP
#define _ARGPARSE_HPP

#include <iostream>
#include <string>
#include "boost/program_options.hpp"
#include "config.hpp"

namespace po = boost::program_options;

inline bool argparse(bool output_flag, int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("nuclear_fric", po::value<double>(&para.nuclear_fric), "nuclear friction")
        ("gamma0", po::value<double>(&para.gamma0), "impurity-bath coupling")
        ("bandwidth", po::value<double>(&para.bandwidth), "bandwith (measured by gamma0)")
        ("Nbath", po::value<int>(&para.Nbath), "orbital number in the bath")
        ("thermal_tau", po::value<double>(&para.thermal_tau), "el thermostat parameter")
        ("Ntraj", po::value<int>(&para.Ntraj), "# trajectories")
        ("dt", po::value<double>(&para.dt), "time interval each step")
        ("Ndtq", po::value<int>(&para.Ndtq), "classical time / quantum time")
        ("Nstep", po::value<size_t>(&para.Nstep), "total simulation steps")
        ("Anastep", po::value<size_t>(&para.Anastep), "# step for each record")
        ("kT0", po::value<double>(&para.kT0), "kT for initialization (measured by kT)")
        ("avgx0", po::value<double>(&para.avgx0), "initial <x>")
        ("initfile", po::value<std::string>(&para.initfile), "init file name")
        ("outfile", po::value<std::string>(&para.outfile), "output file name")
        ("inttablefile", po::value<std::string>(&para.inttablefile), "integral table file name")
        ;
    po::variables_map vm; 
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        if (output_flag)
            std::cout << desc << "\n";
        return false;
    }

    if (vm.count("bandwidth")) {
        para.bandwidth *= para.gamma0;
    }

    if (vm.count("kT0")) {
        para.kT0 *= para.kT;
    }

    return true;
}


#endif // _ARGPARSE_HPP
