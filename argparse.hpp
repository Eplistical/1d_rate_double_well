#ifndef _ARGPARSE_HPP
#define _ARGPARSE_HPP

#include "boost/program_options.hpp"
#include <iostream>
#include <string>

namespace po = boost::program_options;

inline bool argparse(int argc, char** argv, 
                std::string& workdir) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("workdir", po::value<std::string>(&workdir)->default_value("."), "working directory")
        ;
    po::variables_map vm; 
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return false;
    }

    return true;
}


#endif // _ARGPARSE_HPP
