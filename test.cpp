#include <cassert>
#include <cmath>
#include <memory>
#include <algorithm>
#include <string>
#include <iostream>

#include "misc/vector.hpp"
#include "misc/ioer.hpp"
#include "misc/timer.hpp"
#include "misc/randomer.hpp"
#include "misc/MPIer.hpp"

#include "config.hpp"
#include "loadinit.hpp"
#include "argparse.hpp"
#include "particle_1d/cme_particle_1d.hpp"

#include "boost/program_options.hpp"

using ptcl_t = CME_Particle_1D<potential_t>;
using std::string;
using std::vector;
std::string outfile = "cme_mpi.out";

int main(int argc, char** argv) 
{
    MPIer::setup();

    MPIer::finalize();
    return 0;
}
