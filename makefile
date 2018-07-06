INC += -I./dynamics_headers
OPT += -O3 -std=c++11

# Boost
BOOSTHOME = /usr/local
INC += -I$(BOOSTHOME)/include
LIBS += -L$(BOOSTHOME)/lib64 -lboost_program_options

# HDF5
HDF5HOME = /usr/local
INC += -I$(HDF5HOME)/include
LIBS += -L$(HDF5HOME)/lib -lhdf5 -lhdf5_cpp
OPT += -DIOER_WITH_HDF5

# MKL
MKLROOT = $(HOME)/intel/compilers_and_libraries_2017.4.196/linux/mkl
LIBS += -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
LIBS += -Wl,-rpath,$(MKLROOT)/lib/intel64
INC +=  -I$(MKLROOT)/include

# MPI
MPIHOME = $(HOME)/pkgs/Qchem/QC_EXT_LIBS/openmpi_icpc
MPICXX = $(MPIHOME)/bin/mpicxx
MPILIBS += -L$(MPIHOME)/lib -lmpi -lmpi_cxx -Wl,-rpath,$(MPIHOME)/lib
MPIINC += -I$(MPIHOME)/include


default: mpi
mpi: cme_mpi bcme_mpi modbcme_mpi fp_mpi bfp_mpi iesh_mpi

cme_mpi: runCME_mpi.cpp config.hpp
	$(MPICXX) $< $(INC) $(MPIINC) $(LIBS) $(MPILIBS) $(OPT) -o bin/$@

bcme_mpi: runBCME_mpi.cpp config.hpp
	$(MPICXX) $< $(INC) $(MPIINC) $(LIBS) $(MPILIBS) $(OPT) -o bin/$@

modbcme_mpi: runModBCME_mpi.cpp config.hpp
	$(MPICXX) $< $(INC) $(MPIINC) $(LIBS) $(MPILIBS) $(OPT) -o bin/$@

fp_mpi: runFP_mpi.cpp config.hpp
	$(MPICXX) $< $(INC) $(MPIINC) $(LIBS) $(MPILIBS) $(OPT) -o bin/$@

bfp_mpi: runBFP_mpi.cpp config.hpp
	$(MPICXX) $< $(INC) $(MPIINC) $(LIBS) $(MPILIBS) $(OPT) -o bin/$@

iesh_mpi: runIESH_mpi.cpp config.hpp
	$(MPICXX) $< $(INC) $(MPIINC) $(LIBS) $(MPILIBS) $(OPT) -o bin/$@

test: test.cpp 
	$(MPICXX) $< $(INC) $(MPIINC) $(LIBS) $(MPILIBS) $(OPT) -o bin/$@

showsurf: showsurf.cpp config.hpp
	$(MPICXX) $< $(INC) $(LIBS) $(OPT) -o bin/$@

exact_eql_pop: exact_eql_pop.cpp config.hpp
	$(MPICXX) $< $(INC) $(LIBS) $(OPT) -o bin/$@
