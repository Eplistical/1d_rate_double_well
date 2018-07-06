# 1d double well crossing rate, written in C++11

**Algorithms included**

- FP
- BFP
- CME
- BCME
- ModBCME
- IESH 

**Prerequisite**

- [cmake](https://cmake.org/)
- [MPI](https://www.open-mpi.org/)
- [Intel MKL](https://software.intel.com/en-us/mkl)
- [HDF5](https://support.hdfgroup.org/HDF5/)
- [Boost](https://www.boost.org/)
- [misc](https://github.com/Eplistical/misc)


**How to install**

```
git clone https://github.com/Eplistical/1d_rate_double_well.git

cd 1d_rate_double_well

mkdir build

cmake ..

make 
```


**How to use**

```
cd 1d_rate_double_well/build

./bin/iesh_mpi --help
```
