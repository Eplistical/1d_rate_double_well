#!/bin/bash

exelist="
iesh
"

#bfp
#bcme
#cme
#fp

for exe in $exelist; do

cat > .sub${exe}_mpi.sh << EOF
#PBS -N ${exe}_mpi
#PBS -j oe
#PBS -o .sub${exe}_mpi.log
#PBS -l walltime=1214:00:00
#PBS -l mem=60GB
#PBS -l nodes=4:ppn=25

cd \$PBS_O_WORKDIR
mpirun \
    -n 100 \
    ./bin/${exe}_mpi
EOF
#-n \${PBS_NP} \

qsub .sub${exe}_mpi.sh

done
