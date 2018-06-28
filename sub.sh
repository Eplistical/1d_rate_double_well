#!/bin/bash

exelist="
iesh
"
#cme
#bcme
#fp
#bfp

for exe in $exelist; do

cat > .sub${exe}.sh << EOF
#PBS -N ${exe}
#PBS -j oe
#PBS -o .${exe}.log
#PBS -l walltime=1214:00:00
#PBS -l mem=2GB
#PBS -l nodes=1:ppn=1

cd \$PBS_O_WORKDIR
./bin/${exe}
EOF

qsub .sub${exe}.sh

done
