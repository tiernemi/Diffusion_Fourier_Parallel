#!/bin/sh
#SBATCH -n 16           # 8 CPU cores
#SBATCH -t 0-00:30:00   # 1 hour 
#SBATCH -p debug      # partition name
#SBATCH -J FFTPara     # sensible name for the job


# load up the correct modules, if required
module load apps libs cports gcc/4.9.3-gnu fftw/3.3.4_mpi-gnu4.9.3 openmpi/1.8.6-gnu4.9.3 
# launch the code
mpirun -n 16 ./diff -n 400 -m 400

