#!/bin/bash
#PBS -q see
#PBS -N mpi
#PBS -j oe
#PBS -l nodes=4:ppn=8

module load gcc/7.2.0
module load openmpi/3.0.0


cd $PBS_O_WORKDIR

mpiCC mpi.c -o vp #compile

mpiexec ./vp 20 16  >> nodes_4.txt       #run
