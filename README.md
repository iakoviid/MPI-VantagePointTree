# MPI-VantagePointTree
MPI-implementation of nearest neighbor search with vantage point tree.
Random sets of points are used. Dimensions is set to 2 (can change manualy).The inputs are number of points (N) , number of neighbors (K) and number of processes (p).


## Getting Started

Compilation Command

make

Compile

mpiCC mpi.c -o vp #compile

Run

mpiexec -n p ./vp N K


## Running the tests

Command mpiexec -n p ./vp N K 
example inputs p=8,N=20,K=6.

## Files

VantagePointAllK.cpp: main file creates a random set of points makes and searches a vp tree on those points
VantagePointMpiOnePoint.cpp: search function for one point
VP_Report complete explanation: of the method and code 
job.sh: bash script submited to Aristotle University machine

## Authors

Iakovidis Ioannis
