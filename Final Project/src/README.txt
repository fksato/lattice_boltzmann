# CFD Final Project: Group J
# General LBM solver
# Authors: Rohit Abraham Paul, Rajesh Maheswaran, Fukushi Sato

src contains

* a makefile
* parameter files
* headers
* LBM solvers ( Sequential and MPI )

# WARNING:
If a solver was compiled using

	make < target >

and you wish to switch targets, you MUST run
	
	make clean

before running
	
	make < other target >

# Sequential implementation:
sequential implementation of the lattice Boltzmann method

# WARNING:
if you compiled the MPI implementaion prior to compiling
the sequential implementation, you must run:
	
	make clean

make program by using the following command:

	make native

to run the solver, choose a scenario from the dat files
located in the dat directory

available scenarios for the sequential solvers are:

	artery_64
	ca_64
	movingwall_20x20x20
	movingwall_50x50x50
	side_64
	step
	step_over
	straight_artery

run solver with the following:

	./LBSim_SEQ < scenario name >

# MPI Implementation (WIP):
MPI parallel implementation of the lattice Boltzmann method

# WARNING:
if you compiled the sequential implementaion prior to compiling
the MPI implementation, you must run:
	
	make clean

**note: this solver is a work in progress and may not behave properly

make program by using the following command:

	make parallel

to run the solver, run the following command:

	mpirun -n < number of procs > ./LBMSim_MPI

this will default to the parallel moving wall scenario.

**any of the sequential scenarios can be run using the MPI implementation

	mpirun -n < number of procs > ./LBMSim_MPI < scenario name >

# WARNING:
the number of processors used must equal the product of the 
x_proc, y_proc, z_proc configurations.