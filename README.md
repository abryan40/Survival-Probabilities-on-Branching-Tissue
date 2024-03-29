# Survival-Probabilities-on-Branching-Tissue
Simulation code and raw data for [*Survival in Branching Cellular Populations*](https://arxiv.org/abs/2108.04992) by A. S. Bryant and M. O. Lavrentovich from the Department of Physics and Astronomy at the University of Tennessee, Knoxville

The code and data directory files are described below.

## Code Directory

The contents have all of the necessary C++ code to generate the raw data and structure visualizations of the branching simulations. Different programs correspond to different data types produced:

*X.h* : Header files necessary to compile the programs
	
*Perlen.cpp* : Generates a file with orientational correlations for the persistence length calculations. A single branch without collisions is simulated.
	
*BARW_bt.cpp* : Statistics calculations for branched structures as a function of time. Generates raw data file of time-dependent data.
	
*BARW_bs_nd.cpp* : Phase diagram generation for survival probability and other probabilities related to evolutionary dynamics of a mutant as a function of selective advantage *s* and branching rate *b*. Only structures which do not completely die out are considered.
	
*BARW_Nb.cpp* : Phase diagram generation for various branching structure characteristics (number of actively dividing branches, number of terminal branches, etc.) as a function of branch population size *N<sub>0</sub>* and branching rate *b*. Generates raw data file.
	
*modelplot.cpp* : Generates a VTK file for a single simulation run of a branched structure. Can be used to visualize the branched structures.
	
Any of the *X.cpp* files may be compiled via:   
	`g++ -std=c++11 -ggdb -O3 -pthread -o PROGRAM_NAME`

The generated programs require no inputs, so the individual *.cpp* files have to be edited in order to change parameters of the simulation.
	
Check the individual *.cpp* files for additional comments on how to set parameters of the simulation.

## Data

This directory contains subdirectories containing all raw data analyzed in the project. Each subdirectory contains a README describing the format of the raw data. These data were produced using the code in the Code directory.
	
## Branching structure visualization

To visualize the branching structures, the *modelplot.cpp* code has to be compiled, with appropriate choices of parameters within the file (see comments in the code). Running the compiled program will generate a VTK file which can be opened via [Paraview](https://www.paraview.org/download/) for visualization of the structure. The file generated contains 3-vector 1-scalar format data for 3-vector position and a 0.0-1.0 scalar for a fractional representation of the generation number compared to the desired final generation. This quantity is mapped to color to indicate the oldest to newest generations. The VTK file may be read in via the *Particles Reader* in ParaView. The repository contains an *EXAMPLE.vtk* structure file generated by the simulation which can be opened with ParaView as described. Here is the example structure:
	
![EXAMPLE.vtk structure](/example.png)
	
	
	
