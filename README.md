# Survival-Probabilities-on-Branching-Tissue
Simulation code and raw data for "Survival in Branching Cellular Populations" by A. S. Bryant and M. O. Lavrentovich from the Department of Physics and Astronomy at the University of Tennessee, Knoxville.

The code and data directory files are described as follows:

Code
	*.h : header files necessary to compile the programs
	Perlen.cpp : statistics calculations for final branched structures at 1000 generations, generates raw data file (fixed N_0 and b)
	BARW_bt.cpp : statistics calculations for branched structures as a function of time, generates raw data file
	BARW_bs_nd.cpp : Phase diagram generation for survival probability as a function of selective advantage s and branching rate b
	BARW_Nb.cpp : Phase diagram generation for survival probability as a function of branch population size N_0 and branching rate b
	modelplot.cpp : Generates a VTK file for a single simulation run of a branched structure. Can be used to visualize the branched structures.

Data
	This directory contains subdirectories containing all raw data analyzed in the project. Each subdirectory contains a readme describing the format of the raw data.
	