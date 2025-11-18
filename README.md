# README
### The geometry of pMHC-coated nanoparticles and T cell receptor clusters governs the sensitivity-specificity trade-off in T cell response: A modeling investigation
---------------------------------------------------------------------------------------
### Folders

- Capacity Results/
	
	- Folder containing simulation results quantifying the surface capacity of bound
NPs as well as distirbution of covered TCRs corresponding to results presented
in Figure 1.

- Functions/
	
	- Folder containing functions required for MC simulations and analysis.
	To run simulations, make sure this folder is in your MATLAB path with
	the command: 		addpath(Functions/)
	
- LongSims/
	
	- Functions for running simulations: DoseResponseFunction.m
		Takes arguments 1-16 corresponding to index of 16 different NP
		concentration values e.g. "DoseResponseFunction(10)". Simulation
		results will be saved to a temp directory.
		
	- Simulation results are stored using the following file structure:
		NP Radius/TCRs per Cluster/NP valence/koff/filename_NP_concentration.mat
		
- Mutual Information/
	- Results from our mutual information analysis. Mutual information was calculated
	for Kd-, valence- and concentration- (rho) based discrimination. Results presented
	in manuscript Figure 2d and Figure 3c correspond to "kd_mi.mat" and "valence_mi.mat"
	respectively.
	
- Plotting/
	- Scripts and functions for generating manuscript Figures can be found as Figure*.m
	labelled according to order of appearance in the manuscript.
	

All files should be run from the parent folder Code/ as the working directory.
