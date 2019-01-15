
analysis for Micromegas data with APV-FEC readout

data has to be in root-tree format of mmdaq software

additional data from Munich CRF can be evaluated


for the data format have look at the variables in the analysis class in the analysis.h file

this header file provide all needed variables and basic functions


for compilation root and opencv is required (opencv can be removed if two functions are removed manually)

the main functions are in the fitter.C , investigateCRF.C and postprocessor.C files

the executables will display their usage on commandline (for executable names look into the Makefile)


sample parameterfiles, required for all analysis steps, are in the parameterfiles subfolders

the analysis highly depends on these parameterfiles


other scripts are mainly for distribution of slurm jobs on local grid and for making summary plots

