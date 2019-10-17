# cluster-pipeline
README for cluster version of CHILES Pipeline
October 16, 2019
Version 4.5

This is the same as the CHILES pipeline in terms of how it does data reduction 
(see github.com/djpisano3/chiles_pipeline for more details), except it has 
been designed to run on the WVU HPC Spruce Knob.  

The code has been tested and designed for CASA 5.3.0.  It may work with 
older or newer versions, but no guarantee of results is made due to 
CASA changes.

Notes on the different pipeline versions can be found in CHILES_pipe_initial.py.

There are four sets of files in this distribution:

1) Utility code: EVLA_functions.py, lib_EVLApipeutils.py, and chiles_paramters.txt

2) Reduction code: Those codes with "pipe" in their name.  These do the work.

3) Pipeline code:  These codes run some or all of the pipeline and contain "pipeline" in their names.

4) PBS cluster scripts:  Shell scripts for submitting jobs.

In this version, chiles_parameters.txt contains values for running the pipeline
only on epoch 1 & 2 data.  

Note that the "split" code will convert the html output to PDF files, but it will only run on Spruce Knob.  
