#!/bin/bash

#PBS -q comm_mmem_day
#PBS -N batch_split
#PBS -m ae
#PBS -M djpisano@mail.wvu.edu
#PBS -l nodes=1:ppn=1,pvmem=54g
#PBS -n 
#PBS -e split.error
#PBS -o split.output
#PBS -l prologue=/users/djpisano/prologue.sh
#PBS -l epilogue=/users/djpisano/epilogue.sh

# Listing of sessions for testing: 1,2,5,8,9,10,16,17

#PBS -t 1

module load astronomy/casa/5.3.0

module load singularity/2.5.2

/usr/bin/time -v casa --nogui --nologger --agg -c "var=[${PBS_ARRAYID}];execfile('/users/djpisano/cluster_pipeline/CHILES_pipeline_split_batch.py')"
