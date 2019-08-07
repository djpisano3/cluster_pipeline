#!/bin/bash

#PBS -q comm_mmem_day
#PBS -N batch_QA1
#PBS -m ae
#PBS -M djpisano@mail.wvu.edu
#PBS -l nodes=1:ppn=1,pvmem=54g
#PBS -n 
#PBS -e QA1.error
#PBS -o QA1.output
#PBS -l prologue=/users/djpisano/prologue.sh
#PBS -l epilogue=/users/djpisano/epilogue.sh

# Listing of sessions for testing: 1,2,5,8,9,10,16,17

#PBS -t 1,2,5,8,9,10,16,17

module load astronomy/casa/5.3.0

/usr/bin/time -v casa --nogui --nologger --agg -c "var=[${PBS_ARRAYID}];execfile('/users/djpisano/cluster_pipeline/CHILES_pipeline_QA1_batch.py')"
