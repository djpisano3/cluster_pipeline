#!/bin/bash

#PBS -q comm_mmem_day
#PBS -N batch_QA2
#PBS -m ae
#PBS -M djpisano@mail.wvu.edu
#PBS -l nodes=1:ppn=1,pvmem=54g
#PBS -n 
#PBS -e QA2.error
#PBS -o QA2.output
#PBS -l prologue=/users/djpisano/prologue.sh
#PBS -l epilogue=/users/djpisano/epilogue.sh

# Listing of sessions for testing: 1,2,5,8,9,10,16,17

#PBS -t 1

module load astronomy/casa/5.3.0

/usr/bin/time -v mpicasa -n 9 /shared/software/astronomy/casa/5.3.0-143.el6/bin/casa --nogui --nologger --agg -c "var=[${PBS_ARRAYID}];execfile('/users/djpisano/cluster_pipeline/CHILES_pipeline_QA2_batch.py')"
