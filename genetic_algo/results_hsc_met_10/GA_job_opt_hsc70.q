#!/bin/bash
#$ -N job_opt_hsc70
#$ -cwd
#$ -V
#$ -q all.q
#$ -l h_vmem=2G


module load python/3.5.1

python3.5 GA_opt_hsc70.py