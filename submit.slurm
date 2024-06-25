#!/bin/bash

#SBATCH -J scan-xtb
#SBATCH -p 64c512g
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH --exclusive

param=$1
parallel -j 8 -a $1 ./parallel_run.sh {1}

