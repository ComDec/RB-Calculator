#!/bin/bash

#SBATCH -J scan-xtb
#SBATCH -p 64c512g
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH --exclusive

conda activate torch
python extract_xtb_results.py

