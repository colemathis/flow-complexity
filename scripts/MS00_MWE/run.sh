#!/bin/bash
#SBATCH -t 0-1
#SBATCH --job-name=test-array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-100
#SBATCH -e /dev/null
#SBATCH -o /dev/null

module load julia
ID=${SLURM_ARRAY_TASK_ID}
mkdir -p ../../data/sims/MS00_MWE/logs
julia launch.jl sims/MS00_MWE/array.csv ${ID} >> ../../data/sims/MS00_MWE/logs/launch.jl.${ID}.log
