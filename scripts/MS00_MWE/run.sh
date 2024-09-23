#!/bin/bash
#SBATCH --job-name=MS00
#SBATCH --partition=htc
#SBATCH --qos=public
#SBATCH --array=1-1000
#SBATCH --time=01:00:00
#SBATCH --output=./logs/slurm-%A_%a.out
#SBATCH --error=./logs/slurm-%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user="%u@asu.edu"

module load julia
ID=${SLURM_ARRAY_TASK_ID}
mkdir -p ../../data/sims/MS00_MWE/logs
srun julia launch.jl sims/MS00_MWE/array.csv ${ID} >> ../../data/sims/MS00_MWE/logs/launch.jl.${ID}.log
