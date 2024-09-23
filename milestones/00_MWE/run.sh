#!/bin/bash
#SBATCH --job-name=MS00
#SBATCH --partition=htc
#SBATCH --qos=public
#SBATCH --array=1-1000
#SBATCH --time=01:00:00
#SBATCH --output=./data/logs/slurm-%A_%a.out
#SBATCH --error=./data/logs/slurm-%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user="%u@asu.edu"

module load julia
ID=${SLURM_ARRAY_TASK_ID}
mkdir -p ./data/logs
srun julia launch.jl ./data/array.csv ${ID} >> ./data/logs/launch.jl.${ID}.log
