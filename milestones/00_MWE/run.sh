#!/bin/bash
#SBATCH --job-name=flow-complexity
#SBATCH --partition=htc
#SBATCH --qos=public
#SBATCH --array=1-100
#SBATCH --time=01:00:00
#SBATCH --output=./data/logs/slurm_%a.log
#SBATCH --error=./data/logs/slurm_%a.log
#SBATCH --mail-type=NONE
#SBATCH --mail-user="%u@asu.edu"

module load julia
ID=${SLURM_ARRAY_TASK_ID}
mkdir -p ./data/logs
srun julia launch.jl ./data/array.csv ${ID} >> ./data/logs/launch.jl_${ID}.log
