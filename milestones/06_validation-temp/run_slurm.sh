#!/bin/bash

#SBATCH --job-name=flow-complexity
#SBATCH --partition=htc
#SBATCH --qos=public
#SBATCH --array=1-10
#SBATCH --time=04:00:00
#SBATCH --output=./data/logs/slurm_%a.log
#SBATCH --error=./data/logs/slurm_%a.log
#SBATCH --mail-type=NONE
#SBATCH --mail-user="%u@asu.edu"
#SBATCH --mem=4GB

# Load necessary modules
module load julia

# Create logs directory if it doesn't exist
mkdir -p ./data/logs

# Get the task ID from the SLURM array
ID=${SLURM_ARRAY_TASK_ID}

# Run the Julia script with the task ID and redirect output to a log file
srun julia launch.jl ./data/array.csv ${ID} >> ./data/logs/launch.jl_${ID}.log
