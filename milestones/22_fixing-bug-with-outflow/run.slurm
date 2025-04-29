#!/bin/bash

# Job identification
#SBATCH --job-name=flow-sim              # descriptive job name
#SBATCH --array=1-100                    # simulation indices

# Logging (stdout & stderr go to the same file)
#SBATCH --output=./data/logs/slurm_%A_%a.log   # %A = array master ID

# Resources
#SBATCH --partition=general              # 7‑day queue
#SBATCH --qos=public                     # default QoS
#SBATCH --time=24:00:00                  # wall time (hh:mm:ss)
#SBATCH --mem=4G                         # memory

#################################
# SCRIPT BEGINS
#################################

module load julia
mkdir -p ./data/logs
srun flow launch ${SLURM_ARRAY_TASK_ID}
