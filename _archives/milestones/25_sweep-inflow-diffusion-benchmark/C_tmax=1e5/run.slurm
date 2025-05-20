#!/bin/bash

# Job identification
#SBATCH --job-name=MS25C              # descriptive job name
#SBATCH --array=1-100                    # simulation indices

# Logging (stdout & stderr go to the same file)
#SBATCH --output=./data/logs/slurm_%A_%a.log   # %A = array master ID

# Resources
#SBATCH --partition=general              # 7‑day queue
#SBATCH --qos=public                     # default QoS
#SBATCH --time=72:00:00                  # wall time (hh:mm:ss)
#SBATCH --mem=4G                         # memory

#################################
# SCRIPT BEGINS
#################################

# The following line allows launching the job from another job (e.g., vscode)
unset $(printenv | grep SLURM_CPU_BIND | sed -E 's/(.*)=.*/\1/' | xargs)

module load julia
mkdir -p ./data/logs
srun flow launch ${SLURM_ARRAY_TASK_ID}
