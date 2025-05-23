#!/bin/bash

# Parameters which shouldn’t be modified
#SBATCH --job-name=MS18
#SBATCH --output=./data/logs/slurm_%a.log
#SBATCH --error=./data/logs/slurm_%a.log
#SBATCH --mail-type=NONE
#SBATCH --mail-user="%u@asu.edu"
#SBATCH --qos=public

#################################
# CUSTOMIZABLE SLURM PARAMETERS
#################################

# Partition:
#
#   Determines wall time: htc (4 hours, default), general (7 days)

#SBATCH --partition=general

# Array:
#
#   Determines the job array (i.e., simulation numbers)

#SBATCH --array=1-100

# Time:
#
#   Determines wall time, memory determines max available memory

#SBATCH --time=24:00:00
#SBATCH --mem=4GB

#################################
# SCRIPT BEGINS
#################################

# Load required modules
module load julia

# Set log directory and create it
LOGS_DIR=./data/logs
mkdir -p ${LOGS_DIR}

# Get current array id
ID=${SLURM_ARRAY_TASK_ID}

# Launch simulation and redirect output to logfile
srun flow launch ${ID}
