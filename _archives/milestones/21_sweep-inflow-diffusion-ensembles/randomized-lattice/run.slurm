#!/bin/bash

# Parameters which shouldn’t be modified
#SBATCH --output=./data/logs/slurm_%a.log
#SBATCH --error=./data/logs/slurm_%a.log
#SBATCH --mail-type=NONE
#SBATCH --mail-user="%u@asu.edu"
#SBATCH --qos=public

#################################
# CUSTOMIZABLE SLURM PARAMETERS
#################################

# Job name:
#
#   Sets the job name

#SBATCH --job-name=MS21-random-lattice

# Partition:
#
#   Determines wall time: htc (4 hours, default), general (7 days)

#SBATCH --partition=general

# Array:
#
#   Determines the job array (i.e., simulation numbers)

#SBATCH --array=1-1000

# Time:
#
#   Determines wall time, memory determines max available memory

#SBATCH --time=48:00:00
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

# Set logfile
# ID_PADDED=$(printf "%06d" $ID)
# LOG_FN=$(printf flow_${ID_PADDED}.log)
# Note: probably obsolete, disabled second log file. Also removed ref to it below.
# This is captured by SLURM’s logfiles anyway.

# Launch simulation and redirect output to logfile
# srun flow launch ${ID} >> ${LOGS_DIR}/${LOG_FN}
srun flow launch ${ID}
