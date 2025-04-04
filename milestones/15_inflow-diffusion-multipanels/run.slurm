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

module load julia

LOGS_DIR=./data/logs
mkdir -p ${LOGS_DIR}

ID=${SLURM_ARRAY_TASK_ID}
ID_PADDED=$(printf "%06d" $ID)
LOG_FN=$(printf flow_${ID_PADDED}.log)
srun flow launch ${ID} >> ${LOGS_DIR}/${LOG_FN}
