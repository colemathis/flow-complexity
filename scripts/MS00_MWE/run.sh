#!/bin/bash
#SBATCH -t 0-1
#SBATCH --job-name=test-array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000
#SBATCH -e /dev/null
#SBATCH -o /dev/null

module load julia

echo "This is array task ${SLURM_ARRAY_TASK_ID}" >> test-array.sh.log

julia test-array.jl ${SLURM_ARRAY_TASK_ID} >> test-array.jl.log
