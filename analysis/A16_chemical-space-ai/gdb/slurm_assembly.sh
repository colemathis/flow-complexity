#!/bin/bash
#SBATCH --job-name=assembly-gdb17
#SBATCH --partition=htc
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=1-200
#SBATCH --output=/home/achampa4/data/gdb17/samples/logs/slurm_%A_%a.out
#SBATCH --error=/home/achampa4/data/gdb17/samples/logs/slurm_%A_%a.err

# Assembly index computation for GDB-17 molecules — SLURM array job
# Each array task processes one chunk of the manifest.

BINARY="$HOME/repos/assembly-theory/target/release/assembly-theory"
SAMPLES_DIR="$HOME/data/gdb17/samples"
CHUNKS_DIR="$SAMPLES_DIR/chunks"
RESULTS_DIR="$SAMPLES_DIR/results"
TIMEOUT_MS=120000  # 2 minutes per molecule

module load rust 2>/dev/null
module load gcc/15.2.0 2>/dev/null

TASK_ID=$(printf '%04d' $SLURM_ARRAY_TASK_ID)

python3 "$HOME/data/gdb17/worker.py" \
    "$BINARY" \
    "$SAMPLES_DIR/chunks/chunk_${TASK_ID}.csv" \
    "$SAMPLES_DIR/results/batch_${TASK_ID}.csv" \
    "$SAMPLES_DIR" \
    "$TIMEOUT_MS"
