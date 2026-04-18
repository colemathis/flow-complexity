#!/bin/bash
# Full pipeline launcher: sample, chunk, submit SLURM array job
set -e

SAMPLES_DIR="$HOME/data/gdb17/samples"
CHUNKS_DIR="$SAMPLES_DIR/chunks"
RESULTS_DIR="$SAMPLES_DIR/results"
LOGS_DIR="$SAMPLES_DIR/logs"
SAMPLE_PER_BUCKET="${1:-1000}"
CHUNK_SIZE="${2:-50}"

echo "=== Step 1: Sample & convert (${SAMPLE_PER_BUCKET} per bucket) ==="
cd "$HOME/data/gdb17"
python3 sample_and_convert.py "$SAMPLE_PER_BUCKET"

echo ""
echo "=== Step 2: Split manifest into chunks of ${CHUNK_SIZE} ==="
mkdir -p "$CHUNKS_DIR" "$RESULTS_DIR" "$LOGS_DIR"
# Remove old chunks
rm -f "$CHUNKS_DIR"/chunk_*.csv

# Get header
HEADER=$(head -1 "$SAMPLES_DIR/manifest.csv")
TOTAL=$(tail -n +2 "$SAMPLES_DIR/manifest.csv" | wc -l)
NUM_CHUNKS=$(( (TOTAL + CHUNK_SIZE - 1) / CHUNK_SIZE ))

echo "Total molecules: $TOTAL"
echo "Chunk size: $CHUNK_SIZE"
echo "Number of chunks: $NUM_CHUNKS"

# Split (skip header, then chunk)
CHUNK_IDX=1
LINE_IDX=0
CURRENT_CHUNK=""

tail -n +2 "$SAMPLES_DIR/manifest.csv" | while IFS= read -r line; do
    if [ $((LINE_IDX % CHUNK_SIZE)) -eq 0 ]; then
        CURRENT_CHUNK="$CHUNKS_DIR/chunk_$(printf '%04d' $CHUNK_IDX).csv"
        echo "$HEADER" > "$CURRENT_CHUNK"
        CHUNK_IDX=$((CHUNK_IDX + 1))
    fi
    echo "$line" >> "$CURRENT_CHUNK"
    LINE_IDX=$((LINE_IDX + 1))
done

echo "Created $(ls "$CHUNKS_DIR"/chunk_*.csv 2>/dev/null | wc -l) chunk files."

echo ""
echo "=== Step 3: Submit SLURM array job ==="
# Update array size in the SLURM script
SLURM_SCRIPT="$HOME/data/gdb17/slurm_assembly.sh"

# Submit with correct array range
JOB_ID=$(sbatch --array=1-${NUM_CHUNKS} "$SLURM_SCRIPT" | awk '{print $4}')
echo "Submitted job: $JOB_ID (array 1-${NUM_CHUNKS})"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  sacct -j $JOB_ID --format=JobID,State,Elapsed,MaxRSS"
echo ""
echo "When done, collect results with:"
echo "  cat $RESULTS_DIR/batch_*.csv | head -1 > $SAMPLES_DIR/results.csv"
echo "  cat $RESULTS_DIR/batch_*.csv | grep -v '^mol_file' >> $SAMPLES_DIR/results.csv"
echo "  python3 $HOME/data/gdb17/plot_assembly.py"
