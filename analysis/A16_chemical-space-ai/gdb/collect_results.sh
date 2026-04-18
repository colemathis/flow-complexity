#!/bin/bash
# Collect SLURM batch results and regenerate the plot
set -e

SAMPLES_DIR="$HOME/data/gdb17/samples"
RESULTS_DIR="$SAMPLES_DIR/results"
FINAL="$SAMPLES_DIR/results.csv"

echo "Collecting results from $RESULTS_DIR ..."

# Merge all batch CSVs (take header from first, skip headers from rest)
head -1 "$(ls "$RESULTS_DIR"/batch_*.csv | head -1)" > "$FINAL"
for f in "$RESULTS_DIR"/batch_*.csv; do
    tail -n +2 "$f" >> "$FINAL"
done

TOTAL=$(tail -n +2 "$FINAL" | wc -l)
OK=$(grep -c "False" "$FINAL" || true)
TIMEOUT=$(grep -c "True" "$FINAL" || true)

echo "Total: $TOTAL molecules"
echo "Succeeded: $OK"
echo "Timed out: $TIMEOUT"
echo "Results: $FINAL"

echo ""
echo "Generating plot..."
cd "$HOME/data/gdb17"
python3 plot_assembly.py

echo "Done!"
