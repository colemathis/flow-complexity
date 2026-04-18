#!/usr/bin/env python3
"""Compute assembly indices for all molecules in the manifest using the CLI binary."""

import subprocess
import csv
import os
import sys
import time

BINARY = os.path.expanduser("~/repos/assembly-theory/target/release/assembly-theory")
SAMPLES_DIR = "/home/achampa4/data/gdb17/samples"
MANIFEST = os.path.join(SAMPLES_DIR, "manifest.csv")
OUTPUT = os.path.join(SAMPLES_DIR, "results.csv")
TIMEOUT_MS = 120000  # 2 minutes per molecule
TIMEOUT_S = TIMEOUT_MS / 1000 + 5  # subprocess timeout with buffer

with open(MANIFEST, 'r') as f:
    reader = csv.DictReader(f)
    rows = list(reader)

print(f"Computing assembly indices for {len(rows)} molecules...")
print(f"Binary: {BINARY}")
print(f"Timeout: {TIMEOUT_MS}ms per molecule")

results = []
for i, row in enumerate(rows):
    mol_path = os.path.join(SAMPLES_DIR, row['mol_file'])
    smiles = row['smiles']
    n_atoms = row['n_atoms']

    t0 = time.time()
    try:
        proc = subprocess.run(
            [BINARY, "--timeout", str(TIMEOUT_MS), mol_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=TIMEOUT_S
        )
        raw = proc.stdout.strip()
        elapsed = time.time() - t0

        # Parse output: could be "6" or "<= 8 (timed out)"
        if "timed out" in raw:
            idx_str = raw.split()[1] if len(raw.split()) > 1 else raw
            timed_out = True
        else:
            idx_str = raw
            timed_out = False

        try:
            assembly_idx = int(idx_str)
        except ValueError:
            assembly_idx = -1
            timed_out = True

    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        assembly_idx = -1
        timed_out = True
    except Exception as e:
        elapsed = time.time() - t0
        assembly_idx = -1
        timed_out = True
        print(f"  ERROR on {row['mol_file']}: {e}")

    results.append({
        'mol_file': row['mol_file'],
        'smiles': smiles,
        'n_atoms': n_atoms,
        'assembly_index': assembly_idx,
        'timed_out': timed_out,
        'time_s': round(elapsed, 2)
    })

    status = "TIMEOUT" if timed_out else f"a={assembly_idx}"
    print(f"  [{i+1}/{len(rows)}] n={n_atoms:>2s} {status:>10s}  {elapsed:6.1f}s  {smiles[:50]}")

with open(OUTPUT, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['mol_file', 'smiles', 'n_atoms', 'assembly_index', 'timed_out', 'time_s'])
    writer.writeheader()
    writer.writerows(results)

n_ok = sum(1 for r in results if not r['timed_out'])
n_to = sum(1 for r in results if r['timed_out'])
print(f"\nDone. {n_ok} succeeded, {n_to} timed out.")
print(f"Results: {OUTPUT}")
