#!/usr/bin/env python3
"""Worker script: compute assembly indices for one chunk of molecules."""

import subprocess
import csv
import os
import sys
import time

BINARY = sys.argv[1]
CHUNK_FILE = sys.argv[2]
OUTPUT_FILE = sys.argv[3]
SAMPLES_DIR = sys.argv[4]
TIMEOUT_MS = int(sys.argv[5])
TIMEOUT_S = TIMEOUT_MS / 1000 + 5

if not os.path.isfile(CHUNK_FILE):
    print("No chunk file: %s" % CHUNK_FILE)
    sys.exit(0)

with open(CHUNK_FILE, 'r') as f:
    reader = csv.DictReader(f)
    rows = list(reader)

results = []
for row in rows:
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

        if "timed out" in raw:
            parts = raw.split()
            idx_str = parts[1] if len(parts) > 1 else "-1"
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

    results.append({
        'mol_file': row['mol_file'],
        'smiles': smiles,
        'n_atoms': n_atoms,
        'assembly_index': assembly_idx,
        'timed_out': timed_out,
        'time_s': round(elapsed, 2)
    })

with open(OUTPUT_FILE, 'w') as f:
    writer = csv.DictWriter(f, fieldnames=['mol_file', 'smiles', 'n_atoms', 'assembly_index', 'timed_out', 'time_s'])
    writer.writeheader()
    writer.writerows(results)

n_ok = sum(1 for r in results if not r['timed_out'])
print("Done: %d molecules, %d succeeded, %d timed out" % (len(results), n_ok, len(results) - n_ok))
