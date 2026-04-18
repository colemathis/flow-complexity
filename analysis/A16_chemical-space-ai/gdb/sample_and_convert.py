#!/usr/bin/env python3
"""Sample molecules from GDB-17 stratified by heavy atom count, convert to .mol files."""

import os
import re
import random
import sys
import csv

from rdkit import Chem

INPUT = "/home/achampa4/data/gdb17/GDB17.50000000.smi"
OUTPUT_DIR = "/home/achampa4/data/gdb17/samples"
MANIFEST = os.path.join(OUTPUT_DIR, "manifest.csv")

# How many molecules to sample per atom count bucket
# Usage: python3 sample_and_convert.py [sample_per_bucket]
SAMPLE_PER_BUCKET = int(sys.argv[1]) if len(sys.argv) > 1 else 1000

# Threshold: take ALL molecules if bucket has fewer than this many
ALL_THRESHOLD = SAMPLE_PER_BUCKET * 2

SEED = 42

# --- SMILES atom counting (same logic as count_atoms.py) ---
TARGET = frozenset({'C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I', 'c', 'n', 'o', 's'})
ATOM_RE = re.compile(r'Cl|Br|\[([^\]]+)\]|([A-Za-z])')

def count_atoms(smiles):
    count = 0
    for m in ATOM_RE.finditer(smiles):
        if m.group(1):
            inner = m.group(1)
            j = 0
            while j < len(inner) and inner[j].isdigit():
                j += 1
            elem = ''
            if j < len(inner) and inner[j].isalpha():
                elem = inner[j]
                j += 1
                if j < len(inner) and inner[j].islower():
                    elem += inner[j]
            if elem in TARGET:
                count += 1
        else:
            token = m.group(0)
            if token in TARGET:
                count += 1
    return count

# --- Phase 1: Read and bucket all SMILES by atom count ---
print(f"Reading {INPUT} ...")
print(f"Sample per bucket: {SAMPLE_PER_BUCKET}")

buckets = {}  # n_atoms -> list of SMILES
total = 0
with open(INPUT, 'r') as f:
    for line in f:
        smi = line.strip()
        if not smi:
            continue
        smi = smi.split()[0]
        n = count_atoms(smi)
        if n < 1 or n > 17:
            continue
        if n not in buckets:
            buckets[n] = []
        buckets[n].append(smi)
        total += 1
        if total % 5_000_000 == 0:
            print(f"  read {total:,} lines...")

print(f"Total valid molecules: {total:,}")
for n in sorted(buckets):
    print(f"  n={n:2d}: {len(buckets[n]):>12,}")

# --- Phase 2: Sample from each bucket ---
random.seed(SEED)
sampled = {}  # n_atoms -> list of SMILES
for n in sorted(buckets):
    pool = buckets[n]
    if len(pool) <= ALL_THRESHOLD:
        sampled[n] = pool
    else:
        sampled[n] = random.sample(pool, SAMPLE_PER_BUCKET)
    print(f"  n={n:2d}: sampled {len(sampled[n]):,} / {len(pool):,}")

# Free memory
del buckets

# --- Phase 3: Convert to .mol files and write manifest ---
os.makedirs(OUTPUT_DIR, exist_ok=True)

total_written = 0
total_failed = 0

with open(MANIFEST, 'w', newline='') as mf:
    writer = csv.writer(mf)
    writer.writerow(['mol_file', 'smiles', 'n_atoms'])

    for n in sorted(sampled):
        subdir = os.path.join(OUTPUT_DIR, f"n{n:02d}")
        os.makedirs(subdir, exist_ok=True)

        for idx, smi in enumerate(sampled[n], 1):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                total_failed += 1
                continue
            mol_block = Chem.MolToMolBlock(mol)
            fname = f"{idx:05d}.mol"
            fpath = os.path.join(subdir, fname)
            with open(fpath, 'w') as f:
                f.write(mol_block)
            rel_path = os.path.join(f"n{n:02d}", fname)
            writer.writerow([rel_path, smi, n])
            total_written += 1

print(f"\nDone. {total_written:,} .mol files written, {total_failed} SMILES failed to parse.")
print(f"Manifest: {MANIFEST}")
