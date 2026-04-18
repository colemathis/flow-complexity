#!/usr/bin/env python3
"""Count heavy atoms (C,N,O,S,F,Cl,Br,I) per molecule in a SMILES file."""

import re
import sys

INPUT = "/home/achampa4/data/gdb17/GDB17.50000000.smi"
OUTPUT = "/home/achampa4/data/gdb17/atom_count_stats.txt"

TARGET = frozenset({'C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I',
                     'c', 'n', 'o', 's'})

# Tokenize SMILES: try two-letter atoms first, then bracket atoms, then single letters
ATOM_RE = re.compile(r'Cl|Br|\[([^\]]+)\]|([A-Za-z])')


def count_atoms(smiles):
    count = 0
    for m in ATOM_RE.finditer(smiles):
        if m.group(1):  # bracket atom [...]
            inner = m.group(1)
            j = 0
            # skip isotope digits
            while j < len(inner) and inner[j].isdigit():
                j += 1
            # extract element symbol (uppercase + optional lowercase)
            elem = ''
            if j < len(inner) and inner[j].isalpha():
                elem = inner[j]
                j += 1
                if j < len(inner) and inner[j].islower():
                    elem += inner[j]
            if elem in TARGET:
                count += 1
        else:  # bare atom: single letter or Cl/Br
            token = m.group(0)
            if token in TARGET:
                count += 1
    return count


histogram = {}
total = 0

with open(INPUT, 'r') as f:
    for line in f:
        smiles = line.strip()
        if not smiles:
            continue
        smiles = smiles.split()[0]
        n = count_atoms(smiles)
        histogram[n] = histogram.get(n, 0) + 1
        total += 1
        if total % 5_000_000 == 0:
            print(f"  processed {total:,} lines...", flush=True)

with open(OUTPUT, 'w') as f:
    f.write(f"GDB-17 Sample: Heavy Atom Count Distribution\n")
    f.write(f"Atoms counted: C, N, O, S, F, Cl, Br, I\n")
    f.write(f"Total molecules: {total:,}\n")
    f.write(f"{'=' * 50}\n")
    f.write(f"{'Atoms':>6}  {'Count':>12}  {'Percent':>8}\n")
    f.write(f"{'-' * 50}\n")
    for n in sorted(histogram):
        pct = 100.0 * histogram[n] / total
        f.write(f"{n:>6}  {histogram[n]:>12,}  {pct:>7.3f}%\n")

print(f"\nDone. {total:,} molecules processed.")
print(f"Statistics saved to {OUTPUT}")
