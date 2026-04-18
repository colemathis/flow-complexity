#!/usr/bin/env python3
"""Plot assembly index vs reachable chemical space (log10 N) from GDB-17 results."""

import csv
import math
import os

# Attempt matplotlib; if unavailable, just produce the summary CSV
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("WARNING: matplotlib not available, will produce CSV only.")

RESULTS = "/home/achampa4/data/gdb17/samples/results.csv"
OUTDIR = "/home/achampa4/data/gdb17/samples"

# GDB-17 Table 5: reachable chemical space N(n) for each heavy atom count n
TABLE5 = {
    1:  3,
    2:  6,
    3:  14,
    4:  47,
    5:  219,
    6:  1091,
    7:  6029,
    8:  37435,
    9:  243233,
    10: 1670163,
    11: 12219460,
    12: 72051665,
    13: 836687200,
    14: 2921398415,
    15: 15084103347,
    16: 38033661355,
    17: 109481780580,
}

# Read results
with open(RESULTS, 'r') as f:
    reader = csv.DictReader(f)
    rows = list(reader)

# Build plot data
data = []
for row in rows:
    n = int(row['n_atoms'])
    a = int(row['assembly_index'])
    if row['timed_out'] == 'True' or a < 0:
        continue
    if n not in TABLE5:
        continue
    log_n = math.log10(TABLE5[n])
    data.append((n, log_n, a))

# Write enriched CSV
out_csv = os.path.join(OUTDIR, "assembly_vs_space.csv")
with open(out_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['n_atoms', 'log10_N', 'assembly_index'])
    for n, log_n, a in data:
        writer.writerow([n, round(log_n, 4), a])
print("Wrote", out_csv)

if not HAS_MPL:
    print("Install matplotlib to generate the plot.")
    exit(0)

# --- Plot ---
fig, ax = plt.subplots(figsize=(10, 6))

# Color by atom count
import matplotlib.cm as cm
n_vals = sorted(set(d[0] for d in data))
norm = plt.Normalize(min(n_vals), max(n_vals))
cmap = cm.viridis

for n, log_n, a in data:
    ax.scatter(log_n, a, c=[cmap(norm(n))], s=40, alpha=0.7, edgecolors='k', linewidth=0.3)

# Compute and plot mean assembly index per atom count
from collections import defaultdict
by_n = defaultdict(list)
for n, log_n, a in data:
    by_n[n].append(a)

means_x = []
means_y = []
for n in sorted(by_n):
    vals = by_n[n]
    mean_a = sum(vals) / len(vals)
    log_n = math.log10(TABLE5[n])
    means_x.append(log_n)
    means_y.append(mean_a)
    ax.annotate(f"n={n}", (log_n, mean_a), textcoords="offset points",
                xytext=(5, 5), fontsize=7, color='gray')

ax.plot(means_x, means_y, 'k--', alpha=0.5, linewidth=1, label='mean')

sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, pad=0.02)
cbar.set_label('Heavy atom count (n)')

ax.set_xlabel(r'Reachable chemical space $\log_{10}(N)$', fontsize=12)
ax.set_ylabel('Assembly index', fontsize=12)
ax.set_title('Assembly Index vs Reachable Chemical Space (GDB-17 sample)', fontsize=13)
ax.legend(loc='upper left')
ax.grid(True, alpha=0.3)

plt.tight_layout()
out_png = os.path.join(OUTDIR, "assembly_vs_space.png")
fig.savefig(out_png, dpi=150)
print("Wrote", out_png)
