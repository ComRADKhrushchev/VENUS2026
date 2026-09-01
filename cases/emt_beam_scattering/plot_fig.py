#!/usr/bin/env python
"""emt_beam_scattering: C z(t) and r_min(t) for NT=20 EMT-NN trajectories.

Data: fort.1001..fort.1020, one file per trajectory (GWRITE_LEVEL=1). Each
`--- step` block header carries E0/T/H(eV), Q(C), V(C) and Au_nearest/r_min.
Trajectories may terminate early (scattered: Q(3)>=RMAX & receding) or run
the full NS=30000 cycles (=3000 fs, adsorbed). Panels (ps axis):
  (a) C height z(t), 20 trajectories overlaid
  (b) C-Au nearest distance r_min(t)
"""
import math
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NTRJ = 20


def parse(path):
    ts, zs, rmins = [], [], []
    txt = open(path, encoding='utf-8', errors='ignore').read()
    for blk in txt.split('--- step')[1:]:
        m = re.search(r't\(fs\)=\s*(-?[\d.]+)', blk)
        q = re.search(r'Q\(C\)=\s*(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)', blk)
        r = re.search(r'r_min\s+([\d.]+)', blk)
        if not (m and q):
            continue
        z = float(q.group(3))
        if not math.isfinite(z):
            z = np.nan
        ts.append(float(m.group(1)))
        zs.append(z)
        rmins.append(float(r.group(1)) if r else np.nan)
    return np.array(ts), np.array(zs), np.array(rmins)


cmap = plt.cm.viridis
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 9), sharex=True)

stats = []
for i in range(NTRJ):
    ts, zs, rmins = parse(f'fort.{1001 + i}')
    if len(ts) == 0:
        continue
    tps = ts / 1000.0  # fs -> ps
    c = cmap(i / max(NTRJ - 1, 1))
    lw = 0.9 if NTRJ > 10 else 1.2
    ax1.plot(tps, zs, color=c, lw=lw)
    ax2.plot(tps, rmins, color=c, lw=lw)
    full = len(ts) >= 600  # NIP=50 over NS=30000 -> 601 blocks if full run
    stats.append((i + 1, len(ts), zs.min(), full))

ax1.set_ylabel('C height z (Å)')
ax1.set_title(f'EMT-NN C/Au(111) beam scattering, E$_{{rel}}$=14.528 kcal/mol, NT={NTRJ}')
ax2.set_ylabel('C–Au nearest r$_{min}$ (Å)')
ax2.set_xlabel('Time (ps)')
for ax in (ax1, ax2):
    ax.set_xlim(0, 3.05)
    ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig('fig_scattering.png', dpi=150)

full = sum(1 for s in stats if s[3])
early = [s for s in stats if not s[3]]
zmins = [s[2] for s in stats]
print(f'trajectories parsed: {len(stats)}')
print(f'full-NS (adsorbed) runs: {full}; early-terminated: {len(early)}')
for s in early:
    print(f'  trj {s[0]}: {s[1]} blocks, z_min={s[2]:.3f}')
print(f'z_min over all: {min(zmins):.3f} A; mean z_min={np.mean(zmins):.3f} A')
