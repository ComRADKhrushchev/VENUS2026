#!/usr/bin/env python
"""twobody_collision: long-trajectory r_AB(t) overlay (8 trajectories).

Data: fort.1001..fort.1200 (200 trajectories, per-record blocks with header
'--- step N t(fs)=... ---' and C x y z vx vy vz lines for the 3 H atoms).
A = H2 (atoms 1+2), B = atom 3. r_AB = |R_A(com of 1,2) - r_3|.
Note: each per-trajectory file begins with one summary record of the final
state before the history replay; that record is dropped (dedup by step).
Time axis: t(fs)/1000 -> ps.
"""
import os
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse(path):
    recs = []
    with open(path, encoding='utf-8', errors='ignore') as fh:
        for blk in fh.read().split('--- step')[1:]:
            m = re.search(r't\(fs\)=\s*([-\d.]+)', blk)
            coords = []
            for ln in blk.splitlines():
                p = ln.split()
                if len(p) == 7 and p[0] == 'C':
                    try:
                        coords.append([float(x) for x in p[1:4]])
                    except ValueError:
                        pass
            if m and len(coords) >= 3:
                recs.append((float(m.group(1)), np.array(coords[:3])))
    return recs

files = sorted([f for f in os.listdir('.') if re.fullmatch(r'fort\.(1[0-9]{3})', f)],
               key=lambda f: int(f.split('.')[1]))

def load(f):
    recs = parse(f)
    if len(recs) > 1 and recs[0][0] > recs[1][0]:   # drop leading final-state record
        recs = recs[1:]
    ts = np.array([t for t, _ in recs]) / 1000.0
    rab = np.array([np.linalg.norm((xyz[0] + xyz[1]) / 2 - xyz[2]) for _, xyz in recs])
    return ts, rab

# pick the trajectories that actually reach small r_AB (real close encounters)
scored = []
for f in files:
    ts, rab = load(f)
    if len(ts) >= 50:
        scored.append((rab.min(), f, ts, rab))
scored.sort(key=lambda s: s[0])
# prefer longer records among the closest encounters
close = [s for s in scored if s[0] < 1.5 and len(s[2]) > 2000][:8]
sel = close if len(close) >= 5 else scored[:8]
print('selected:', [(s[1], round(s[0], 2), len(s[2])) for s in sel])

fig, ax = plt.subplots(figsize=(10, 5.5))
cmap = plt.cm.viridis(np.linspace(0, 0.9, len(sel)))
for (_, f, ts, rab), c in zip(sel, cmap):
    ax.plot(ts, rab, lw=1.2, color=c, label=f'{f} (min {rab.min():.2f} A)')

ax.set_ylim(0, 15)
ax.set_xlabel('time t [ps]')
ax.set_ylabel(r'$r_{AB}$ = |R(H$_2$ com) - H$_B$| [$\AA$]')
ax.set_title('Two-body collision H$_2$ + H (LEPS batch, 8 closest encounters of 200 trajectories):\n'
             'approach - collision (r$_{AB}$ dip) - scattering (axis clipped at 15 A)', fontsize=12)
ax.legend(fontsize=8, ncol=2, loc='upper right')
fig.tight_layout()
fig.savefig('fig_long_trajectory.png', dpi=150)
print('saved fig_long_trajectory.png')
for _, f, ts, rab in sel:
    print(f'{f}: n={len(ts)} t=[{ts[0]:.3f},{ts[-1]:.3f}] ps r_AB min={rab.min():.3f} A')
