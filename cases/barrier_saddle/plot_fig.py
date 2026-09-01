#!/usr/bin/env python
"""barrier_saddle: long-trajectory r_min(t) for 2 saddle-crossing trajectories.

Data: fort.1001, fort.1002 (3 ps each, 30001 records of 3 atoms).
r_min(t) = distance of atom 1 to its nearest other atom.
Time axis: t(fs)/1000 -> ps.
"""
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

fig, ax = plt.subplots(figsize=(10, 5.5))
for f, c in [('fort.1001', '#4878CF'), ('fort.1002', '#EE854A')]:
    recs = parse(f)
    if len(recs) > 1 and recs[0][0] > recs[1][0]:
        recs = recs[1:]
    ts = np.array([t for t, _ in recs]) / 1000.0
    rmin = np.array([min(np.linalg.norm(x[0] - x[1]), np.linalg.norm(x[0] - x[2]))
                     for _, x in recs])
    ax.plot(ts, rmin, lw=0.9, color=c, alpha=0.9,
            label=f'{f}: min {rmin.min():.3f} A at t={ts[rmin.argmin()]:.3f} ps')
    print(f'{f}: n={len(ts)} t=[{ts[0]:.3f},{ts[-1]:.3f}] ps r_min range=[{rmin.min():.4f},{rmin.max():.4f}]')

ax.set_xlabel('time t [ps]')
ax.set_ylabel('r$_{min}$: atom 1 to nearest atom [$\\AA$]')
ax.set_title('Barrier/saddle trajectories (LEPS H$_3$, 2 trajectories, 3 ps):\n'
             'r$_{min}$(t) - saddle crossing when r$_{min}$ dips', fontsize=12)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig('fig_long_trajectory.png', dpi=150)
print('saved fig_long_trajectory.png')
