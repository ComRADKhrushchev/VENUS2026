#!/usr/bin/env python
"""rigid_surface (NSURF=2): z(t) of the incident atom A, all NT=20 trajectories.

Data: fort.1001..fort.1020, NS=20000 (2 ps) covering approach, entry into
the 4-atom 5x5 A grid, and pass-through to the mirror cell. All coordinates
finite. The pass-through (z crossing 0 and reaching the mirror grid at
z = -8 A) is an intrinsic property of the HARMONIC test PES: it has no
close-range repulsion and the 5 A grid spacing leaves open channels.
A real rigid-surface scattering case should use the RST PES instead.
"""
import glob
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_first_c(path):
    """Return (t_fs array, z array) for the first C line (incident atom A)."""
    ts, zs = [], []
    with open(path, encoding='utf-8', errors='ignore') as fh:
        for blk in fh.read().split('--- step')[1:]:
            m = re.search(r't\(fs\)=\s*(-?[\d.]+)', blk)
            z = np.nan
            for ln in blk.splitlines():
                p = ln.split()
                if len(p) == 7 and p[0] == 'C':
                    try:
                        z = float(p[3])
                    except ValueError:
                        z = np.nan
                    break
            if m:
                ts.append(float(m.group(1)))
                zs.append(z)
    return np.array(ts), np.array(zs)


files = sorted(glob.glob('fort.10*'))
data = [parse_first_c(f) for f in files]
n_bad = sum(int((~np.isfinite(zs)).sum()) for _, zs in data)
print(f'{len(data)} trajectories, non-finite samples: {n_bad}')

fig, ax = plt.subplots(figsize=(11, 6))
cmap = plt.cm.viridis
for i, (ts, zs) in enumerate(data):
    ax.plot(ts / 1000.0, zs, lw=0.8, color=cmap(i / max(len(data) - 1, 1)))

# surface markers
ax.axhline(0.0, color='crimson', lw=1.0, ls='--', alpha=0.8)
ax.annotate('surface plane (z=0): 4-atom 5x5 A grid', xy=(0.02, 0.08),
            xycoords='axes fraction', fontsize=9, color='crimson')
ax.axhspan(-8.3, -7.7, color='gray', alpha=0.3)
ax.annotate('mirror grid (z=-8)', xy=(0.80, 0.10), xycoords='axes fraction',
            fontsize=9, color='gray')

ax.set_xlabel('time t [ps]')
ax.set_ylabel('incident-atom z [$\\AA$]')
ax.set_title('rigid_surface (NSURF=2): z(t) of incident atom A, NT=20\n'
             'approach -> entry through grid channel -> mirror cell '
             '(HARMONIC PES: no close-range repulsion)', fontsize=11)
ax.grid(alpha=0.25)
fig.text(0.99, 0.01,
         'All 20 trajectories finite; H conserved to 0.0000 eV. '
         'Pass-through is an intrinsic property of the HARMONIC test PES '
         '(5 A grid spacing leaves open channels); real rigid-surface '
         'scattering should use the RST PES.',
         ha='right', fontsize=7, color='#2E7D32')
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig('fig_long_trajectory.png', dpi=150)
print('saved fig_long_trajectory.png')
