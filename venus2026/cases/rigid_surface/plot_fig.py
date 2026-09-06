#!/usr/bin/env python
"""rigid_surface (GLO surface oscillator): H z(t) + surface-Au oscillation.

Data: fort.1001..fort.1020, NS=20000 (2 ps). The surface atom (Au) is
coupled to a dissipative ghost (WS1/WG1 = 52.1 cm^-1, FCG friction) and
is thermally sampled by GLOSELECT at TVIB_B; the H atom starts 8 A above
the surface. All trajectories finite; the Au oscillator is excited by the
impact. The H passes through z=0 to the mirror position - an intrinsic
property of the HARMONIC test PES (no close-range repulsion); real
surface scattering should use the RST PES.
"""
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse(path):
    ts, z_h, z_au = [], [], []
    t = zh = za = None
    with open(path, encoding='utf-8', errors='ignore') as fh:
        for blk in fh.read().split('--- step')[1:]:
            m = re.search(r't\(fs\)=\s*(-?[\d.]+)', blk)
            got_h = got_a = False
            for ln in blk.splitlines():
                p = ln.split()
                if len(p) == 7 and not got_h and p[0] == 'C':
                    try:
                        zh = float(p[3]); got_h = True
                    except ValueError:
                        pass
                elif len(p) == 7 and not got_a and p[0] == 'Au':
                    try:
                        za = float(p[3]); got_a = True
                    except ValueError:
                        pass
                if got_h and got_a:
                    break
            if m:
                ts.append(float(m.group(1)))
                z_h.append(zh)
                z_au.append(za)
    return np.array(ts), np.array(z_h), np.array(z_au)


data = [parse(f'fort.1{k:03d}') for k in range(1, 21)]
n_bad = sum(int((~np.isfinite(np.array(zh)[np.array(zh) is not None])).sum())
            for _, zh, _ in data if zh is not None and len(zh))
print(f'20 trajectories, non-finite: {n_bad}')

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
cmap = plt.cm.viridis
for i, (ts, zh, za) in enumerate(data):
    c = cmap(i / 19)
    ax1.plot(ts / 1000.0, zh, lw=0.8, color=c)
    ax2.plot(ts / 1000.0, za, lw=0.8, color=c)

ax1.axhline(0.0, color='crimson', lw=1.0, ls='--', alpha=0.8)
ax1.annotate('surface plane (z=0)', xy=(0.02, 0.06),
             xycoords='axes fraction', fontsize=9, color='crimson')
ax1.set_ylabel('H z [$\\AA$]')
ax1.set_title('rigid_surface GLO surface oscillator, NT=20\n'
              'H atom (top) + surface Au with ghost spring (bottom)', fontsize=11)
ax1.grid(alpha=0.25)

ax2.set_xlabel('time t [ps]')
ax2.set_ylabel('surface Au z [$\\AA$]')
ax2.grid(alpha=0.25)
ax2.annotate('oscillator excited by impact (WS1=52 cm$^{-1}$)',
             xy=(0.02, 0.85), xycoords='axes fraction', fontsize=9, color='gray')

fig.text(0.99, 0.01,
         'All 20 trajectories finite. H passes z=0 to the mirror position: '
         'intrinsic to the HARMONIC test PES (no close-range repulsion); '
         'real surface scattering should use the RST PES.',
         ha='right', fontsize=7, color='#2E7D32')
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig('fig_long_trajectory.png', dpi=150)
print('saved fig_long_trajectory.png')
