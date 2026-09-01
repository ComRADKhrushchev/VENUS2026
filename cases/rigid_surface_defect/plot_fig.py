#!/usr/bin/env python
"""rigid_surface_defect: DEFECT EVIDENCE - coordinates are Infinity/NaN from step 0.

All 20 trajectory files (fort.1001..fort.1020) emit `Infinity` (step 0) and
`NaN` (step >= 1) in Q(C)/C lines from the very first integration step, so
z(t) and r_min(t) are undefined. This plot records that failure: x(t) of the
incident atom is shown where finite, with NaN/Inf samples marked by vertical
markers - the figure is the defect evidence, not a physical trajectory.
"""
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_first_c(path):
    """Return (t_fs array, x array with NaN where non-finite) for first C line."""
    ts, xs = [], []
    with open(path, encoding='utf-8', errors='ignore') as fh:
        for blk in fh.read().split('--- step')[1:]:
            m = re.search(r't\(fs\)=\s*(-?[\d.]+)', blk)
            x = np.nan
            for ln in blk.splitlines():
                p = ln.split()
                if len(p) == 7 and p[0] == 'C':
                    try:
                        x = float(p[1])
                    except ValueError:
                        x = np.nan
                    break
            if m:
                ts.append(float(m.group(1)))
                xs.append(x)
    return np.array(ts), np.array(xs)

files = [f'fort.{n}' for n in range(1001, 1021)]
fig, ax = plt.subplots(figsize=(10, 6))
n_inf_nan = 0
for i, f in enumerate(files):
    ts, xs = parse_first_c(f)
    bad = ~np.isfinite(xs)
    n_inf_nan += int(bad.sum())
    ax.plot(ts / 1000.0, xs, lw=0.8, alpha=0.7, color='#EE854A')
    tb = ts[bad]
    if len(tb):
        ax.plot(tb / 1000.0, np.zeros_like(tb), '|', ms=6, color='crimson',
                alpha=0.3)
print(f'total non-finite (Inf/NaN) samples across 20 files: {n_inf_nan}')

ax.set_ylim(-10, 10)
ax.set_xlabel('time t [ps]')
ax.set_ylabel('incident-atom x [$\\AA$]')
ax.set_title('rigid_surface_defect - DEFECT: non-finite coordinates from step 0\n'
             '(20 trajectories; crimson marks = Infinity/NaN samples, '
             'no finite coordinate is ever produced)', fontsize=11)
ax.grid(alpha=0.25)
fig.text(0.99, 0.01, 'All samples non-finite: integrator produced Infinity at step 0 '
         'and NaN thereafter (H(eV)=Infinity).', ha='right', fontsize=7, color='crimson')
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig('fig_long_trajectory.png', dpi=150)
print('saved fig_long_trajectory.png')
