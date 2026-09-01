#!/usr/bin/env python
"""microcanonical_normalmodes: NT=1000 final internal-energy statistics.

Data: stdout_nt1000.txt per-trajectory 'INTERNAL ENERGY = <val> KCAL/MOL'
(final of the iterative orthant rescaling echoes per trajectory) and the
'CHOSEN: EVIBA' final vibrational energies.
Theory: microcanonical ensemble -> final internal energy pinned at
ENMT = 1.0 kcal/mol within the acceptance tolerance of +-0.1%.
"""
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ENMT = 1.0
TOL = 0.001  # +-0.1 %

txt = open('stdout_nt1000.txt', encoding='utf-8', errors='ignore').read()
# last INTERNAL ENERGY echo per trajectory (final accepted value):
# split file at trajectory boundaries (SELECT lines) and take last per block
sel_split = re.split(r'\s+SELECT:', txt)
eint = []
for blk in sel_split[1:]:
    m = re.findall(r'INTERNAL ENERGY =\s+([-\d.E+]+)', blk)
    if m:
        eint.append(float(m[-1]))
eint = np.array(eint)
eviba = np.array([float(m) for m in re.findall(
    r'CHOSEN:\s+EVIBA =\s*([-\d.E+]+)', txt)])
print(f'parsed {len(eint)} final INTERNAL ENERGY, {len(eviba)} EVIBA')

fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

ax = axes[0]
ax.hist(eint, bins=60, color='#4878CF', edgecolor='k', alpha=0.85,
        label=f'sampled final E_int (N={len(eint)})')
ax.axvline(ENMT, color='r', ls='--', lw=2, label='ENMT = 1.0 kcal/mol')
ax.axvline(ENMT * (1 - TOL), color='orange', ls=':', lw=1.8,
           label='+-0.1% tolerance')
ax.axvline(ENMT * (1 + TOL), color='orange', ls=':', lw=1.8)
ax.set_xlabel('final internal energy [kcal/mol]')
ax.set_ylabel('trajectories')
rel = (eint - ENMT) / ENMT
frac_in = float(np.mean(np.abs(rel) <= TOL))
ax.set_title(f'E_int: mean={eint.mean():.6f}, sd={eint.std():.2e}, '
             f'max|dev|={np.abs(rel).max():.2e}\n{frac_in*100:.1f}% within +-0.1% band')
ax.legend(fontsize=9)

ax = axes[1]
ax.hist(rel * 100, bins=60, color='#6ACC64', edgecolor='k', alpha=0.85,
        label='relative deviation')
ax.axvline(0, color='r', ls='--', lw=2)
ax.axvline(-0.1, color='orange', ls=':', lw=1.8)
ax.axvline(+0.1, color='orange', ls=':', lw=1.8)
ax.set_xlabel('(E_int - ENMT)/ENMT  [%]')
ax.set_ylabel('trajectories')
ax.set_title(f'Relative deviation from ENMT (mean={rel.mean()*100:+.4f}%)')

fig.suptitle('Microcanonical normal-mode sampling (LEPS H$_3$, ENMT=1.0 kcal/mol, NT=1000):\n'
             'final internal energy pinned at target within tolerance', fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.86))
fig.savefig('fig_sampling_stats.png', dpi=150)
print('saved fig_sampling_stats.png')
print(f'E_int mean={eint.mean():.6f} sd={eint.std():.3e} '
      f'in-band={frac_in*100:.1f}%')
