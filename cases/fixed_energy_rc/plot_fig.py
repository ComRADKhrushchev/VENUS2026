#!/usr/bin/env python
"""fixed_energy_rc: NT=1000 fixed-energy sampling incl. reaction coordinate.

Data: stdout_nt1000.txt (per-trajectory 'REACTION COORDINATE ENERGY' and
'CHOSEN: EVIBA/EROTA' echoes, NT=1000).
Theory: E_RC uniform on (0, ENMT) would give flat marginal; the program
implements Beyer-Swinehart grid state counting, so marginal follows the grid
PMF (README: joint chi2=6.2/12df). We overlay the uniform reference and mark
the KS statistic of E_RC against Uniform(0, ENMT_A=20 kcal/mol).
Energy closure: EVIBA + E_RC == 20.0 kcal/mol.
"""
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

txt = open('stdout_nt1000.txt', encoding='utf-8', errors='ignore').read()
ERC = np.array([float(m) for m in re.findall(r'REACTION COORDINATE ENERGY =\s*([-\d.E+]+)', txt)])
EVIBA = np.array([float(m) for m in re.findall(r'CHOSEN:\s+EVIBA =\s*([-\d.E+]+)', txt)])
N = len(ERC)
ENMT = 20.0

ks = stats.kstest(ERC / ENMT, 'uniform')
closure = (EVIBA + ERC) - ENMT

fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
ax = axes[0]
ax.hist(ERC, bins=25, density=True, color='#4878CF', edgecolor='k', alpha=0.85,
        label='sampled E$_{RC}$')
ax.axhline(1.0 / ENMT, color='r', ls='--', lw=2, label='uniform 1/ENMT reference')
ax.set_xlabel('reaction-coordinate energy E$_{RC}$ [kcal/mol]')
ax.set_ylabel('probability density')
ax.set_title(f'E$_{{RC}}$ marginal (N={N})')
ax.text(0.98, 0.95, f'KS vs U(0,20): D={ks.statistic:.4f}, p={ks.pvalue:.3f}',
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        bbox=dict(fc='wheat', alpha=0.6))
ax.legend(fontsize=8)

ax = axes[1]
ax.hist(closure, bins=30, color='#EE854A', edgecolor='k')
ax.set_xlabel('EVIBA + E$_{RC}$ - ENMT [kcal/mol]')
ax.set_ylabel('trajectories')
ax.set_title(f'Energy closure: max|dev|={np.abs(closure).max():.5f} kcal/mol')

ax = axes[2]
ax.scatter(ERC, EVIBA, s=4, alpha=0.4, color='#6ACC64', edgecolors='none')
xx = np.linspace(0, ENMT, 10)
ax.plot(xx, ENMT - xx, 'r--', lw=2, label='EVIBA = ENMT - E$_{RC}$')
ax.set_xlabel('E$_{RC}$ [kcal/mol]')
ax.set_ylabel('EVIBA [kcal/mol]')
ax.set_title('Microcanonical split: EVIBA vs E$_{RC}$')
ax.legend(fontsize=8)

fig.suptitle('Fixed-energy sampling incl. reaction coordinate (LEPS H$_3$, ENMT=20 kcal/mol, NT=1000)',
             fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig('fig_sampling_stats.png', dpi=150)
print('saved; N=', N, 'KS D=', ks.statistic, 'p=', ks.pvalue,
      'closure max', np.abs(closure).max())
