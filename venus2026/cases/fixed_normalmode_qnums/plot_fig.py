#!/usr/bin/env python
"""fixed_normalmode_qnums: NT=1000 fixed normal-mode quantum numbers.

Data: stdout_nt1000.txt (CHOSEN EVIBA echoes, NT=1000; quantum numbers are
fixed input n=(1,2,3) echoed as ANQ each trajectory).
Theory: fixed state (n1,n2,n3)=(1,2,3) -> harmonic energy
E = sum (n_i+1/2)*h*nu_i exactly constant; sampled EVIBA histogram should be
a delta-like spike at the analytic value with only integration-rounding
spread (README: 29.636 vs 29.642 kcal/mol, spread 0.059).
"""
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

txt = open('stdout_nt1000.txt', encoding='utf-8', errors='ignore').read()
EVIBA = np.array([float(m) for m in re.findall(r'CHOSEN:\s+EVIBA =\s*([-\d.E+]+)', txt)])
N = len(EVIBA)

# harmonic closure for n=(1,2,3), nu=(1013.47,1013.47,1803.83) cm-1
NU = np.array([1013.47, 1013.47, 1803.83])
nq = np.array([1, 2, 3])
CM2KCAL = 0.002859144
E_theory = np.sum((nq + 0.5) * NU) * CM2KCAL

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
ax = axes[0]
ax.hist(EVIBA, bins=30, color='#4878CF', edgecolor='k')
ax.axvline(E_theory, color='r', ls='--', lw=2,
           label=f'harmonic E[(1,2,3)]={E_theory:.3f} kcal/mol')
ax.set_xlabel('EVIBA [kcal/mol]')
ax.set_ylabel('trajectories')
ax.set_title(f'Fixed quantum numbers: mean={EVIBA.mean():.3f}, sd={EVIBA.std():.4f}')
ax.legend(fontsize=9)

ax = axes[1]
ax.plot(EVIBA, np.linspace(1, N, N) / N, lw=1.2)
ax.set_xlabel('EVIBA sorted value [kcal/mol]')
ax.set_ylabel('ECDF')
ax.set_title(f'ECDF: tight spike at fixed state\n(spread {EVIBA.max()-EVIBA.min():.4f} kcal/mol)')

fig.suptitle('Fixed normal-mode quantum numbers (LEPS H$_3$, n=(1,2,3), NT=1000):\n'
             'state energy conservation across trajectories', fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.88))
fig.savefig('fig_sampling_stats.png', dpi=150)
print('saved; mean/sd', EVIBA.mean(), EVIBA.std(), 'theory', E_theory)
