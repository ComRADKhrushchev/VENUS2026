#!/usr/bin/env python
"""orthant_sampling: NT=1000 orthant-sampling energy statistics.

Data: stdout_nt1000.txt per-trajectory 'CHOSEN: EVIBA' (internal vibrational
energy, target HSCALE=10.185254 kcal/mol) and 'CHOSEN: EROTA' (target
EROTT=0.185254 kcal/mol).
Theory: orthant rescaling pins E_int at HSCALE and Erot at EROTT, so the
sampled values concentrate on those fixed values.
"""
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HSCALE = 10.185254   # kcal/mol (README value; run input used HSCALE=10.0, see note)
EROTT = 0.185254     # kcal/mol, fixed rotational-energy target

txt = open('stdout_nt1000.txt', encoding='utf-8', errors='ignore').read()
eviba = np.array([float(m) for m in re.findall(r'CHOSEN:\s+EVIBA =\s*([-\d.E+]+)', txt)])
erota = np.array([float(m) for m in re.findall(r'CHOSEN:\s+EROTA =\s*([-\d.E+]+)', txt)])
print(f'parsed {len(eviba)} EVIBA, {len(erota)} EROTA')

fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

ax = axes[0]
ax.hist(eviba, bins=40, color='#4878CF', edgecolor='k', alpha=0.85,
        label=f'sampled E_int (N={len(eviba)})')
ax.axvline(HSCALE, color='r', ls='--', lw=2, label=f'HSCALE = {HSCALE}')
ax.set_xlabel('internal energy EVIBA [kcal/mol]')
ax.set_ylabel('trajectories')
ax.set_title(f'E_int: mean={eviba.mean():.4f}, sd={eviba.std():.4f} '
             f'(input HSCALE=10.0)\n'
             f'max dev from HSCALE = {np.abs(eviba - HSCALE).max():.4f}')
ax.legend(fontsize=9)

ax = axes[1]
ax.hist(erota, bins=40, color='#EE854A', edgecolor='k', alpha=0.85,
        label=f'sampled E_rot (N={len(erota)})')
ax.axvline(EROTT, color='r', ls='--', lw=2, label=f'EROTT = {EROTT}')
ax.set_xlabel('rotational energy EROTA [kcal/mol]')
ax.set_ylabel('trajectories')
ax.set_title(f'E_rot: mean={erota.mean():.4f}, sd={erota.std():.4f}\n'
             f'max dev from EROTT = {np.abs(erota - EROTT).max():.4f}')
ax.legend(fontsize=9)

fig.suptitle('Orthant sampling (LEPS H$_3$, NT=1000):\n'
             'internal energy pinned at HSCALE, rotational energy at EROTT', fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.86))
fig.savefig('fig_sampling_stats.png', dpi=150)
print('saved fig_sampling_stats.png')
print(f'EVIBA mean={eviba.mean():.4f} sd={eviba.std():.4f} range=[{eviba.min():.3f},{eviba.max():.3f}]')
print(f'EROTA mean={erota.mean():.4f} sd={erota.std():.4f} range=[{erota.min():.3f},{erota.max():.3f}]')
