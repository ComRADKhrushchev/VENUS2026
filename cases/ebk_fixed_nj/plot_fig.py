#!/usr/bin/env python
"""ebk_fixed_nj: NT=1000 EBK fixed (n=3, J=2) sampling statistics.

Data: stdout_nt1000.txt (CHOSEN EROTA/EVIBA echoes, NT=1000) + per-trajectory
bond length r (fort.1001..fort.2000 step-0 'C' coordinate lines, bohr).
Theory: fixed EBK state -> EVIBA exactly constant (22.313 kcal/mol, harmonic
closure); EROTA = Erot(J=2) fixed, but each trajectory draws a different
rotational phase so the instantaneous centrifugal stretch r varies within the
classically allowed range [r-, r+] of the J=2 effective potential.
"""
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BOHR = 0.529177
txt = open('stdout_nt1000.txt', encoding='utf-8', errors='ignore').read()
EROTA = np.array([float(m) for m in re.findall(r'CHOSEN:\s+EROTA =\s*([-\d.E+]+)', txt)])
EVIBA = np.array([float(m) for m in re.findall(r'CHOSEN:\s+EVIBA =\s*([-\d.E+]+)', txt)])
N = len(EROTA)

# bond lengths at step 0 from per-trajectory files
rs = []
for it in range(1, N + 1):
    with open(f'fort.{1000+it}', encoding='utf-8', errors='ignore') as fh:
        blk = fh.read(2000)
    coords = []
    for ln in blk[:blk.find('--- step', 10)].splitlines():
        parts = ln.split()
        if len(parts) == 7 and parts[0] == 'C':
            try:
                coords.append([float(x) for x in parts[1:4]])
            except ValueError:
                pass
    if len(coords) >= 2:
        rs.append(np.linalg.norm(np.array(coords[0]) - np.array(coords[1])) * BOHR)
rs = np.array(rs)

fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
# 1) EVIBA fixed value
ax = axes[0]
ax.hist(EVIBA, bins=30, color='#EE854A', edgecolor='k')
ax.axvline(22.313, color='r', ls='--', lw=2, label='EBK fixed state E(n=3)')
ax.set_xlabel('EVIBA [kcal/mol]')
ax.set_ylabel('trajectories')
ax.set_title(f'Vibrational energy: sd={EVIBA.std():.4f} (fixed state, no spread)')
ax.legend(fontsize=8)
# 2) EROTA distribution
ax = axes[1]
ax.hist(EROTA, bins=30, color='#4878CF', edgecolor='k')
ax.axvline(EROTA.mean(), color='r', ls='--', lw=2,
           label=f'mean={EROTA.mean():.3f}')
ax.set_xlabel('EROTA [kcal/mol]')
ax.set_ylabel('trajectories')
ax.set_title(f'Rotational energy (J=2 state, phase-sampled): [{EROTA.min():.3f},{EROTA.max():.3f}]')
ax.legend(fontsize=8)
# 3) bond-length sampling vs classical turning band
ax = axes[2]
ax.hist(rs, bins=40, density=True, color='#6ACC64', edgecolor='k', alpha=0.85,
        label=f'sampled r (N={len(rs)})')
rmin, rmax = rs.min(), rs.max()
ax.axvspan(rmin, rmax, color='r', alpha=0.08)
ax.set_xlabel('bond length r at t=0 [Angstrom]')
ax.set_ylabel('probability density')
ax.set_title(f'Instantaneous bond length: r in [{rmin:.4f},{rmax:.4f}] A\n'
             f'(re=0.7414 A; centrifugal stretch of J=2 orbit)')
ax.legend(fontsize=8)

fig.suptitle('EBK fixed (n=3, J=2) sampling (MORSE H$_2$, NT=1000):\n'
             'fixed-state energies + classical orbit phase sampling', fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.88))
fig.savefig('fig_sampling_stats.png', dpi=150)
print('saved; EVIBA mean/sd', EVIBA.mean(), EVIBA.std(),
      'EROTA mean/sd', EROTA.mean(), EROTA.std(), 'r range', rmin, rmax)
