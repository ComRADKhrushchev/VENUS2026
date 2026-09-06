#!/usr/bin/env python
"""rst_beam_scattering: C z(t) and r_min(t) for NT=20 RST trajectories.

Data: fort.1001..fort.1020, one file per trajectory (GWRITE_LEVEL=1). Each
`--- step` block header carries E0/T/H(eV), Q(C), V(C) and Au_nearest/r_min.
All 20 trajectories run the full NS=30000 cycles (=3000 fs) -> adsorbed.
Parsed line-by-line with split() (no cross-line regex). Panels (ps axis):
  (a) C height z(t), 20 trajectories overlaid + adsorption-well annotations
  (b) C-Au nearest distance r_min(t)
"""
import math

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NTRJ = 20
Z_EQ = 2.90      # atop adsorption well position (A), from smoke test
E_ADS = -2.81    # adsorption well depth (eV)
E_INC = 0.63     # incident kinetic energy (eV) = 14.528 kcal/mol
START_Z = 6.5    # asymptotic-region start height (A)


def parse(path):
    """Line-by-line parse of one fort.1NN trajectory file."""
    ts, zs, rmins = [], [], []
    with open(path, encoding='utf-8', errors='ignore') as fh:
        for ln in fh:
            tok = ln.split()
            if ln.startswith('--- step'):
                # ['---','step','0','t(fs)=','0.000','---']
                if len(tok) >= 5:
                    ts.append(float(tok[4]))
            elif ln.startswith('Q(C)='):
                # ['Q(C)=','20.36738','0.16956','2.50000'] -> z = last value
                z = float(tok[-1])
                zs.append(z if math.isfinite(z) else np.nan)
            elif 'Au_nearest' in ln and 'r_min' in ln:
                rmins.append(float(ln.split('r_min')[1]))
    n = min(len(ts), len(zs), len(rmins))
    return np.array(ts[:n]), np.array(zs[:n]), np.array(rmins[:n])


cmap = plt.cm.viridis
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 9), sharex=True)
n_scat = n_ads = 0  # filled in the loop below

stats = []
for i in range(NTRJ):
    ts, zs, rmins = parse(f'fort.{1001 + i}')
    if len(ts) == 0:
        continue
    tps = ts / 1000.0  # fs -> ps
    full = len(ts) >= 601  # NIP=50 over NS=30000 -> 601 blocks if full run
    c = 'crimson' if full else cmap(i / max(NTRJ - 1, 1))
    ax1.plot(tps, zs, color=c, lw=0.9)
    ax2.plot(tps, rmins, color=c, lw=0.9)
    if full: n_ads += 1
    else: n_scat += 1
    stats.append((i + 1, len(ts), zs.min(), zs[-1], rmins.min(), full))

# adsorption-well annotations on the z panel
ax1.axhline(Z_EQ, color='crimson', ls='--', lw=1.0,
            label=f'atop well z$_{{eq}}$={Z_EQ:.2f} Å')
ax1.text(0.02, 0.96,
         f'E$_{{ads}}$ = {E_ADS:.2f} eV  vs  E$_{{inc}}$ = {E_INC:.2f} eV\n'
         f'scattered {n_scat}/{NTRJ}, adsorbed {n_ads}/{NTRJ}\n'
         '(phonon channel of elastic slab carries away energy)',
         transform=ax1.transAxes, va='top', ha='left', fontsize=9,
         bbox=dict(boxstyle='round', fc='wheat', alpha=0.85))
ax1.set_ylabel('C height z (Å)')
ax1.set_title('RST C/Au(111) beam scattering (asymptotic start z=6.5 A), '
              f'E$_{{rel}}$=14.528 kcal/mol, NT={NTRJ}; red=adsorbed')
ax1.legend(loc='upper right', fontsize=9)
ax2.set_ylabel('C–Au nearest r$_{min}$ (Å)')
ax2.set_xlabel('Time t (ps)')
for ax in (ax1, ax2):
    ax.set_xlim(0, 3.05)
    ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig('fig_scattering.png', dpi=150)

full = sum(1 for s in stats if s[5])
early = [s for s in stats if not s[5]]
print(f'trajectories parsed: {len(stats)}')
print(f'full-NS (adsorbed) runs: {full}; early-terminated: {len(early)}')
zmins = [s[2] for s in stats]
zends = [s[3] for s in stats]
rmin_all = [s[4] for s in stats]
print(f'z_min over all: {min(zmins):.3f}-{max(zmins):.3f} A')
print(f'final z over all: {min(zends):.3f}-{max(zends):.3f} A')
print(f'r_min over all: {min(rmin_all):.3f}-{max(rmin_all):.3f} A')
