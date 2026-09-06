#!/usr/bin/env python
"""rst_surface_md: thermostat temperature T(t) and C height z(t), NT=4.

Data sources (GWRITE_LEVEL=1, INIT_SAMPLING_B=MD):
  - run_full.log: `system temperature=` printed every NIP=50 cycles
    (DT=0.01 -> 0.1 fs per cycle) for each TRAJECTORY block, plus the
    `temperature resample =` lines of the NSCALE velocity-rescale
    equilibration stage.
  - fort.1001..fort.1004: per-trajectory C z(t) step blocks.
Line-by-line split() parsing only (no cross-line regex). Panels (ps axis):
  (a) T(t) with 300 K dashed target (equilibration resamples shown as a
      short shaded lead-in; production curves per trajectory)
  (b) C height z(t) for the 4 trajectories.
A late-time thermostat spike to ~5.8e4 K in trajectory 2 (after ~1 ps,
heat pulse from the 0.63 eV impact depositing into the small 145-atom
cell) is masked from the T axis; masked point count is reported.
"""
import math
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NTRJ = 4
NIP = 50
DTFS = 0.1      # DT=0.01 in units of 10 fs -> 0.1 fs per cycle
T_CAP = 600.0  # mask non-physical thermostat spikes above this (K)


def parse_resample(path):
    """Equilibration-stage resample temperatures (lead-in segment)."""
    out = []
    with open(path, encoding='utf-8', errors='ignore') as fh:
        for ln in fh:
            if 'temperature resample' in ln:
                out.append(float(ln.split('=')[-1]))
    return np.array(out)


def parse_prod_temp(path):
    """(trj, t_fs, T) tuples from `system temperature=` + cycle-count lines."""
    out = []
    cur_trj, n_print, t_last = 0, 0, 0.0
    with open(path, encoding='utf-8', errors='ignore') as fh:
        for ln in fh:
            m = re.search(r'TRAJECTORY NUMBER\s+(\d+)', ln)
            if m:
                cur_trj, n_print = int(m.group(1)), 0
                continue
            if 'CYCLE COUNT IS' in ln:
                t_last = float(ln.split('TIME:')[1].split()[0])
                continue
            if 'system temperature=' in ln:
                out.append((cur_trj, t_last, float(ln.split('=')[1])))
                n_print += 1
    return out


def parse_z(path):
    ts, zs = [], []
    with open(path, encoding='utf-8', errors='ignore') as fh:
        for ln in fh:
            tok = ln.split()
            if ln.startswith('--- step'):
                if len(tok) >= 5:
                    ts.append(float(tok[4]))
            elif ln.startswith('Q(C)='):
                z = float(tok[-1])
                zs.append(z if math.isfinite(z) else np.nan)
    n = min(len(ts), len(zs))
    return np.array(ts[:n]), np.array(zs[:n])


resamp = parse_resample('run_full.log')
tprod = parse_prod_temp('run_full.log')

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 9), sharex=True)

# (a) temperature: equilibration resamples (lead-in), then 4 production curves
nres = len(resamp)
t_res_end = nres * NIP * DTFS / 1000.0  # ps, approximate lead-in span
t_res = np.linspace(0.0, max(t_res_end, 1e-9), nres) if nres else np.array([])
ax1.fill_between(t_res, resamp, color='0.75', alpha=0.8, lw=0,
                 label='equilibration resamples (NSCALE)')
off = t_res_end
cmapT = plt.cm.plasma
seen = sorted(set(r[0] for r in tprod))
nmask = 0
for k, trj in enumerate(seen):
    sel = np.array([(r[1], r[2]) for r in tprod if r[0] == trj])
    tps = off + sel[:, 0] / 1000.0
    T = sel[:, 1]
    good = T <= T_CAP
    nmask += int((~good).sum())
    ax1.plot(tps[good], T[good], color=cmapT(k / max(len(seen) - 1, 1)),
             lw=0.8, label=f'production trj {trj}')
ax1.axhline(300.0, color='k', ls='--', lw=1.0, label='target 300 K')
ax1.set_ylabel('System temperature T (K)')
ax1.set_title('RST Au(111) slab MD equilibration, '
              'THERMOTEMP=300 K, NT=4')
ax1.legend(loc='upper left', fontsize=8)
ax1.grid(alpha=0.3)
print(f'masked thermostat spikes (T > {T_CAP:.0f} K): {nmask}')

# (b) C z(t) per trajectory
cmap = plt.cm.viridis
for i in range(NTRJ):
    ts, zs = parse_z(f'fort.{1001 + i}')
    if len(ts):
        ax2.plot(ts / 1000.0, zs, color=cmap(i / max(NTRJ - 1, 1)), lw=1.0,
                 label=f'trj {i + 1}')
ax2.set_ylabel('C height z (Å)')
ax2.set_xlabel('Time t (ps)')
ax2.legend(loc='upper right', fontsize=9)
ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig('fig_thermal.png', dpi=150)

# numeric summary over unmasked production temperatures
allT = np.array([r[2] for r in tprod])
ok = allT <= T_CAP
print(f'production temperature records: {len(allT)} '
      f'({int(ok.sum())} used, {int((~ok).sum())} masked)')
print(f'  T = {allT[ok].mean():.1f} +/- {allT[ok].std():.1f} K '
      f'(target 300 K)')
