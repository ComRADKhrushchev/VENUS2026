#!/usr/bin/env python
"""emt_surface_md: thermostat temperature T(t) and C height z(t), NT=4.

Data sources (GWRITE_LEVEL=1, NSEL_B=1):
  - fort.30: equilibration-stage SYSTEM TEMPERATURE history (one line per
    equilibration MD step; NSCALE=5000 velocity rescale, NEQUAL=5000).
  - run.log: production-stage `system temperature=` printed every NIP=50
    cycles (DT=0.01 -> 0.5 fs per cycle) per trajectory, together with the
    fort.10NN z(t) of the C atom.
Panels (ps axis): (a) T(t) + 300 K dashed target, (b) C z(t).
"""
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NTRJ = 4
NIP = 50
DTFS = 0.1  # DT=0.01 in units of 10 fs -> 0.1 fs per cycle


def parse_temp_log(path):
    """Return (stage, T, t_fs) lists: 'equil' from fort.30, 'prod' from run.log."""
    teq = []
    for ln in open('fort.30', encoding='utf-8', errors='ignore'):
        m = re.search(r'SYSTEM TEMPERATURE=\s*([\d.]+)', ln)
        if m:
            teq.append(float(m.group(1)))
    tprod = []
    cur_trj, n_print = 0, 0
    for ln in open(path, encoding='utf-8', errors='ignore'):
        m = re.search(r'TRAJECTORY NUMBER\s+(\d+)', ln)
        if m:
            cur_trj = int(m.group(1))
            n_print = 0
        elif 'system temperature=' in ln:
            m = re.search(r'system temperature=\s*([\d.]+)', ln)
            tprod.append((cur_trj, n_print * NIP * DTFS, float(m.group(1))))
            n_print += 1
    return np.array(teq), np.array(tprod)


def parse_z(path):
    ts, zs = [], []
    txt = open(path, encoding='utf-8', errors='ignore').read()
    for blk in txt.split('--- step')[1:]:
        m = re.search(r't\(fs\)=\s*(-?[\d.]+)', blk)
        q = re.search(r'Q\(C\)=\s*(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)', blk)
        if m and q:
            ts.append(float(m.group(1)))
            zs.append(float(q.group(3)))
    return np.array(ts), np.array(zs)


teq, tprod = parse_temp_log('run.log')

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 9), sharex=True)

# (a) temperature: equilibration phase then production phases (4 trajectories)
# equilibration time axis: assume same 0.1 fs/cycle step density
t_eq_ps = np.arange(len(teq)) * DTFS / 1000.0
ax1.plot(t_eq_ps, teq, color='tab:blue', lw=0.8, label='equilibration (fort.30)')
off = t_eq_ps[-1] if len(t_eq_ps) else 0.0
cmapT = plt.cm.plasma
ntrj_seen = sorted(set(int(r[0]) for r in tprod))
for k, trj in enumerate(ntrj_seen):
    sel = tprod[tprod[:, 0] == trj]
    ax1.plot(off + sel[:, 1] / 1000.0, sel[:, 2],
             color=cmapT(k / max(len(ntrj_seen) - 1, 1)), lw=0.8,
             label=f'production trj {trj}')
ax1.axhline(300.0, color='k', ls='--', lw=1.0, label='target 300 K')
ax1.set_ylabel('System temperature (K)')
ax1.set_title('EMT-NN Au(111) slab MD equilibration, THERMOTEMP=300 K, NT=4')
ax1.legend(loc='upper right', fontsize=8, ncol=2)
ax1.grid(alpha=0.3)

# (b) C z(t) per trajectory
cmap = plt.cm.viridis
for i in range(NTRJ):
    ts, zs = parse_z(f'fort.{1001 + i}')
    if len(ts):
        ax2.plot(ts / 1000.0, zs, color=cmap(i / max(NTRJ - 1, 1)), lw=1.0,
                 label=f'trj {i + 1}')
ax2.set_ylabel('C height z (Å)')
ax2.set_xlabel('Time (ps)')
ax2.legend(loc='upper right', fontsize=9)
ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig('fig_thermal.png', dpi=150)

# numeric summary
T_eq = np.array(teq)
last = T_eq[len(T_eq) // 2:] if len(T_eq) else T_eq
Tp = tprod[:, 2]
print(f'fort.30 equilibration samples: {len(T_eq)}')
print(f'  equil T: mean={T_eq.mean():.1f} sd={T_eq.std():.1f} K; '
      f'last-half mean={last.mean():.1f} K')
print(f'production samples: {len(Tp)}')
print(f'  prod T: mean={Tp.mean():.1f} sd={Tp.std():.1f} K; '
      f'last-quarter mean={Tp[len(Tp) * 3 // 4:].mean():.1f} K')
