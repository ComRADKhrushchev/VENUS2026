#!/usr/bin/env python
"""morse_bootstrap: analytic Morse curve vs fort.1001 total energy H(eV).

V(r) = De*(exp(-a*(r-re)) - 1)^2 - De, with De=4.746 eV, re=1.401 A, a=1.028/A,
plotted for r in [0.8, 8] A. Overlay: fort.1001 records of the pair distance
r(t) (two C atoms) vs H(eV) (constant -0.06838 eV).
"""
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

De, re_, a = 4.746, 1.401, 1.028

r = np.linspace(0.8, 8.0, 1200)
em = np.exp(-a * (r - re_))
V = De * (em - 1.0) ** 2 - De

# --- fort.1001: pair distance and H(eV) per step ---
ts, rr, HH = [], [], []
with open('fort.1001', encoding='utf-8', errors='ignore') as fh:
    for blk in fh.read().split('--- step')[1:]:
        m = re.search(r't\(fs\)=\s*(-?[\d.]+)', blk)
        h = re.search(r'H\(eV\)\s*(-?[\d.]+)', blk)
        cs = []
        for ln in blk.splitlines():
            p = ln.split()
            if len(p) == 7 and p[0] == 'C':
                try:
                    cs.append([float(x) for x in p[1:4]])
                except ValueError:
                    pass
        if m and h and len(cs) >= 2:
            ts.append(float(m.group(1)))
            rr.append(np.linalg.norm(np.array(cs[0]) - np.array(cs[1])))
            HH.append(float(h.group(1)))
ts, rr, HH = np.array(ts), np.array(rr), np.array(HH)
print(f'fort.1001: {len(ts)} steps, r=[{rr.min():.4f},{rr.max():.4f}] A, '
      f'H const={HH[0]:.5f} eV (unique: {len(set(np.round(HH,6)))==1})')
r_pred = r[np.argmin(np.abs(V - HH[0]))]
consistent = abs(r_pred - rr.mean()) < 0.1
print(f'Morse V(r)={HH[0]:.5f} eV at r={r_pred:.4f} A vs data r~{rr.mean():.4f} A '
      f'-> {"consistent" if consistent else "MISMATCH"}')

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(r, V, lw=1.8, color='#4878CF', label='Morse: De=4.746 eV, r$_e$=1.401 $\\AA$, a=1.028 $\\AA^{-1}$')
if consistent:
    ax.plot(rr, HH, '.', ms=3, color='#EE854A', alpha=0.5,
            label=f'fort.1001 H(eV) vs pair r(t), const {HH[0]:.5f} eV')
else:
    ax.plot([], [])  # keep color cycle stable
ax.axhline(0.0, color='gray', lw=0.6, ls=':')
ax.set_xlabel('r [$\\AA$]')
ax.set_ylabel('V(r), H [eV]')
ax.set_title('morse_bootstrap: analytic Morse potential vs trajectory total energy '
             '(fort.1001, 0.5 ps)', fontsize=11)
ax.legend(fontsize=9)
ax.grid(alpha=0.25)
if not consistent:
    fig.text(0.99, 0.01,
             f'fort.1001 H(eV) is constant {HH[0]:.5f} eV over r=[{rr.min():.2f},{rr.max():.2f}] A '
             f'and does not follow V(r) (V={HH[0]:.5f} eV only at r={r_pred:.2f} A): '
             'scatter omitted, analytic curve only.', ha='right', fontsize=7, color='crimson')
fig.tight_layout(rect=(0, 0.05, 1, 1))
fig.savefig('fig_potential_curve.png', dpi=150)
print('saved fig_potential_curve.png')
