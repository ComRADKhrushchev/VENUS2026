#!/usr/bin/env python
"""mb_thermal: NT=1000 Maxwell-Boltzmann thermal (n,J) sampling statistics.

Data: stdout_nt1000.txt (per-trajectory 'NNA =' vibrational quantum number n and
'JA =' rotational quantum number J, TVIB=TROT=2000 K on MORSE H2).
Theory:
  n: geometric distribution P(n) = (1-q) q^n, q = exp(-1.43878*nu/T) = 0.1876
     with nu = 2326.4 cm-1, T = 2000 K  ->  <n> = q/(1-q) = 0.2308
  J: P(J) prop. WGT(J) (2J+1) exp(-B J(J+1)/kT), B/kT given via B=0.012259
     (rotational constant in units matched to kT at 2000 K), nuclear spin
     weights WGT(odd J)=0.75, WGT(even J)=0.25 for homonuclear H2 (ortho/para 3:1).
"""
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

T = 2000.0
NU = 2326.4            # cm-1, H2 Morse frequency
q = float(np.exp(-1.43878 * NU / T))          # 0.1876
B = 0.012259                                   # rotational constant (units of B*J(J+1)/kT)
WGT = {0: 0.25, 1: 0.75}                       # even/odd J spin weights

txt = open('stdout_nt1000.txt', encoding='utf-8', errors='ignore').read()
ns = np.array([int(m) for m in re.findall(r'NNA =\s+(\d+)', txt)])
js = np.array([int(m) for m in re.findall(r'JA =\s+(\d+)', txt)])
N = min(len(ns), len(js))
print(f'parsed N={len(ns)} n-values, {len(js)} J-values')

fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

# ---- panel 1: vibrational quantum number n vs geometric ----
ax = axes[0]
nmax = ns.max()
counts = np.bincount(ns, minlength=nmax + 1)
nn = np.arange(0, nmax + 1)
p = (1 - q) * q ** nn
p = p / p.sum()
expc = len(ns) * p
mask = expc >= 5
chi2 = float(np.sum((counts[mask] - expc[mask]) ** 2 / expc[mask]))
dof = int(mask.sum() - 1)
pval = float(1 - stats.chi2.cdf(chi2, dof))
ax.bar(nn, counts / len(ns), color='#4878CF', edgecolor='k', alpha=0.85,
       label=f'sampled n (N={len(ns)})')
ax.plot(nn, p, 'ro-', ms=6, lw=1.8, label=f'theory geometric, q={q:.4f}')
ax.set_xlabel('vibrational quantum number n')
ax.set_ylabel('probability')
ax.set_title(f'Vibrational n: <n>={ns.mean():.3f} vs theory q/(1-q)={q/(1-q):.3f}')
ax.text(0.97, 0.60, f'$\\chi^2$={chi2:.2f} (dof={dof})\np={pval:.3f}',
        transform=ax.transAxes, ha='right', va='top', fontsize=10,
        bbox=dict(fc='wheat', alpha=0.6))
ax.legend(fontsize=9)

# ---- panel 2: rotational quantum number J vs Boltzmann P(J) ----
ax = axes[1]
jmax = js.max()
jj = np.arange(0, jmax + 1)
wgt = np.where(jj % 2 == 1, WGT[1], WGT[0])
pj = wgt * (2 * jj + 1) * np.exp(-B * jj * (jj + 1))
pj = pj / pj.sum()
countsJ = np.bincount(js, minlength=jmax + 1)
expcJ = len(js) * pj
maskJ = expcJ >= 5
chi2J = float(np.sum((countsJ[maskJ] - expcJ[maskJ]) ** 2 / expcJ[maskJ]))
dofJ = int(maskJ.sum() - 1)
pvalJ = float(1 - stats.chi2.cdf(chi2J, dofJ))
meanJ_theory = float((jj * pj).sum())
ax.bar(jj, countsJ / len(js), color='#EE854A', edgecolor='k', alpha=0.85,
       label=f'sampled J (N={len(js)})')
ax.plot(jj, pj, 'ks-', ms=5, lw=1.8, mfc='w',
        label=r'theory $\propto$ WGT$(2J+1)e^{-BJ(J+1)}$')
ax.set_xlabel('rotational quantum number J')
ax.set_ylabel('probability')
ax.set_title(f'Rotational J: <J>={js.mean():.2f} vs theory {meanJ_theory:.2f}'
             '\n(spin weights: odd 0.75 / even 0.25)')
ax.text(0.97, 0.60, f'$\\chi^2$={chi2J:.2f} (dof={dofJ})\np={pvalJ:.3f}',
        transform=ax.transAxes, ha='right', va='top', fontsize=10,
        bbox=dict(fc='wheat', alpha=0.6))
ax.legend(fontsize=9)

fig.suptitle('Maxwell-Boltzmann thermal sampling (MORSE H$_2$, TVIB=TROT=2000 K, NT=1000):\n'
             'vibrational n vs geometric & rotational J vs Boltzmann with spin weights',
             fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.86))
fig.savefig('fig_sampling_stats.png', dpi=150)
print('saved fig_sampling_stats.png')
print(f'n: mean={ns.mean():.4f} theory={q/(1-q):.4f} chi2={chi2:.2f}/{dof} p={pval:.3f}')
print(f'J: mean={js.mean():.4f} theory={meanJ_theory:.4f} chi2={chi2J:.2f}/{dofJ} p={pvalJ:.3f}')
