#!/usr/bin/env python
"""boltzmann_vib: NT=1000 Boltzmann vibrational sampling statistics.

Data: fort.9 (per-trajectory normal-mode quantum numbers n1,n2,n3 for the
three LEPS H3 modes) from the NT=1000 run (stdout_nt1000.txt also kept).
Theory: geometric distribution P(n) = (1-q) q^n with
q = exp(-1.43878 * nu_i / T), nu = (1013.47, 1013.47, 1803.83) cm-1, T = 2000 K.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

NU = np.array([1013.47, 1013.47, 1803.83])   # cm-1 (mode 7/8 E' pair, mode 9 A1')
T = 2000.0                                    # K
LABELS = ['mode1 (E\', 1013.5)', 'mode2 (E\', 1013.5)', 'mode3 (A1\', 1803.8)']

nq = np.loadtxt('fort.9')                     # (1000, 3) quantum numbers
N = nq.shape[0]

fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharey=False)
for k, ax in enumerate(axes):
    n = nq[:, k].astype(int)
    q = np.exp(-1.43878 * NU[k] / T)
    nmax = n.max()
    counts = np.bincount(n, minlength=nmax + 1)
    ks = stats.kstest(n / (nmax + 1.0), lambda x: 1 - q ** np.maximum(1, np.ceil(x * (nmax + 1))))
    # chi2 on observed bins vs geometric
    nn = np.arange(0, nmax + 1)
    p = (1 - q) * q ** nn
    p = p / p.sum()
    expc = N * p
    mask = expc >= 5
    chi2 = np.sum((counts[mask] - expc[mask]) ** 2 / expc[mask])
    dof = mask.sum() - 1
    pval = 1 - stats.chi2.cdf(chi2, dof)
    ax.bar(nn, counts / N, color='#4878CF', edgecolor='k', alpha=0.8, label='sampled (NT=1000)')
    ax.plot(nn, p, 'ro-', ms=5, lw=1.5, label=f'theory geom q={q:.3f}')
    ax.set_xlabel(f'quantum number n  [{LABELS[k]}]')
    if k == 0:
        ax.set_ylabel('probability')
    ax.set_title(f'mode {k+1}: <n>={n.mean():.3f} vs q/(1-q)={q/(1-q):.3f}')
    ax.text(0.98, 0.60, f'$\\chi^2$={chi2:.2f} (dof={dof})\np={pval:.3f}',
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(fc='wheat', alpha=0.6))
    ax.legend(fontsize=8)
fig.suptitle('Boltzmann vibrational sampling (LEPS H$_3$, TVIB=2000 K, NT=1000):\n'
             'mode quantum numbers vs geometric distribution P(n)=(1-q)q$^n$', fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.90))
fig.savefig('fig_sampling_stats.png', dpi=150)
print('saved fig_sampling_stats.png')
print('means:', nq.mean(axis=0), 'theory:', np.exp(-1.43878*NU/T)/(1-np.exp(-1.43878*NU/T)))
for k in range(3):
    q = np.exp(-1.43878*NU[k]/T)
    n = nq[:, k].astype(int); nn = np.arange(0, n.max()+1)
    counts = np.bincount(n, minlength=n.max()+1)
    p = (1-q)*q**nn; p/=p.sum(); expc = N*p; m = expc>=5
    chi2 = ((counts[m]-expc[m])**2/expc[m]).sum()
    print(f'mode{k+1}: chi2={chi2:.2f} dof={m.sum()-1} p={1-stats.chi2.cdf(chi2,m.sum()-1):.3f}')
