#!/usr/bin/env python
"""qm_micro: NT=1000 microcanonical quantum-state degenerate-pair occupancy.

Data: fort.9 (per-trajectory normal-mode quantum-number triples) - the two
degenerate E' modes of the LEPS H3 sample n1/n2 from a twofold degenerate pair.
Theory: equiprobable 0.5/0.5 occupation between the degenerate partners.
Verified earlier: 507 / 493 split.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

nq = np.loadtxt('fort.9')          # (1000, 3) quantum numbers
N = nq.shape[0]

# which mode carries the excitation: mode index (0-based) of max quantum number
lead = np.argmax(nq, axis=1)
counts = np.bincount(lead, minlength=3)
labels = ["mode 1 (E')", "mode 2 (E')", "mode 3 (A1')"]
print('lead-mode counts:', counts, 'N =', N)

# degenerate pair = modes 1&2 (E'): occupancy vs 50/50 theory
pair_counts = counts[:2]
chi2 = float(np.sum((pair_counts - N / 2) ** 2 / (N / 2)))
pval = float(1 - stats.chi2.cdf(chi2, 1))

fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

ax = axes[0]
x = np.arange(2)
ax.bar(x, pair_counts / N, width=0.5, color=['#4878CF', '#EE854A'],
       edgecolor='k', alpha=0.85,
       label=[f'mode1 occupied: {pair_counts[0]}', f'mode2 occupied: {pair_counts[1]}'])
ax.axhline(0.5, color='r', ls='--', lw=2, label='theory 0.5 / 0.5')
for xi, c in zip(x, pair_counts):
    ax.text(xi, c / N + 0.01, f'{c}/{N}\n({c/N:.3f})', ha='center', fontsize=10)
ax.set_xticks(x)
ax.set_xticklabels(["mode 1 (E')", "mode 2 (E')"])
ax.set_ylabel('fraction of trajectories')
ax.set_ylim(0, 0.62)
ax.set_title(f"Degenerate E' pair occupancy: {pair_counts[0]}/{pair_counts[1]}")
ax.text(0.97, 0.95, f'$\\chi^2$={chi2:.2f} (dof=1)\np={pval:.3f}',
        transform=ax.transAxes, ha='right', va='top', fontsize=10,
        bbox=dict(fc='wheat', alpha=0.6))
ax.legend(fontsize=9, loc='upper center')

ax = axes[1]
x3 = np.arange(3)
ax.bar(x3, counts / N, width=0.55, color=['#4878CF', '#EE854A', '#6ACC64'],
       edgecolor='k', alpha=0.85, label='sampled lead mode (NT=1000)')
for xi, c in zip(x3, counts):
    ax.text(xi, c / N + 0.01, f'{c}', ha='center', fontsize=10)
ax.set_xticks(x3)
ax.set_xticklabels(labels)
ax.set_ylabel('fraction of trajectories')
ax.set_title(f'Leading-mode distribution (means n = {nq.mean(axis=0).round(3)})')
ax.legend(fontsize=9)

fig.suptitle('Quantum microcanonical sampling (LEPS H$_3$ degenerate E\' pair, NT=1000):\n'
             'degenerate-pair occupancy vs equiprobable 0.5/0.5', fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.86))
fig.savefig('fig_sampling_stats.png', dpi=150)
print('saved fig_sampling_stats.png')
print(f'pair occupancy {pair_counts}, chi2={chi2:.3f} p={pval:.3f}')
