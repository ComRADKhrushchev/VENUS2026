"""LEPS H3 potential-surface verification figures.

Reproduces the CGM LEPS triatomic surface (test_potentials.f90 leps_vg:
V = 0.5*sum(Q_p) - 0.5*sqrt(sum((A_p - Q_p)^2)), pair branches from
De/re/a/Sato-Delta) on two standard triatomic cuts:
  (a) collinear H+H2 (r13 = r12 + r23) - reaction-path surface
  (b) fixed r13 = 1.7599 A cut through the equilateral minimum
Anchors: equilateral well r=1.760 A / -6.381 eV (matches program VZERO
-147.147 kcal/mol to 5 digits); symmetric-path barrier 11.8 kcal/mol above
the H2+H asymptote (-4.05 eV).
"""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

De, RE, A, DELTA = 4.746, 1.401, 1.028, 0.164


def leps_branch(r):
    # exact transcription of leps_branch in test_potentials.f90
    # (note the +Delta*De constant in the A branch)
    em = math.exp(-A * (r - RE))
    qe = 0.5 * De * ((3 + DELTA) * em * em - (2 + 2 * DELTA) * em)
    ae = 0.5 * De * ((3 - DELTA) * em * em + (2 - 2 * DELTA) * em) + DELTA * De
    return qe, ae


def leps_v3(r12, r23, r13):
    qq = np.array([leps_branch(x) for x in (r12, r23, r13)])
    q, aa = qq[:, 0], qq[:, 1]
    return 0.5 * q.sum() - 0.5 * np.sqrt(((aa - q) ** 2).sum())


n = 300
r = np.linspace(0.5, 4.0, n)
R12, R23 = np.meshgrid(r, r)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.6))

V1 = np.array([[leps_v3(x, y, x + y) for x in r] for y in r])
im1 = ax1.contourf(R12, R23, V1, levels=np.linspace(-6.5, -1.0, 40),
                   cmap='viridis', extend='both')
cs1 = ax1.contour(R12, R23, V1, levels=np.arange(-6.4, -0.9, 0.4),
                  colors='w', linewidths=0.6)
ax1.clabel(cs1, fmt='%.1f', fontsize=7)
ax1.plot(1.607, 1.607, 'r*', ms=14, label='collinear symmetric min (1.607, -6.38 eV)')
ax1.plot(0.7414, 3.6, 'w^', ms=9)
ax1.plot(3.6, 0.7414, 'w^', ms=9, label='H2+H asymptote (-4.05 eV)')
ax1.plot([0.5, 4.0], [0.5, 4.0], 'r--', lw=0.8, alpha=0.7, label='symmetric path')
ax1.set_xlabel('r12 [A]'); ax1.set_ylabel('r23 [A]')
ax1.set_title('(a) Collinear H+H2: r13 = r12+r23\nbarrier 11.8 kcal/mol on symmetric path')
ax1.legend(fontsize=8, loc='upper right')
plt.colorbar(im1, ax=ax1, label='V [eV]')

r13_fix = 1.7599
V2 = np.array([[leps_v3(x, y, r13_fix) for x in r] for y in r])
im2 = ax2.contourf(R12, R23, V2, levels=np.linspace(-6.5, -1.0, 40),
                   cmap='viridis', extend='both')
cs2 = ax2.contour(R12, R23, V2, levels=np.arange(-6.4, -0.9, 0.4),
                  colors='w', linewidths=0.6)
ax2.clabel(cs2, fmt='%.1f', fontsize=7)
ax2.plot(r13_fix, r13_fix, 'r*', ms=14, label='equilateral min (1.760, -6.38 eV)')
ax2.set_xlabel('r12 [A]'); ax2.set_ylabel('r23 [A]')
ax2.set_title('(b) Fixed r13 = 1.7599 A cut\nMorse walls + equilateral well')
ax2.legend(fontsize=8, loc='upper right')
plt.colorbar(im2, ax=ax2, label='V [eV]')

fig.suptitle('LEPS H3 potential surfaces (CGM: De=4.746 eV, re=1.401 A, '
             'a=1.028 1/A, Sato D=0.164)', fontsize=11)
plt.tight_layout()
plt.savefig('fig_leps_surfaces.png', dpi=150, bbox_inches='tight',
            facecolor='white')
print('saved fig_leps_surfaces.png')
