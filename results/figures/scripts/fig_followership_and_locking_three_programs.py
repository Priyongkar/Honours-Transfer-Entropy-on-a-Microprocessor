"""Figure 5.x (lever L1): the followership penalty on five data sets, its gradient with the target's locking
on three recorded programs, and the same gradient in synthetic networks. Reads the hash-logged outputs
(A9, A11, A12, A15, R2); writes thesis/images/new/fig_c5_followership_locking.{pdf,png}."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os; H = Path(os.environ.get('HONOURS_DATA_ROOT', Path.home() / 'Documents/Honours 2')) / 'analyses'
NT, RP = H / 'new_tests/output', H / 'replication/output'
OUT = H.parent / 'thesis/images/new/fig_c5_followership_locking'
ORANGE, INDIGO, INK, MUTED = '#e8692d', '#4a5fc1', '#222222', '#666666'
PROG = {'Donkey Kong': '#222222', 'Pitfall': '#4a5fc1', 'Space Invaders': '#2a9d8f'}
RAMP = {0.3: '#f0955f', 0.6: '#d4561c', 0.9: '#8f3509'}; MARK = {0.3: 'o', 0.6: 's', 0.9: 'D'}

a9 = json.load(open(NT / 'A9_results.json')); a11 = json.load(open(NT / 'A11_results.json'))
a12 = json.load(open(NT / 'A12_results.json'))['A12']
r2 = {g: json.load(open(RP / f'R2_{f}_results.json'))['clock_aligned'] for g, f in (('Pitfall', 'Pitfall'), ('Space Invaders', 'SpaceInvaders'))}

plt.rcParams.update({'font.size': 11})
fig = plt.figure(figsize=(10.0, 8.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1])
ax1 = fig.add_subplot(gs[0, :]); ax2 = fig.add_subplot(gs[1, 0]); ax3 = fig.add_subplot(gs[1, 1])
axes = [ax1, ax2, ax3]
# (a) activity vs lag-1 TE within targets, five data sets
sets = [('Donkey Kong\n(recording)', a9['recording']['tests']), ('Pitfall', r2['Pitfall']['R2.4']),
        ('Space\nInvaders', r2['Space Invaders']['R2.4']), ('Chip,\nprogram 1', a9['chip_program1']['tests']),
        ('Chip,\nprogram 2', a9['chip_program2']['tests'])]
for x, (lab, t) in enumerate(sets):
    for name, col, dx, mk in (('desc_activity_within_target', INDIGO, -0.13, 'o'), ('desc_TE_within_target', ORANGE, 0.13, 's')):
        r = t[name]
        ax1.errorbar(x + dx, r['auc'], yerr=[[r['auc'] - r['lo']], [r['hi'] - r['auc']]], fmt=mk, color=col,
                     markersize=7.5, capsize=3, elinewidth=1.4, markeredgecolor='white')
    ax1.annotate('', xy=(x + 0.13, t['desc_TE_within_target']['auc']), xytext=(x - 0.13, t['desc_activity_within_target']['auc']),
                 arrowprops=dict(arrowstyle='->', color=MUTED, lw=1))
ax1.axvline(2.5, color='0.8', lw=0.8)
ax1.axhline(0.5, color='0.55', ls='--', lw=1); ax1.text(2.42, 0.503, 'chance', fontsize=8.5, color=MUTED, va='bottom', ha='right')
ax1.set_xticks(range(5)); ax1.set_xticklabels([s[0] for s in sets], fontsize=10.5)
ax1.set_xlim(-0.5, 4.5); ax1.set_ylim(0.3, 0.68)
ax1.set_ylabel('Within-target AUC,\nwired vs unwired')
ax1.set_title('(a) Followed sources lose rank', loc='left', fontsize=11)
ax1.plot([], [], 'o', color=INDIGO, label='activity alone'); ax1.plot([], [], 's', color=ORANGE, label='lag-1 TE')
ax1.legend(loc='upper right', frameon=False, fontsize=10)
# (b) A12 / R2.5: within-target AUC by the target's locking bin, three recorded programs
bins = ['mixed', 'locked', 'strict']; xs = [0, 1, 2]
dk = [a12[f'lag1_{b}'] for b in bins]
series = [('Donkey Kong', dk, -0.12)] + [(g, [r2[g]['R2.5'][b] for b in bins], dx) for g, dx in (('Pitfall', 0.0), ('Space Invaders', 0.12))]
for g, rs, dx in series:
    y = [r['auc'] for r in rs]
    ax2.errorbar([x + dx for x in xs], y, yerr=[[r['auc'] - r['lo'] for r in rs], [r['hi'] - r['auc'] for r in rs]],
                 fmt='-o', color=PROG[g], markersize=6.5, capsize=3, lw=1.6, elinewidth=1.2, markeredgecolor='white', label=g)
ax2.axhline(0.5, color='0.55', ls='--', lw=1)
ax2.set_xticks(xs); ax2.set_xticklabels(['mixed\n($L<0.95$)', 'locked\n($0.95\\leq L<1$)', 'strict\n($L=1$)'], fontsize=10)
ax2.set_xlim(-0.4, 2.4); ax2.set_ylim(0.3, 0.75)
ax2.set_ylabel('Within-target AUC of lag-1 TE')
ax2.set_xlabel("Target's locking index $L$")
ax2.set_title('(b) Graded by locking', loc='left', fontsize=11)
ax2.legend(loc='upper right', frameon=False, fontsize=9.5)
# (c) A11 and A15: synthetic networks
lams = [0.5, 0.7, 0.85, 0.95, 1.0]
for rho in (0.3, 0.6, 0.9):
    rs = [a11['conditions'][f'lambda={l},rho={rho}']['fast_lag1'] for l in lams]
    ax3.fill_between(lams, [r['lo'] for r in rs], [r['hi'] for r in rs], color=RAMP[rho], alpha=0.2, lw=0)
    ax3.plot(lams, [r['auc'] for r in rs], '-', marker=MARK[rho], color=RAMP[rho], lw=2, markersize=6.5,
             markeredgecolor='white', label=f'coupling {rho}')
ax3.axhline(0.5, color='0.55', ls='--', lw=1)
ax3.set_xlim(0.48, 1.03); ax3.set_ylim(0, 0.6)
ax3.set_xticks(lams); ax3.set_xticklabels(['0.5\n(none)', '0.7', '0.85', '0.95', '1.0\n(strict)'], fontsize=10)
ax3.set_xlabel('Share of events in preferred phase')
ax3.set_ylabel('AUC of lag-1 TE,\ndrivers vs non-drivers')
ax3.set_title('(c) Synthetic networks', loc='left', fontsize=11)
ax3.legend(loc='center right', bbox_to_anchor=(1.0, 0.68), frameon=False, fontsize=9.5)
for ax in axes:
    ax.grid(axis='y', alpha=0.3)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(str(OUT) + '.pdf'); fig.savefig(str(OUT) + '.png', dpi=150)
print('written', OUT)
