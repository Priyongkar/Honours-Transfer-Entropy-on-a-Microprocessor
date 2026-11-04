"""Headline figure: the followership penalty (A9) and its generality (A11).
Reads output/A9_results.json and output/A11_results.json;
writes thesis/images/new/fig_c6_followership.{pdf,png}."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / 'thesis/images/new/fig_c6_followership'
ORANGE, INDIGO = '#e8692d', '#4a5fc1'                    # lag-1 TE is orange throughout the thesis
RAMP = {0.3: '#f0955f', 0.6: '#d4561c', 0.9: '#8f3509'}   # coupling strength, light to dark
MARK = {0.3: 'o', 0.6: 's', 0.9: 'D'}
INK, MUTED = '#222222', '#666666'

a9 = json.load(open(HERE / 'output/A9_results.json'))
a11 = json.load(open(HERE / 'output/A11_results.json'))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.0, 3.9), gridspec_kw=dict(width_ratios=[1, 1.25]))

# (a) A9: within-target AUC of activity and of lag-1 TE, wired against unwired sources
sets = [('recording', 'Recording'), ('chip_program1', 'Chip,\nprogram 1'), ('chip_program2', 'Chip,\nprogram 2')]
for x, (key, _) in enumerate(sets):
    t = a9[key]['tests']
    for name, colour, dx, marker in [('desc_activity_within_target', INDIGO, -0.12, 'o'),
                                     ('desc_TE_within_target', ORANGE, 0.12, 's')]:
        r = t[name]
        ax1.errorbar(x + dx, r['auc'], yerr=[[r['auc'] - r['lo']], [r['hi'] - r['auc']]],
                     fmt=marker, color=colour, markersize=8, capsize=3, elinewidth=1.5, markeredgecolor='white')
    ax1.annotate('', xy=(x + 0.12, t['desc_TE_within_target']['auc']),
                 xytext=(x - 0.12, t['desc_activity_within_target']['auc']),
                 arrowprops=dict(arrowstyle='->', color=MUTED, lw=1))
ax1.axhline(0.5, color='0.55', linestyle='--', linewidth=1)
ax1.set_xticks(range(3)); ax1.set_xticklabels([s[1] for s in sets], fontsize=9)
ax1.set_xlim(-0.5, 2.5); ax1.set_ylim(0.3, 0.65)
ax1.set_ylabel('Within-target AUC,\nwired against unwired sources', color=INK)
ax1.set_title('(a) Being followed lowers wired sources', loc='left', fontsize=11)
ax1.text(-0.47, 0.503, 'chance', ha='left', va='bottom', fontsize=8.5, color=MUTED)
ax1.plot([], [], 'o', color=INDIGO, markersize=7, label='activity alone')
ax1.plot([], [], 's', color=ORANGE, markersize=7, label='lag-1 TE')
ax1.legend(loc='upper right', frameon=False, fontsize=8.5)

# (b) A11: lag-1 TE for within-sample drivers against locking strength, one line per coupling
lams = [0.5, 0.7, 0.85, 0.95, 1.0]
for rho in [0.3, 0.6, 0.9]:
    rs = [a11['conditions'][f'lambda={l},rho={rho}']['fast_lag1'] for l in lams]
    auc = [r['auc'] for r in rs]; lo = [r['lo'] for r in rs]; hi = [r['hi'] for r in rs]
    ax2.fill_between(lams, lo, hi, color=RAMP[rho], alpha=0.2, linewidth=0)
    ax2.plot(lams, auc, '-', marker=MARK[rho], color=RAMP[rho], linewidth=2, markersize=7,
             markeredgecolor='white', label=f'coupling {rho}')
    ax2.text(1.012, auc[-1], f'{rho}', va='center', fontsize=8.5, color=INK)
ax2.axhline(0.5, color='0.55', linestyle='--', linewidth=1)
ax2.text(0.99, 0.508, 'chance', ha='right', va='bottom', fontsize=8.5, color=MUTED)
ax2.set_xlim(0.48, 1.06); ax2.set_ylim(0, 0.6)
ax2.set_xticks(lams); ax2.set_xticklabels(['0.5\n(none)', '0.7', '0.85', '0.95', '1.0\n(strict)'], fontsize=9)
ax2.set_xlabel('Locking strength: share of events in the preferred phase', color=INK)
ax2.set_ylabel('AUC of lag-1 TE, within-sample\ndrivers against non-drivers', color=INK)
ax2.set_title('(b) Synthetic networks: the penalty grows with locking', loc='left', fontsize=11)
ax2.legend(loc='center right', bbox_to_anchor=(1.0, 0.72), frameon=False, fontsize=8.5)

for ax in (ax1, ax2):
    ax.grid(axis='y', alpha=0.3)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK)
fig.tight_layout()
fig.savefig(str(OUT) + '.pdf')
fig.savefig(str(OUT) + '.png', dpi=150)
print('written', OUT)
