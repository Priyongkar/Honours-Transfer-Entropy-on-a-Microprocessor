"""Thesis figure: H6 on two independent samples of 100 targets (R1, the run of record; R5) and pooled.
Group precision (95% target-cluster intervals) of the main run's lag-1 network on the same targets, the
lags-{0,1} multivariate search, lag-0 CMI and source-group size allocated as the search allocated.
Reads analyses/supplementary_recompute/output/h6_second_draw_R5.json; no new statistics."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[2]
d = json.loads((ROOT / 'analyses/supplementary_recompute/output/h6_second_draw_R5.json').read_text())
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#d9d8d4'
plt.rcParams.update({'font.size': 8.5, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': INK, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
rows = [('Main run, lag-1 sources only', 'main_run_lag1_same_targets'), ('Search admitting lag-0 sources', 'search_lags_0_1'),
        ('Lag-0 CMI, same allocation', 'lag0_cmi_same_allocation'), ('Source-group size, same allocation', 'source_size_same_allocation')]
sets = [('First sample (R1)', 'R1', BLUE, 'o'), ('Second sample (R5)', 'R5', ORANGE, 's'), ('Pooled', 'pooled', AQUA, '^')]
fig, ax = plt.subplots(figsize=(6.3, 2.5))
for k, (lab, key, col, mk) in enumerate(sets):
    for y, (_, r) in enumerate(rows):
        v = d[key]['networks'][r]; lo, hi = v['group_precision_ci']; yy = y + (k - 1) * 0.22
        ax.plot([lo, hi], [yy, yy], color=col, lw=1.6, solid_capstyle='round', zorder=2)
        ax.plot(v['group_precision'], yy, mk, color=col, ms=5.5, mec='white', mew=0.8, zorder=3, label=lab if y == 0 else None)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows]); ax.invert_yaxis()
ax.set_xlim(0, 37); ax.set_xlabel('Group-level precision against the wiring (%)')
ax.grid(axis='x', color=GRID, lw=0.6); ax.set_axisbelow(True); ax.tick_params(axis='y', length=0)
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=3, frameon=False, fontsize=7.5, handletextpad=0.3)
fig.tight_layout()
OUTD = ROOT / 'thesis/images/new'
fig.savefig(OUTD / 'fig_c6_h6_two_samples.pdf'); fig.savefig(OUTD / 'fig_c6_h6_two_samples.png', dpi=200)
print('written')
