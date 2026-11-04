"""Figure F5: lag profile of oriented latch pairs and pull-down wired pairs on the
simulated chip, both programs (A5, A7-A5). Reads output/chip_results.json and
output/A7_results.json; writes thesis/images/new/fig_c5_latch_profile.{pdf,png}."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / 'thesis/images/new/fig_c5_latch_profile'
BLUE, ORANGE = '#2b6fd6', '#e8692d'

fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6), sharey=True)
for ax, (f, title) in zip(axes, [('chip_results.json', '(a) Program 1 (Fibonacci and LFSR)'),
                                 ('A7_results.json', '(b) Program 2 (held out)')]):
    prof = json.load(open(HERE / 'output' / f))['A5_chip_delay_profile']
    u = [p['delay'] for p in prof]
    for key, colour, label in [('latch', BLUE, 'latch pairs, driven to storage'),
                               ('pulldown', ORANGE, 'pull-down wired pairs')]:
        auc = [p[key]['auc'] for p in prof]
        lo = [p[key]['lo'] for p in prof]
        hi = [p[key]['hi'] for p in prof]
        ax.fill_between(u, lo, hi, color=colour, alpha=0.18, linewidth=0)
        ax.plot(u, auc, 'o-', color=colour, label=label, markersize=5, linewidth=1.5)
    ax.axhline(0.5, color='0.4', linestyle='--', linewidth=1)
    ax.set_title(title, loc='left', fontsize=11)
    ax.set_xlabel('Lag $u$ of the source (half-steps)')
    ax.set_xticks(u)
    ax.set_ylim(0, 1.02)
    ax.grid(alpha=0.3)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
axes[0].set_ylabel('Stratified AUC against\npairs with no measured influence')
axes[0].text(8.0, 0.52, 'chance', ha='right', va='bottom', fontsize=9, color='0.4')
axes[1].legend(loc='upper right', frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(str(OUT) + '.pdf')
fig.savefig(str(OUT) + '.png', dpi=150)
print('written', OUT)
