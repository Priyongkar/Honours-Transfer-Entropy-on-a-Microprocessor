"""Figure (Ch5, replication): the four headline rankings on the three recorded programs. Donkey Kong values from
the main analyses (lag-1 TE and lag-0 CMI clock-aligned, Table 5.1; latch pairs and layout distance whole recording,
which for distance is activity-free); Pitfall and Space Invaders from the R2 outputs (clock-aligned).
Writes thesis/images/new/fig_c5_replication.{pdf,png}."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os; H = Path(os.environ.get('HONOURS_DATA_ROOT', Path.home() / 'Documents/Honours 2')) / 'analyses'
OUT = H.parent / 'thesis/images/new/fig_c5_replication'
rp = {g: json.load(open(H / f'replication/output/R2_{f}_results.json'))['clock_aligned'] for g, f in (('Pitfall', 'Pitfall'), ('Space Invaders', 'SpaceInvaders'))}
def ci(d): return (d['auc'], d['lo'], d['hi'])
rows = ['Lag-1 TE\nagainst wiring', 'Lag-0 CMI\nagainst wiring', 'Lag-1 TE\non latch pairs', 'Layout distance\nagainst wiring']
data = {
 'Donkey Kong': [(0.396, 0.372, 0.421), (0.561, 0.538, 0.582), (0.985, 0.975, 0.993), (0.928, 0.915, 0.940)],
 'Pitfall': [ci(rp['Pitfall']['R2.1_group_auc_lag1']), ci(rp['Pitfall']['group_auc_lag0']), ci(rp['Pitfall']['R2.6_latch_lag1']), ci(rp['Pitfall']['group_auc_distance'])],
 'Space Invaders': [ci(rp['Space Invaders']['R2.1_group_auc_lag1']), ci(rp['Space Invaders']['group_auc_lag0']), ci(rp['Space Invaders']['R2.6_latch_lag1']), ci(rp['Space Invaders']['group_auc_distance'])],
}
col = {'Donkey Kong': '#222222', 'Pitfall': '#4a5fc1', 'Space Invaders': '#2a9d8f'}
mk = {'Donkey Kong': 'o', 'Pitfall': 's', 'Space Invaders': 'D'}
plt.rcParams.update({'font.size': 11})
fig, ax = plt.subplots(figsize=(9.0, 4.2))
for j, (g, vals) in enumerate(data.items()):
    for i, (a, lo, hi) in enumerate(vals):
        x = i + (j - 1) * 0.18
        ax.errorbar(x, a, yerr=[[a - lo], [hi - a]], fmt=mk[g], color=col[g], markersize=8, capsize=3,
                    elinewidth=1.4, markeredgecolor='white', label=g if i == 0 else None)
ax.axhline(0.5, color='0.55', ls='--', lw=1); ax.text(-0.35, 0.505, 'chance', ha='left', va='bottom', fontsize=9.5, color='#666666')
ax.set_xticks(range(4)); ax.set_xticklabels(rows, fontsize=10.5)
ax.set_ylim(0.3, 1.02); ax.set_ylabel('AUC (95% interval)')
ax.legend(loc='center right', bbox_to_anchor=(1.0, 0.62), frameon=False, fontsize=10)
ax.grid(axis='y', alpha=0.3)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(str(OUT) + '.pdf'); fig.savefig(str(OUT) + '.png', dpi=150); print('written', OUT)
