"""Thesis figure: what the lag-1 multivariate network (the main run's IDTxl network, 1,443 edges) would tell a
reader about the chip, beside the wiring. No new statistics: composition of its edges by phase relation, and
its ten highest-out-degree units against their wired out-degree, with the clock group for comparison.
Writes the figure and analyses/replication/output/published_network_summary.json."""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'src'))
from te6502.paths import DATA_ROOT, RESULTS, load_analysis  # noqa: E402
R = load_analysis('09_recomputing_the_supplementary_results/recompute_supplementary_results_pass1.py')
ROOT = DATA_ROOT

BLUE, ORANGE = '#2a78d6', '#eb6834'
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#d9d8d4'
plt.rcParams.update({'font.size': 8.5, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': INK,
                     'ytick.color': INK, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})

X, K, TE1, TE0, MI, p2g, NET = R.load(); G = X.shape[0]
ref = R.references(p2g, G); Wd, size, dist = ref['Wd'], ref['size'], ref['dist']
Lat = np.load(ROOT / 'analyses/new_tests/output/recording_reference.npz')['latch']
ph = R.phase_classes(X); rate = X.mean(1); clock = int(np.argmax(rate)); off = ~np.eye(G, dtype=bool)
E = NET & off; ii, jj = np.nonzero(E); wi, wj = np.nonzero(Wd & off)
def rel(a, b):
    same = (ph[a] == ph[b]) & (ph[a] < 2); opp = (ph[a] != ph[b]) & (ph[a] < 2) & (ph[b] < 2)
    return np.array([same.mean(), opp.mean(), 1 - same.mean() - opp.mean()]) * 100
comp = {'inferred': rel(ii, jj), 'wiring': rel(wi, wj)}
outd = E.sum(1); wout = (Wd & off).sum(1); top = np.argsort(-outd)[:10]
summary = dict(edges=int(len(ii)), wired=int(Wd[ii, jj].sum()), latch=int(Lat[ii, jj].sum()),
               composition_percent={k: dict(same=v[0], opposite=v[1], mixed=v[2]) for k, v in comp.items()},
               opposite_0_to_1=int(((ph[ii] == 0) & (ph[jj] == 1)).sum()), opposite_1_to_0=int(((ph[ii] == 1) & (ph[jj] == 0)).sum()),
               median_distance=dict(inferred=float(np.median(dist[ii, jj])), wiring=float(np.median(dist[wi, wj])), all_pairs=float(np.median(dist[off]))),
               top_sources=[dict(unit=int(u), inferred_out=int(outd[u]), wired_out=int(wout[u]), phase_class=int(ph[u])) for u in top],
               clock_group=dict(unit=clock, inferred_out=int(outd[clock]), inferred_in=int(E[:, clock].sum()), wired_out=int(wout[clock])),
               locked_pairs_opposite_share_percent=dict(inferred=float(100 * (((ph[ii] != ph[jj]) & (ph[ii] < 2) & (ph[jj] < 2)).sum() / ((ph[ii] < 2) & (ph[jj] < 2)).sum())),
                                                        wiring=float(100 * (((ph[wi] != ph[wj]) & (ph[wi] < 2) & (ph[wj] < 2)).sum() / ((ph[wi] < 2) & (ph[wj] < 2)).sum()))),
               top50_by_te=dict(wired=int(Wd[ii, jj][np.argsort(-np.nan_to_num(TE1[ii, jj], nan=-1))[:50]].sum()), latch=int(Lat[ii, jj][np.argsort(-np.nan_to_num(TE1[ii, jj], nan=-1))[:50]].sum()),
                                opposite=int((((ph[ii] != ph[jj]) & (ph[ii] < 2) & (ph[jj] < 2))[np.argsort(-np.nan_to_num(TE1[ii, jj], nan=-1))[:50]]).sum())),
               recall_percent=float(100 * Wd[ii, jj].sum() / Wd[off].sum()))
(RESULTS / 'figures/published_network_summary.json').write_text(json.dumps(summary, indent=1))

fig, (a, b) = plt.subplots(1, 2, figsize=(6.3, 3.1), gridspec_kw=dict(width_ratios=[1, 1.35]))
labs = ['Same\nphase', 'Opposite\nphases', 'Mixed']; x = np.arange(3); w = 0.36
for k, (name, col, dx) in enumerate((('Inferred lag-1 network', BLUE, -w / 2), ('The wiring', ORANGE, w / 2))):
    v = comp['inferred' if k == 0 else 'wiring']
    bars = a.bar(x + dx, v, w, color=col, label=name, edgecolor='white', lw=0.8)
    for xx, vv in zip(x + dx, v): a.text(xx, vv + 1.2, f'{vv:.0f}%', ha='center', va='bottom', fontsize=7.5, color=INK)
a.set_xticks(x); a.set_xticklabels(labs); a.set_ylabel('Share of pairs (%)'); a.set_ylim(0, 100)
a.set_title('a  Phase relation of the pairs', loc='left', fontsize=8.5); a.tick_params(axis='x', length=0)
a.legend(frameon=False, fontsize=7.5, loc='upper left', bbox_to_anchor=(0.0, 1.0), handlelength=1.2); a.grid(axis='y', color=GRID, lw=0.6); a.set_axisbelow(True)
units = [int(u) for u in top] + [clock]; y = np.arange(len(units)); h = 0.38
inf = [int(outd[u]) for u in units]; tru = [int(wout[u]) for u in units]
b.barh(y - h / 2, inf, h, color=BLUE, label='Inferred out-degree', edgecolor='white', lw=0.8)
b.barh(y + h / 2, tru, h, color=ORANGE, label='Wired out-degree', edgecolor='white', lw=0.8)
for yy, v1, v2 in zip(y, inf, tru):
    b.text(v1 + 3, yy - h / 2, str(v1), va='center', fontsize=7, color=INK); b.text(v2 + 3, yy + h / 2, str(v2), va='center', fontsize=7, color=INK)
b.set_yticks(y); b.set_yticklabels([f'unit {u}' for u in units[:-1]] + ['clock group']); b.invert_yaxis()
b.set_xlim(0, 320); b.set_xlabel('Out-degree'); b.set_title('b  The ten inferred hubs, and the real one', loc='left', fontsize=8.5)
b.axhline(len(units) - 1.5, color=GRID, lw=0.8); b.tick_params(axis='y', length=0); b.grid(axis='x', color=GRID, lw=0.6); b.set_axisbelow(True)
b.legend(frameon=False, fontsize=7.5, loc='upper right', handlelength=1.2)
fig.tight_layout()
OUTD = Path(sys.argv[1]) if len(sys.argv) > 1 else RESULTS / 'figures'; OUTD.mkdir(parents=True, exist_ok=True)
fig.savefig(OUTD / 'fig_c7_published_network.pdf'); fig.savefig(OUTD / 'fig_c7_published_network.png', dpi=200)
print(json.dumps(summary['composition_percent']), summary['clock_group'])
