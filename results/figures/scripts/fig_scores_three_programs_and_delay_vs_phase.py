"""Two thesis figures (Figures 5.1 and 5.7) from saved results in this repository (no new statistics).
Usage: python fig_scores_three_programs_and_delay_vs_phase.py [output_dir]   (default results/figures/)

(a) fig_c5_scores_three_programs: clock-aligned group AUC of each score on Donkey Kong, Pitfall and
    Space Invaders (R3 results; Donkey Kong from posthoc_align_offset.json and recomputed baselines).
(b) fig_c5_delay_vs_phase: R4, AUC at lags 0, 1, 2 for each phase relation and physical delay."""
import json, sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[3]
RES = REPO / 'results'
OUTD = Path(sys.argv[1]) if len(sys.argv) > 1 else RES / 'figures'
OUTD.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#d9d8d4'
plt.rcParams.update({'font.size': 8.5, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': INK, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})

def load(p): return json.loads(Path(p).read_text())

# ---------- (a)
R3 = {g: load(RES / f'10_remedy_held_out_R3/R3_{g}_results.json')['clock_aligned'] for g in ('Pitfall', 'SpaceInvaders')}
DKa = load(RES / '10_remedy_held_out_R3/posthoc_align_offset.json')['1408']
DKb = load(RES / '09_recomputed_supplementary/recomputed_supplementary.json')['B_baselines_all']
rows = [('Lag-1 transfer entropy', 'lag1', 'lag1_te', None),
        ('Lag-0 CMI', 'lag0', 'lag0_cmi', None),
        ('Zero-lag MI given the phase', 'mi_parity', 'mi_parity', None),
        ('Phase relation alone', 'phase_class', 'phase_class', None),
        ('Source-group size', None, 'source_size', 'source_group_size'),
        ('Layout distance', None, 'distance', 'layout_distance')]
def dk(r):
    _, a, _, b = r
    v = DKa[a] if a else DKb[b]
    return v['group_auc'], v['group_auc_ci']
def held(g, r):
    v = R3[g][f'group_auc_{r[2]}']; return v['auc'], [v['lo'], v['hi']]
progs = [('Donkey Kong', BLUE, 'o', lambda r: dk(r)), ('Pitfall', ORANGE, 's', lambda r: held('Pitfall', r)),
         ('Space Invaders', AQUA, '^', lambda r: held('SpaceInvaders', r))]
fig, ax = plt.subplots(figsize=(6.3, 3.0))
ys = [0, 1, 2, 3, 4.6, 5.6]
for k, (name, col, mk, get) in enumerate(progs):
    for y, r in zip(ys, rows):
        a, (lo, hi) = get(r); yy = y + (k - 1) * 0.22
        ax.plot([lo, hi], [yy, yy], color=col, lw=1.6, solid_capstyle='round', zorder=2)
        ax.plot(a, yy, mk, color=col, ms=5.5, mec='white', mew=0.8, zorder=3, label=name if y == 0 else None)
ax.axvline(0.5, color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
ax.text(0.506, 3.6, 'chance', color=MUTED, fontsize=7.5, va='bottom')
ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows]); ax.invert_yaxis()
ax.set_xlim(0.35, 1.0); ax.set_xlabel('Group-level AUC against the wiring (clock-aligned samples)')
ax.grid(axis='x', color=GRID, lw=0.6); ax.set_axisbelow(True); ax.tick_params(axis='y', length=0)
ax.axhline(3.8, color=GRID, lw=0.8)
ax.text(0.995, 3.6, 'Information measures', color=MUTED, fontsize=7.5, ha='right', va='bottom')
ax.text(0.995, 3.84, 'Using no activity', color=MUTED, fontsize=7.5, ha='right', va='top')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=3, frameon=False, fontsize=7.5, handletextpad=0.3)
fig.tight_layout(); fig.savefig(OUTD / 'fig_c5_scores_three_programs.pdf'); fig.savefig(OUTD / 'fig_c5_scores_three_programs.png', dpi=200)

# ---------- (b)
R4 = load(RES / '11_delay_vs_phase_R4/R4_delay_vs_phase_results.json')
cells = [('same_immediate', 'Same phase, immediate'), ('same_one_step', 'Same phase, one-step'),
         ('opposite_immediate', 'Opposite phases, immediate'), ('opposite_one_step', 'Opposite phases, one-step'),
         ('mixed_immediate', 'Mixed, immediate'), ('mixed_one_step', 'Mixed, one-step')]
lags = [('lag0_cmi', 'Lag 0', ORANGE, 's'), ('lag1_te', 'Lag 1', BLUE, 'o'), ('lag2_te', 'Lag 2', AQUA, '^')]
fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.9), sharey=True)
for ax, prog, title in zip(axs, ('program1', 'program2'), ('Program 1', 'Program 2')):
    for y, (c, lab) in enumerate(cells):
        for k, (key, name, col, mk) in enumerate(lags):
            v = R4[prog][c][key]; yy = y + (k - 1) * 0.22
            ax.plot([v['lo'], v['hi']], [yy, yy], color=col, lw=1.6, solid_capstyle='round', zorder=2)
            ax.plot(v['auc'], yy, mk, color=col, ms=5, mec='white', mew=0.8, zorder=3, label=name if y == 0 else None)
    ax.axvline(0.5, color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.set_xlim(0.38, 0.70); ax.set_title(title, fontsize=8.5, color=INK, loc='left')
    ax.grid(axis='x', color=GRID, lw=0.6); ax.set_axisbelow(True); ax.tick_params(axis='y', length=0)
    for yb in (1.5, 3.5): ax.axhline(yb, color=GRID, lw=0.8)
axs[0].set_yticks(range(len(cells))); axs[0].set_yticklabels([c[1] for c in cells]); axs[0].invert_yaxis()
fig.supxlabel('AUC against pairs with no influence (stratified by event rate)', fontsize=8.5, color=INK)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc='upper center', bbox_to_anchor=(0.62, 1.0), ncol=3, frameon=False, fontsize=7.5, handletextpad=0.3)
fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(OUTD / 'fig_c5_delay_vs_phase.pdf'); fig.savefig(OUTD / 'fig_c5_delay_vs_phase.png', dpi=200)
print('written')
