"""H6 on its run of record (R1): the stratified 100-target IDTxl run with source lags 0 and 1, whose target
list was recorded. Recomputes the H6 table (edges, wired, group and physical precision and recall with
target-cluster intervals; the paired difference against the main run's lag-1 network on the same targets,
under the target and two-way bootstraps; composition by lag; lag-0 CMI allocated the same way) from the
saved per-target IDTxl outputs. Supplementary: the original H6 run's target list was not kept.
Writes output/h6_run_of_record_R1.json."""
import json, sys, glob
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'src'))
from te6502.paths import load_analysis, RESULTS, DATA_ROOT  # noqa: E402
R = load_analysis('09_recomputing_the_supplementary_results/recompute_supplementary_results_pass1.py')                                    # noqa: E402
tk = R.tk
NT = DATA_ROOT / 'analyses/new_tests/output'


def main():
    X, K, TE1, TE0, MI, p2g, NET = R.load(); G = X.shape[0]
    ref = R.references(p2g, G); Wd, Wphys, Nphys = ref['Wd'], ref['Wphys'], ref['Nphys']
    targets = json.load(open(NT / 'idtxl_results.json'))['targets']; T = np.array(targets)
    out = {}
    for name, d in (('R1', 'R1_targets'), ('A17', 'A17_targets')):
        summ = [json.load(open(f)) for f in sorted(glob.glob(str(NT / d / '*.json')))]
        for s in summ:
            s['selected_vars_sources'] = [v for v in s['selected_vars_sources'] if v[0] < G]
        keep = tk.network_fdr_omnibus(summ)
        A = np.zeros((G, G), bool); n_sel = 0; lag_w = {0: [0, 0], 1: [0, 0]}
        for s in summ:
            if keep.get(s['target']):
                for u, lag in s['selected_vars_sources']:
                    A[u, s['target']] = True; n_sel += 1
                    lag_w[int(lag)][0] += 1; lag_w[int(lag)][1] += int(Wd[u, s['target']])
        out[name] = dict(selected_sources=n_sel, edges=int(A.sum()), wired=int((A & Wd).sum()),
                         by_lag={f'lag{l}': dict(sources=v[0], wired=v[1]) for l, v in lag_w.items()}, A=A)
    A = out['R1']['A']
    M = np.zeros((G, G), bool); M[:, T] = NET[:, T]                    # main run's lag-1 network, same targets
    alloc = np.zeros(G, int); alloc[T] = A[:, T].sum(0)
    L0 = np.zeros((G, G), bool)                                         # lag-0 CMI, same allocation
    for t in T:
        col = np.where(np.isfinite(TE0[:, t]), TE0[:, t], -np.inf); col[t] = -np.inf
        L0[np.argsort(-col, kind='stable')[:alloc[t]], t] = True
    rng = np.random.default_rng(R.SEED); mt = rng.multinomial(len(T), np.full(len(T), 1 / len(T)), size=R.DRAWS)
    ms = np.random.default_rng(R.SEED + 1).multinomial(G, np.full(G, 1 / G), size=R.DRAWS)
    wired_t = Wd[:, T].sum(0); wphys_t = Wphys[:, T].sum(0)

    def stats(B, wt=None, ws=None):
        wt = np.ones(len(T)) if wt is None else wt; ws = np.ones(G) if ws is None else ws
        Bt = B[:, T]; W = ws[:, None] * wt[None, :]
        e = (Bt * W).sum(); h = ((Bt & Wd[:, T]) * W).sum()
        pe = ((Wphys[:, T] + Nphys[:, T]) * Bt * W).sum(); ph = (Wphys[:, T] * Bt * W).sum()
        return dict(group_precision=100 * h / e if e else np.nan, group_recall=100 * h / (Wd[:, T] * W).sum(),
                    physical_precision=100 * ph / pe if pe else np.nan, physical_recall=100 * ph / (Wphys[:, T] * W).sum())
    res = {}
    for nm, B in (('R1_lags_0_1', A), ('main_run_lag1_same_targets', M), ('lag0_cmi_same_allocation', L0)):
        r = stats(B); r['edges'] = int(B[:, T].sum()); r['wired'] = int((B & Wd)[:, T].sum())
        for k in ('group_precision', 'physical_precision'):
            r[k + '_ci'], _ = R.ci(lambda m: stats(B, wt=m.astype(float))[k], mt)
        res[nm] = r
    diffs = {}
    for nm, (a, b) in (('R1_minus_lag1', (A, M)), ('R1_minus_lag0_cmi', (A, L0))):
        dd = {}
        for k in ('group_precision', 'physical_precision'):
            d0 = stats(a)[k] - stats(b)[k]
            b1 = np.array([stats(a, wt=m.astype(float))[k] - stats(b, wt=m.astype(float))[k] for m in mt])
            b2 = np.array([stats(a, wt=mt[i].astype(float), ws=ms[i].astype(float))[k] - stats(b, wt=mt[i].astype(float), ws=ms[i].astype(float))[k] for i in range(R.DRAWS)])
            b1, b2 = b1[np.isfinite(b1)], b2[np.isfinite(b2)]
            dd[k] = dict(diff=d0, target_ci=[float(np.quantile(b1, .025)), float(np.quantile(b1, .975))], target_positive=float((b1 > 0).mean()),
                         two_way_ci=[float(np.quantile(b2, .025)), float(np.quantile(b2, .975))], two_way_positive=float((b2 > 0).mean()))
        diffs[nm] = dd
    for nm in ('R1', 'A17'):
        out[nm].pop('A')
    out['networks'] = res; out['differences'] = diffs
    (RESULTS / '09_recomputed_supplementary/h6_run_of_record_R1.json').write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
