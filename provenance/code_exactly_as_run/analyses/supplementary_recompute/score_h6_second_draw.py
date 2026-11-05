"""R5 scoring: the second independent H6 draw, scored as the run of record (h6_run_of_record_R1.py), plus
the source-size baseline at the same allocation and the two draws pooled (design: ../replication/PREDICTIONS_R5.md).
Writes output/h6_second_draw_R5.json."""
import json, sys, glob
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import recompute_supplementary as R                                    # noqa: E402
tk = R.tk
NT = HERE.parents[0] / 'new_tests/output'; RP = HERE.parents[0] / 'replication/output'


def network(folder, G, Wd):
    summ = [json.load(open(f)) for f in sorted(glob.glob(str(folder / 'target_*.json')))]
    for s in summ:
        s['selected_vars_sources'] = [v for v in s['selected_vars_sources'] if v[0] < G]
    keep = tk.network_fdr_omnibus(summ)
    A = np.zeros((G, G), bool); lag_w = {0: [0, 0], 1: [0, 0]}
    for s in summ:
        if keep.get(s['target']):
            for u, lag in s['selected_vars_sources']:
                A[u, s['target']] = True; lag_w[int(lag)][0] += 1; lag_w[int(lag)][1] += int(Wd[u, s['target']])
    return A, dict(targets=len(summ), targets_kept=int(sum(keep.values())),
                   by_lag={f'lag{l}': dict(sources=v[0], wired=v[1]) for l, v in lag_w.items()})


def allocate(score, A, T, G):
    alloc = np.zeros(G, int); alloc[T] = A[:, T].sum(0); B = np.zeros((G, G), bool)
    for t in T:
        col = np.where(np.isfinite(score[:, t]), score[:, t], -np.inf); col[t] = -np.inf
        B[np.argsort(-col, kind='stable')[:alloc[t]], t] = True
    return B


def score_set(name, A, T, G, Wd, Wphys, Nphys, TE0, size, NET, seed):
    M = np.zeros((G, G), bool); M[:, T] = NET[:, T]
    L0 = allocate(TE0, A, T, G); SZ = allocate(np.repeat(size[:, None], G, 1).astype(float), A, T, G)
    mt = np.random.default_rng(seed).multinomial(len(T), np.full(len(T), 1 / len(T)), size=R.DRAWS)
    ms = np.random.default_rng(seed + 1).multinomial(G, np.full(G, 1 / G), size=R.DRAWS)
    def stats(B, wt=None, ws=None):
        wt = np.ones(len(T)) if wt is None else wt; ws = np.ones(G) if ws is None else ws
        Bt = B[:, T]; W = ws[:, None] * wt[None, :]
        e = (Bt * W).sum(); h = ((Bt & Wd[:, T]) * W).sum()
        pe = ((Wphys[:, T] + Nphys[:, T]) * Bt * W).sum(); ph = (Wphys[:, T] * Bt * W).sum()
        return dict(group_precision=100 * h / e if e else np.nan, group_recall=100 * h / max((Wd[:, T] * W).sum(), 1e-9),
                    physical_precision=100 * ph / pe if pe else np.nan)
    res = {}
    for nm, B in (('search_lags_0_1', A), ('main_run_lag1_same_targets', M), ('lag0_cmi_same_allocation', L0), ('source_size_same_allocation', SZ)):
        r = stats(B); r['edges'] = int(B[:, T].sum()); r['wired'] = int((B & Wd)[:, T].sum())
        r['group_precision_ci'], _ = R.ci(lambda m: stats(B, wt=m.astype(float))['group_precision'], mt); res[nm] = r
    diffs = {}
    for nm, (a, b) in (('search_minus_lag1', (A, M)), ('search_minus_lag0_cmi', (A, L0)), ('size_minus_search', (SZ, A))):
        dd = {}
        for k in ('group_precision', 'physical_precision'):
            d0 = stats(a)[k] - stats(b)[k]
            b1 = np.array([stats(a, wt=m.astype(float))[k] - stats(b, wt=m.astype(float))[k] for m in mt])
            b2 = np.array([stats(a, wt=mt[i].astype(float), ws=ms[i].astype(float))[k] - stats(b, wt=mt[i].astype(float), ws=ms[i].astype(float))[k] for i in range(R.DRAWS)])
            b1, b2 = b1[np.isfinite(b1)], b2[np.isfinite(b2)]
            dd[k] = dict(diff=d0, target_ci=[float(np.quantile(b1, .025)), float(np.quantile(b1, .975))],
                         two_way_ci=[float(np.quantile(b2, .025)), float(np.quantile(b2, .975))])
        diffs[nm] = dd
    return dict(networks=res, differences=diffs)


def main():
    X, K, TE1, TE0, MI, p2g, NET = R.load(); G = X.shape[0]
    ref = R.references(p2g, G); Wd, Wphys, Nphys, size = ref['Wd'], ref['Wphys'], ref['Nphys'], ref['size']
    T5 = np.array(json.load(open(RP / 'R5_targets.json'))['targets']); T1 = np.array(json.load(open(NT / 'idtxl_results.json'))['targets'])
    assert not set(T5) & set(T1)
    A5, meta5 = network(NT / 'R5_targets', G, Wd); A1, meta1 = network(NT / 'R1_targets', G, Wd)
    out = dict(R5_meta=meta5, R1_meta=meta1, R5=score_set('R5', A5, T5, G, Wd, Wphys, Nphys, TE0, size, NET, R.SEED),
               R1=score_set('R1', A1, T1, G, Wd, Wphys, Nphys, TE0, size, NET, R.SEED),
               pooled=score_set('pooled', A5 | A1, np.concatenate([T1, T5]), G, Wd, Wphys, Nphys, TE0, size, NET, R.SEED))
    r5 = out['R5']; lag0 = meta5['by_lag']['lag0']['sources']; tot = lag0 + meta5['by_lag']['lag1']['sources']
    out['predictions_met'] = {'R5.1': r5['differences']['search_minus_lag1']['group_precision']['target_ci'][0] > 0,
                              'R5.2': lag0 / max(tot, 1) > 0.5,
                              'R5.3': r5['differences']['size_minus_search']['group_precision']['target_ci'][0] > 0,
                              'R5.4_descriptive_includes_0': r5['differences']['search_minus_lag0_cmi']['group_precision']['target_ci'][0] <= 0 <= r5['differences']['search_minus_lag0_cmi']['group_precision']['target_ci'][1]}
    (HERE / 'output/h6_second_draw_R5.json').write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps({k: out[k] for k in ('R5_meta', 'predictions_met')}, indent=1, default=float))
    for nm in ('R5', 'R1', 'pooled'):
        n = out[nm]['networks']; print(nm, {k: round(v['group_precision'], 2) for k, v in n.items()})


if __name__ == '__main__':
    main()
