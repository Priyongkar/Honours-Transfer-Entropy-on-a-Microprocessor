"""Details of the H6 run of record (R1) for the thesis table: precision of the lag-0 and lag-1 selections,
allocation-matched lag-1 TE, source-group size and layout distance, edge overlap with lag-0 CMI, and the
transistor-path relation classes of each network. Supplementary. Writes output/h6_run_of_record_R1_details.json."""
import json, sys, glob
from pathlib import Path
import numpy as np
from scipy import sparse
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import recompute_supplementary as R                                    # noqa: E402
tk = R.tk
NT = HERE.parents[0] / 'new_tests/output'


def main():
    X, K, TE1, TE0, MI, p2g, NET = R.load(); G = X.shape[0]
    ref = R.references(p2g, G); Wd, Wphys, Nphys, size, dist = ref['Wd'], ref['Wphys'], ref['Nphys'], ref['size'], ref['dist']
    T = np.array(json.load(open(NT / 'idtxl_results.json'))['targets'])
    summ = [json.load(open(f)) for f in sorted(glob.glob(str(NT / 'R1_targets/*.json')))]
    for s in summ:
        s['selected_vars_sources'] = [v for v in s['selected_vars_sources'] if v[0] < G]
    keep = tk.network_fdr_omnibus(summ)
    A = np.zeros((G, G), bool); Al = {0: np.zeros((G, G), bool), 1: np.zeros((G, G), bool)}
    for s in summ:
        if keep.get(s['target']):
            for u, lag in s['selected_vars_sources']:
                A[u, s['target']] = True; Al[int(lag)][u, s['target']] = True
    alloc = np.zeros(G, int); alloc[T] = A[:, T].sum(0)

    def matched(S):
        B = np.zeros((G, G), bool)
        for t in T:
            col = np.where(np.isfinite(S[:, t]), S[:, t], -np.inf); col[t] = -np.inf
            B[np.argsort(-col, kind='stable')[:alloc[t]], t] = True
        return B
    sz = np.repeat(size[:, None].astype(float), G, 1)
    nets = {'R1_all': A, 'R1_lag0': Al[0], 'R1_lag1': Al[1], 'lag0_cmi': matched(TE0), 'lag1_te': matched(TE1),
            'source_group_size': matched(sz), 'layout_distance': matched(-dist.astype(float))}
    mt = np.random.default_rng(R.SEED).multinomial(len(T), np.full(len(T), 1 / len(T)), size=R.DRAWS)

    def st(B, wt=None):
        wt = np.ones(len(T)) if wt is None else wt; Bt = B[:, T] * wt[None, :]
        e = Bt.sum(); h = (Bt * Wd[:, T]).sum(); pe = (Bt * (Wphys[:, T] + Nphys[:, T])).sum(); ph = (Bt * Wphys[:, T]).sum()
        return 100 * h / e if e else np.nan, 100 * ph / pe if pe else np.nan, 100 * h / (Wd[:, T] * wt[None, :]).sum(), 100 * ph / (Wphys[:, T] * wt[None, :]).sum()
    # transistor-path classes
    gate, c1, c2, E = ref['gate'], ref['c1'], ref['c2'], ref['E']; nT = len(gate)
    At = sparse.csr_matrix((np.ones(len(E)), (E[:, 0], E[:, 1])), shape=(nT, nT)); act = p2g >= 0
    P = sparse.csr_matrix((np.ones(act.sum()), (np.flatnonzero(act), p2g[act])), shape=(nT, G))
    def grp(M):
        Q = (P.T @ (M > 0).astype(float) @ P).toarray() > 0; np.fill_diagonal(Q, False); return Q
    w, r, two, sp, sc = grp(At), grp(At.T), grp(At @ At), grp(At.T @ At), grp(At @ At.T)
    C = np.full((G, G), 5)
    for k, M in ((4, sc), (3, sp), (2, two), (1, r & ~w), (0, w)):
        C[M] = k
    names = ['wired', 'reverse_only', 'two_step', 'shared_parent', 'shared_child', 'unrelated']
    out = {}
    for nm, B in nets.items():
        g, ph, gr, pr = st(B)
        cg, _ = R.ci(lambda m: st(B, m.astype(float))[0], mt); cp, _ = R.ci(lambda m: st(B, m.astype(float))[1], mt)
        cls = np.bincount(C[B], minlength=6)
        out[nm] = dict(edges=int(B[:, T].sum()), wired_group=int((B & Wd)[:, T].sum()), group_precision=g, group_ci=cg,
                       group_recall=gr, physical_precision=ph, physical_ci=cp, physical_recall=pr,
                       transistor_path_classes={n: int(c) for n, c in zip(names, cls)})
    L0 = nets['lag0_cmi']
    out['overlap_R1_lag0cmi'] = dict(shared_edges=int((A & L0).sum()), wired_only_in_R1=int((A & ~L0 & Wd).sum()),
                                     wired_only_in_lag0cmi=int((L0 & ~A & Wd).sum()))
    (HERE / 'output/h6_run_of_record_R1_details.json').write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
