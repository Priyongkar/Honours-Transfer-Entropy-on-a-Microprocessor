"""A8: lag-1 TE on the synthetic NMOS circuits' known latches (design: PREDICTIONS.md).
Usage: python nmos.py [--quick]"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / 'src'))  # repository library
from te6502.paths import ANALYSES, DATA_ROOT, RESULTS, HELPERS_FILE, load_analysis  # noqa: E402,F401
import argparse
import numpy as np
from te6502.test_helpers import tk, log, start, save, change, phase_class, auc_with_ci, verdict, quintile
from te6502 import switch_level_simulator as sim6502

SEEDS = range(500, 510)


def circuit(seed, T):
    defs, meta = sim6502.make_nmos_circuit(n_gates=80, n_inputs=6, n_latches=16, seed=seed)
    S, chip = sim6502.simulate_nmos(defs, meta, T + 200, seed=seed + 13)
    S = S[200:].T.astype(np.uint8)
    nodes = np.asarray(chip.used_nodes); idx = {int(v): i for i, v in enumerate(nodes)}
    g_, c1_, c2_ = chip.transistor_matrix()
    A = tk.node_reference_from_transistors(g_, c1_, c2_, nodes, supply={chip.gnd, chip.pwr}).toarray() > 0
    inp, rdr = set(), set()
    for s_, src, _ in meta['latches']:
        if s_ in idx and src in idx:
            inp.add((idx[src], idx[s_]))
        for g in meta['gates']:
            if s_ in g['inputs'] and s_ in idx and g['output'] in idx:
                rdr.add((idx[s_], idx[g['output']]))
    return S, A, inp, rdr


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); a = ap.parse_args()
    sha = start('nmos')
    lengths = [2000] if a.quick else [6000, 20000]
    seeds = list(SEEDS)[:3] if a.quick else list(SEEDS)
    res = dict(predictions_sha256=sha, quick=a.quick)
    for T in lengths:
        cols = {k: [] for k in ['lag1', 'lag0', 'inp', 'rdr', 'neg', 'strat', 'circ']}
        for c, seed in enumerate(seeds):
            S, A, inp, rdr = circuit(seed, T)
            D = change(S)
            reps, ug = tk.group_identical(D); Xg = D[reps]
            ks = np.asarray([tk.ais_select_k(x, 4) for x in Xg])
            L1 = tk.te_matrix(Xg, ks, delay=1)[0]; L0 = tk.te_matrix(Xg, ks, delay=0)[0]
            n = len(S); ph, _ = phase_class(D); rq = quintile(D.mean(1))
            I = np.zeros((n, n), bool); R = np.zeros((n, n), bool)
            for p in inp: I[p] = True
            for p in rdr: R[p] = True
            rel = A | A.T | I | I.T | R | R.T
            valid = (ug[:, None] >= 0) & (ug[None, :] >= 0) & (ug[:, None] != ug[None, :])
            ii, jj = np.nonzero(valid)
            cols['lag1'].append(L1[ug[ii], ug[jj]]); cols['lag0'].append(L0[ug[ii], ug[jj]])
            cols['inp'].append(I[ii, jj]); cols['rdr'].append(R[ii, jj]); cols['neg'].append(~rel[ii, jj])
            cols['strat'].append((((c * 3 + ph[ii]) * 3 + ph[jj]) * 5 + rq[ii]) * 5 + rq[jj])
            cols['circ'].append(np.full(len(ii), c))
            log(f'T={T} circuit {seed}: groups {len(reps)}, input->storage {I[ii, jj].sum()}, storage->reader {R[ii, jj].sum()}')
        C = {k: np.concatenate(v) for k, v in cols.items()}
        mult = np.random.default_rng(2026).multinomial(len(seeds), np.full(len(seeds), 1 / len(seeds)), size=50 if a.quick else 1000)
        out = {}
        for pname in ['inp', 'rdr']:
            for sname in ['lag1', 'lag0']:
                m = C[pname] | C['neg']
                r, bs = auc_with_ci(C[sname][m], C[pname][m], C['strat'][m], C['circ'][m], mult)
                r['verdict'] = verdict(r); out[f'{pname}_{sname}'] = r
        res[f'T{T}'] = out
        log(f"T={T}: input->storage lag1 {out['inp_lag1']['auc']:.3f}, storage->reader lag1 {out['rdr_lag1']['auc']:.3f}")
    save('nmos_results' + ('_quick' if a.quick else ''), res)
    log('done')


if __name__ == '__main__':
    main()
