"""A10: A8 replicated on 50 new NMOS circuits, seeds 600-649 (design: PREDICTIONS_A9.md, section A10).
Circuit generation and scoring are those of nmos.py (A8), imported unchanged.

Usage:
  python a10_nmos.py --selftest   times one circuit (seed 9999, outside the test) at both lengths; prints seconds only
  python a10_nmos.py              full run
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / 'src'))  # repository library
from te6502.paths import ANALYSES, DATA_ROOT, RESULTS, HELPERS_FILE, load_analysis  # noqa: E402,F401
import argparse, hashlib, time
from pathlib import Path
import numpy as np
from te6502.test_helpers import HERE, tk, log, save, change, phase_class, auc_with_ci, verdict, diff_ci, quintile
circuit = load_analysis('02_latch_transfer_tests_A1_to_A8/nmos_circuit_latch_tests_A8.py').circuit

PRED = ANALYSES / '04_followership_mechanism_A9_to_A11/design/predictions_A9_A10.md'
NEW = list(range(600, 650))
OLD = list(range(500, 510))          # A8's circuits, re-simulated for the pooled, descriptive total
LENGTHS = [6000, 20000]
SEED, DRAWS = 2026, 1000


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def score_circuit(seed, T):
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
    return dict(lag1=L1[ug[ii], ug[jj]], lag0=L0[ug[ii], ug[jj]], inp=I[ii, jj], rdr=R[ii, jj],
                neg=~rel[ii, jj], ph_i=ph[ii], ph_j=ph[jj], rq_i=rq[ii], rq_j=rq[jj])


def analyse(per_circuit, draws):
    cols = {k: [] for k in ['lag1', 'lag0', 'inp', 'rdr', 'neg', 'strat', 'circ']}
    for c, r in enumerate(per_circuit):
        for k in ['lag1', 'lag0', 'inp', 'rdr', 'neg']:
            cols[k].append(r[k])
        cols['strat'].append((((c * 3 + r['ph_i']) * 3 + r['ph_j']) * 5 + r['rq_i']) * 5 + r['rq_j'])
        cols['circ'].append(np.full(len(r['lag1']), c))
    C = {k: np.concatenate(v) for k, v in cols.items()}
    nc = len(per_circuit)
    mult = np.random.default_rng(SEED).multinomial(nc, np.full(nc, 1 / nc), size=draws)
    out, bs = {}, {}
    for pname in ['inp', 'rdr']:
        for sname in ['lag1', 'lag0']:
            m = C[pname] | C['neg']
            r, b = auc_with_ci(C[sname][m], C[pname][m], C['strat'][m], C['circ'][m], mult)
            r['verdict'] = verdict(r); out[f'{pname}_{sname}'] = r; bs[f'{pname}_{sname}'] = b
    out['inp_lag1_minus_lag0'] = diff_ci(bs['inp_lag1'], bs['inp_lag0'], out['inp_lag1']['auc'], out['inp_lag0']['auc'])
    out['rdr_lag0_minus_lag1'] = diff_ci(bs['rdr_lag0'], bs['rdr_lag1'], out['rdr_lag0']['auc'], out['rdr_lag1']['auc'])
    out['n_circuits'] = nc
    out['n_input_storage_pairs'] = int(C['inp'].sum()); out['n_storage_reader_pairs'] = int(C['rdr'].sum())
    out['circuits_with_input_storage_pairs'] = int(sum(r['inp'].any() for r in per_circuit))
    return out


def predictions(o):
    return {'A10.1_inp_lag1_above_0.5': bool(o['inp_lag1']['lo'] > 0.5),
            'A10.2_inp_lag1_minus_lag0_above_0': bool(o['inp_lag1_minus_lag0']['lo'] > 0),
            'A10.3_rdr_lag0_minus_lag1_above_0': bool(o['rdr_lag0_minus_lag1']['lo'] > 0),
            'A10.4_rdr_lag1_not_above_0.5': bool(not (o['rdr_lag1']['lo'] > 0.5))}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--selftest', action='store_true'); a = ap.parse_args()
    if a.selftest:
        for T in LENGTHS:
            t0 = time.time(); score_circuit(9999, T)
            log(f'one circuit at {T} half-steps: {time.time() - t0:.1f} s')
        log('SELF-TEST PASSED')
        return
    hashes = {p.name: sha(p) for p in (PRED, Path(__file__).resolve(), HELPERS_FILE, ANALYSES / '02_latch_transfer_tests_A1_to_A8/nmos_circuit_latch_tests_A8.py')}
    for k, v in hashes.items():
        log(f'A10: {k} sha256 {v}')
    res = dict(hashes=hashes, seeds_new=NEW, seeds_A8=OLD, draws=DRAWS)
    for T in LENGTHS:
        new, old = [], []
        for seed in NEW:
            new.append(score_circuit(seed, T)); log(f'T={T} circuit {seed} scored')
        for seed in OLD:
            old.append(score_circuit(seed, T))
        o = analyse(new, DRAWS); o['predictions_met'] = predictions(o)
        res[f'T{T}'] = dict(new_50=o, pooled_60_descriptive=analyse(new + old, DRAWS))
        save('A10_results', res)
        log(f'T={T} done')
    log('done')


if __name__ == '__main__':
    main()
