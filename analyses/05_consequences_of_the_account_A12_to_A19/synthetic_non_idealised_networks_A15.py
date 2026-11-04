"""A15: non-idealised synthetic networks (hidden co-drivers, spurious events) scored as events and as
states (design: PREDICTIONS_A12.md). Network generator extends a11_synthetic.network.

Usage:
  python a15_synthetic.py --selftest   one network (seed outside the test): timing only
  python a15_synthetic.py              full run, 80 networks
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / 'src'))  # repository library
from te6502.paths import ANALYSES, DATA_ROOT, RESULTS, HELPERS_FILE, load_analysis  # noqa: E402,F401
import argparse, hashlib, time
from pathlib import Path
import numpy as np
from te6502.test_helpers import HERE, tk, log, save, auc_with_ci, verdict, diff_ci

PRED = ANALYSES / '05_consequences_of_the_account_A12_to_A19/design/predictions_A12_to_A19.md'
LAMBDAS = [0.5, 0.7, 0.85, 1.0]
RHO, NETS, N, T = 0.6, 20, 100, 50000
L0, L1 = 40, 30
HIDDEN, HID_FAN, HID_P, NOISE = 10, 4, 0.5, 0.01
SEED, DRAWS = 2026, 1000


def network(lam, seed):
    rng = np.random.default_rng(seed)
    phase = np.repeat([0, 1], N // 2); rng.shuffle(phase)
    layer = np.r_[np.zeros(L0, int), np.ones(L1, int), np.full(N - L0 - L1, 2)]
    rate = rng.uniform(0.02, 0.10, N)
    fast = np.zeros((N, N), bool); slow = np.zeros((N, N), bool)
    for j in range(L0, N):
        pool = np.flatnonzero((layer < layer[j]) & (phase == phase[j]))
        for i in rng.choice(pool, size=min(2, len(pool)), replace=False):
            fast[i, j] = True
        if rng.random() < 0.5:
            slow[rng.choice(np.flatnonzero((layer == 0) & (phase != phase[j]))), j] = True
    par = np.arange(T) % 2
    p_on = np.where(par[None, :] == phase[:, None], 2 * lam, 2 * (1 - lam)) * rate[:, None]
    E = rng.random((N, T)) < p_on
    # hidden co-drivers: each drives HID_FAN observed units of its phase within the sample
    hph = rng.integers(0, 2, HIDDEN)
    H = rng.random((HIDDEN, T)) < np.where(par[None, :] == hph[:, None], 2 * lam, 2 * (1 - lam)) * 0.06
    targets_of = [rng.choice(np.flatnonzero((phase == hph[h]) & (layer > 0)), HID_FAN, replace=False) for h in range(HIDDEN)]
    for h in range(HIDDEN):
        for j in targets_of[h]:
            E[j] |= H[h] & (rng.random(T) < HID_P)
    for j in range(L0, N):
        for i in np.flatnonzero(fast[:, j]):
            E[j] |= E[i] & (rng.random(T) < RHO)
        for i in np.flatnonzero(slow[:, j]):
            E[j, 1:] |= E[i, :-1] & (rng.random(T - 1) < RHO)
    E |= rng.random((N, T)) < NOISE                    # spurious events in either phase
    return E.astype(np.uint8), phase, fast, slow


def score(E):
    par = (np.arange(T) % 2).astype(np.uint8)
    ks = np.asarray([tk.ais_select_k(x, 6) for x in E])
    ev1 = tk.te_matrix(E, ks, delay=1)[0]
    ev0p = tk.te_matrix(E, ks, delay=0, cond=lambda t: par[None, :])[0]
    S = (np.cumsum(E, axis=1) % 2).astype(np.uint8)    # states whose changes are the events
    kss = np.asarray([tk.ais_select_k(x, 6) for x in S])
    st1 = tk.te_matrix(S, kss, delay=1)[0]
    return ev1, ev0p, st1


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--selftest', action='store_true'); a = ap.parse_args()
    if a.selftest:
        t0 = time.time(); E, *_ = network(1.0, 999_999); t1 = time.time(); score(E)
        log(f'one network: simulate {t1 - t0:.1f} s, score {time.time() - t1:.1f} s; full run about {80 * (time.time() - t0) / 60:.0f} min')
        log('SELF-TEST PASSED'); return
    hashes = {p.name: sha(p) for p in (PRED, Path(__file__).resolve(), HELPERS_FILE)}
    for k, v in hashes.items():
        log(f'A15: {k} sha256 {v}')
    res = dict(hashes=hashes, conditions={}, draws=DRAWS); boots = {}
    for ci, lam in enumerate(LAMBDAS):
        cols = {k: [] for k in ['ev1', 'ev0p', 'st1', 'pos', 'neg', 'unit']}
        for rep in range(NETS):
            E, phase, fast, slow = network(lam, 20_000 + 100 * ci + rep)
            ev1, ev0p, st1 = score(E)
            off = ~np.eye(N, dtype=bool); wired = fast | slow | fast.T | slow.T
            same = phase[:, None] == phase[None, :]
            ii, jj = np.nonzero(off)
            cols['ev1'].append(ev1[ii, jj]); cols['ev0p'].append(ev0p[ii, jj]); cols['st1'].append(st1[ii, jj])
            cols['pos'].append(fast[ii, jj]); cols['neg'].append((same & ~wired)[ii, jj]); cols['unit'].append(rep * N + jj)
        C = {k: np.concatenate(v) for k, v in cols.items()}
        mult = np.random.default_rng(SEED + ci).multinomial(NETS * N, np.full(NETS * N, 1 / (NETS * N)), size=DRAWS)
        m = C['pos'] | C['neg']; out = {}
        for name in ['ev1', 'ev0p', 'st1']:
            r, b = auc_with_ci(C[name][m], C['pos'][m], C['unit'][m], C['unit'][m], mult)
            r['verdict'] = verdict(r); out[name] = r; boots[(lam, name)] = b
        res['conditions'][f'lambda={lam}'] = out; save('A15_results', res)
        log(f'A15 condition lambda={lam} done')
    o1 = res['conditions']['lambda=1.0']; o0 = res['conditions']['lambda=0.5']
    d_lock = diff_ci(boots[(1.0, 'ev1')], boots[(0.5, 'ev1')], o1['ev1']['auc'], o0['ev1']['auc'])
    d_code = diff_ci(boots[(1.0, 'st1')], boots[(1.0, 'ev1')], o1['st1']['auc'], o1['ev1']['auc'])
    res['A15.2_lambda1_minus_lambda0.5_events'] = d_lock; res['A15.3_states_minus_events_lambda1'] = d_code
    res['predictions_met'] = {'A15.1': bool(o1['ev1']['hi'] < 0.5), 'A15.2': bool(d_lock['hi'] < 0),
                              'A15.3': bool(d_code['lo'] > 0),
                              'A15.4': bool(all(res['conditions'][f'lambda={l}']['ev0p']['lo'] > 0.5 for l in LAMBDAS))}
    save('A15_results', res); log(f"done: {res['predictions_met']}")


if __name__ == '__main__':
    main()
