"""A11: generality of the followership penalty on synthetic phase-locked networks
(design: PREDICTIONS_A11.md).

Usage:
  python a11_synthetic.py --selftest   one network (seed outside the test): identity check and timing only
  python a11_synthetic.py              full run, 300 networks
"""
import argparse, hashlib, time
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from common import HERE, tk, log, save, phase_class, auc_with_ci, verdict, diff_ci, bound_all_sources
from a9_decompose import decompose, checks

PRED = HERE / 'PREDICTIONS_A11.md'
LAMBDAS = [0.5, 0.7, 0.85, 0.95, 1.0]
RHOS = [0.3, 0.6, 0.9]
NETS, N, T = 20, 100, 50000
L0, L1 = 40, 30                      # L2 = N - L0 - L1
SEED, DRAWS = 2026, 1000


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def network(lam, rho, seed):
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
    for j in range(L0, N):            # layer order: drivers are final before they propagate
        for i in np.flatnonzero(fast[:, j]):
            E[j] |= E[i] & (rng.random(T) < rho)
        for i in np.flatnonzero(slow[:, j]):
            E[j, 1:] |= E[i, :-1] & (rng.random(T - 1) < rho)
    return E.astype(np.uint8), phase, fast, slow


def score(E):
    par = (np.arange(T) % 2).astype(np.uint8)
    ks = np.asarray([tk.ais_select_k(x, 6) for x in E])
    L1s = tk.te_matrix(E, ks, delay=1)[0]
    L0s = tk.te_matrix(E, ks, delay=0)[0]
    L0p = tk.te_matrix(E, ks, delay=0, cond=lambda t: par[None, :])[0]
    return ks, L1s, L0s, L0p


def identity_check(E, ks, L1s):
    """Remark 3 on strictly locked same-phase pairs of one network."""
    ph, strict = phase_class(E)
    T_ = E.shape[1]; par = (np.arange(T_) % 2).astype(np.int64)
    Xsp = tk.as_sparse_sources(E); off = ~np.eye(len(E), dtype=bool)
    S = (ph[:, None] == ph[None, :]) & (ph[:, None] < 2) & strict[:, None] & strict[None, :] & off
    tg = np.flatnonzero(S.any(0))
    Bv = np.full_like(L1s, np.nan)
    for j in tg:
        Bv[:, j] = bound_all_sources(Xsp, E[j], ks[j], par)
    D, _ = decompose(E, Xsp, ks, ph, tg, par)
    return checks(L1s, Bv, D, S), D, S


def locking_index(E):
    n = E.sum(1).astype(float)
    even = E[:, 0::2].sum(1) / np.maximum(n, 1)
    return float(np.mean(np.maximum(even, 1 - even)[n > 0]))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--selftest', action='store_true'); a = ap.parse_args()
    if a.selftest:
        t0 = time.time(); E, phase, fast, slow = network(1.0, 0.9, 999_999); t_sim = time.time() - t0
        t0 = time.time(); ks, L1s, L0s, L0p = score(E); t_sc = time.time() - t0
        ck, _, _ = identity_check(E, ks, L1s)
        log(f'one network: simulate {t_sim:.1f} s, score {t_sc:.1f} s; identity checks pass: {ck["all_checks_pass"]}')
        log(f'estimated full run: {300 * (t_sim + t_sc) / 60:.0f} min plus bootstrap')
        assert ck['all_checks_pass']
        log('SELF-TEST PASSED')
        return
    hashes = {p.name: sha(p) for p in (PRED, Path(__file__).resolve(), HERE / 'common.py', HERE / 'a9_decompose.py')}
    for k, v in hashes.items():
        log(f'A11: {k} sha256 {v}')
    res = dict(hashes=hashes, conditions={}, draws=DRAWS)
    boot = {}
    for ci, (lam, rho) in enumerate([(l, r) for l in LAMBDAS for r in RHOS]):
        cols = {k: [] for k in ['l1', 'l0', 'l0p', 'kap', 'fpos', 'fneg', 'spos', 'sneg', 'unit']}
        lock, ident = [], []
        for rep in range(NETS):
            E, phase, fast, slow = network(lam, rho, 10_000 + 100 * ci + rep)
            ks, L1s, L0s, L0p = score(E)
            lock.append(locking_index(E))
            kap = np.full((N, N), np.nan)
            if lam == 1.0:
                ck, D, _ = identity_check(E, ks, L1s); ident.append(ck); kap = D['kap']
            off = ~np.eye(N, dtype=bool)
            wired = fast | slow | fast.T | slow.T
            same = phase[:, None] == phase[None, :]
            ii, jj = np.nonzero(off)
            cols['l1'].append(L1s[ii, jj]); cols['l0'].append(L0s[ii, jj]); cols['l0p'].append(L0p[ii, jj])
            cols['kap'].append(kap[ii, jj])
            cols['fpos'].append(fast[ii, jj]); cols['fneg'].append((same & ~wired)[ii, jj])
            cols['spos'].append(slow[ii, jj]); cols['sneg'].append((~same & ~wired)[ii, jj])
            cols['unit'].append(rep * N + jj)
        C = {k: np.concatenate(v) for k, v in cols.items()}
        mult = np.random.default_rng(SEED + ci).multinomial(NETS * N, np.full(NETS * N, 1 / (NETS * N)), size=DRAWS)
        out = {}

        def auc(name, s, pos, neg):
            m = pos | neg
            r, bs = auc_with_ci(s[m], pos[m], C['unit'][m], C['unit'][m], mult)
            r['verdict'] = verdict(r); out[name] = r
            return bs
        boot[(lam, rho)] = auc('fast_lag1', C['l1'], C['fpos'], C['fneg'])
        auc('fast_lag0', C['l0'], C['fpos'], C['fneg'])
        auc('fast_lag0_given_parity', C['l0p'], C['fpos'], C['fneg'])
        auc('delayed_lag1', C['l1'], C['spos'], C['sneg'])
        if lam == 1.0:
            m = ~np.isnan(C['kap'])
            r, _ = auc_with_ci(C['kap'][m], C['fpos'][m], C['unit'][m], C['unit'][m], mult) if m.any() else ({}, None)
            out['fast_followership_kappa'] = r
            out['identity_checks_pass_all_networks'] = bool(all(c['all_checks_pass'] for c in ident))
            out['identity_max_abs_te_minus_formula'] = float(max(c['C1_max_abs_te_minus_formula'] for c in ident))
        out['locking_index_mean'] = float(np.mean(lock))
        out['n_fast_pos'] = int(C['fpos'].sum()); out['n_delayed_pos'] = int(C['spos'].sum())
        res['conditions'][f'lambda={lam},rho={rho}'] = out
        save('A11_results', res)
        log(f'condition {ci + 1}/15 (lambda {lam}, rho {rho}) done')
    pred = {}
    for rho in RHOS:
        o1 = res['conditions'][f'lambda=1.0,rho={rho}']; o0 = res['conditions'][f'lambda=0.5,rho={rho}']
        d = diff_ci(boot[(1.0, rho)], boot[(0.5, rho)], o1['fast_lag1']['auc'], o0['fast_lag1']['auc'])
        res[f'diff_lambda1_minus_lambda0.5_rho={rho}'] = d
        if rho in (0.6, 0.9):
            pred[f'A11.1_rho={rho}'] = bool(o1['fast_lag1']['hi'] < 0.5)
            pred[f'A11.2_diff_rho={rho}'] = bool(d['hi'] < 0)
        pred[f'A11.2_blind_at_0.5_rho={rho}'] = bool(o0['fast_lag1']['lo'] <= 0.5 <= o0['fast_lag1']['hi'])
    conds = res['conditions'].values()
    pred['A11.3_all_15'] = bool(all(c['fast_lag0_given_parity']['lo'] > 0.5 for c in conds))
    pred['A11.4_all_15'] = bool(all(c['delayed_lag1']['lo'] > 0.5 for c in conds))
    res['predictions_met'] = pred
    li = [c['locking_index_mean'] for c in conds]; au = [c['fast_lag1']['auc'] for c in conds]
    rs = spearmanr(li, au)
    res['diagnostic_spearman_locking_vs_fast_lag1_auc'] = dict(rho=float(rs.statistic), p=float(rs.pvalue))
    save('A11_results', res)
    log('done')


if __name__ == '__main__':
    main()
