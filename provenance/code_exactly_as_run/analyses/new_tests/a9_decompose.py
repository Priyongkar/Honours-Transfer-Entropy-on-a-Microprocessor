"""A9: exact decomposition of lag-1 TE for strictly locked same-phase pairs
(design: PREDICTIONS_A9.md).

Usage:
  python a9_decompose.py --selftest     synthetic data only: checks the identities and times the run
  python a9_decompose.py --check-load   loads the real inputs and prints array shapes only
  python a9_decompose.py                full run: recording, chip program 1, chip program 2
"""
import argparse, hashlib, time
from pathlib import Path
import numpy as np
from common import (HERE, RUN, ROOT, OUT, CHIP_DIR, tk, log, save, change, phase_class, boot_mult,
                    auc_with_ci, verdict, diff_ci, quintile, bound_all_sources)

PRED = HERE / 'PREDICTIONS_A9.md'
TOL = 1e-9
SEED, DRAWS = 2026, 1000


# ---------------------------------------------------------------- provenance
def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def start(tag):
    for p in (PRED, Path(__file__).resolve(), HERE / 'common.py'):
        log(f'{tag}: {p.name} sha256 {sha(p)}')
    return {p.name: sha(p) for p in (PRED, Path(__file__).resolve(), HERE / 'common.py')}


# ---------------------------------------------------------------- the identity
def h(p):
    p = np.clip(np.asarray(p, float), 0.0, 1.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        v = -p * np.log2(p) - (1 - p) * np.log2(1 - p)
    return np.where((p <= 0) | (p >= 1), 0.0, v)


def F(pi, c, a):
    """Lag-1 TE within the quiet cell (bits)."""
    m = 1 - (1 - pi) * np.asarray(a, float)
    x = np.divide(pi * c, m, out=np.zeros_like(m), where=m > 0)
    return h(pi * c) - m * h(x)


def G(pi, a):
    """Proposition-2 bound within the quiet cell (bits)."""
    m = 1 - (1 - pi) * np.asarray(a, float)
    x = np.divide(pi, m, out=np.zeros_like(m), where=m > 0)
    return h(pi) - m * h(x)


def target_terms(Xsp, y, k, p, par):
    """Quiet-cell terms for target y, taken as locked to parity p, with history k, for every
    source row of Xsp. Uses exactly the samples t = max(k, 1) .. T-1 of the lag-1 TE estimator."""
    y = np.asarray(y).astype(np.int64); T = len(y); k = int(k); st = max(k, 1)
    t = np.arange(st, T); n = T - st
    Q = tk._encode_history(y, k, st, T) == 0
    inph = par[t] == p
    nQ = int(Q.sum()); nQp = int((Q & inph).sum())
    q0 = nQ / n
    pi = nQp / nQ if nQ else np.nan
    c = float(y[t][Q & inph].mean()) if nQp else 0.0
    s = t - 1                                   # source sample
    Hm = par[s] == p
    for j in range(1, k):
        Hm &= y[s - j] == 0
    yH = y[s] == 1
    assert int((Hm & ~yH).sum()) == int((Q & ~inph).sum()), 'quiet-cell bookkeeping'
    V = np.zeros((T, 3))
    V[s[Hm], 0] = 1.0; V[s[Hm & yH], 1] = 1.0; V[s[Hm & ~yH], 2] = 1.0
    C = np.asarray(Xsp @ V)
    nH = int(Hm.sum()); nH0 = int((Hm & ~yH).sum())
    rows = C.shape[0]
    rt = C[:, 0] / nH if nH else np.zeros(rows)
    kap = np.divide(C[:, 1], C[:, 0], out=np.zeros(rows), where=C[:, 0] > 0)
    ct = float(yH[Hm].mean()) if nH else 0.0
    a = C[:, 2] / nH0 if nH0 else np.zeros(rows)
    with np.errstate(divide='ignore', invalid='ignore'):
        a_dec = np.where(ct < 1, rt * (1 - kap) / (1 - ct), 0.0)
    if nQ == 0:
        te = np.zeros(rows); B = np.zeros(rows)
    else:
        te = q0 * F(pi, c, a); B = q0 * G(pi, a)
    return dict(q0=q0, pi=pi, c=c, ct=ct, a=a, a_dec=a_dec, rt=rt, kap=kap, te=te, B=B)


def decompose(X, Xsp, K, ph, targets, par):
    """Formula terms for every (source, target) with target in `targets`."""
    Gn = X.shape[0]
    out = {key: np.full((Gn, Gn), np.nan) for key in ('te', 'B', 'a', 'a_dec', 'rt', 'kap')}
    tgt = {key: np.full(Gn, np.nan) for key in ('q0', 'pi', 'c', 'ct')}
    for j in targets:
        r = target_terms(Xsp, X[j], K[j], int(ph[j]), par)
        for key in out:
            out[key][:, j] = r[key]
        for key in tgt:
            tgt[key][j] = r[key]
    return out, tgt


def checks(TE, Bv, D, S):
    """C1-C3 over the boolean pair mask S (source, target)."""
    si, sj = np.nonzero(S)
    c1_te = float(np.nanmax(np.abs(D['te'][si, sj] - TE[si, sj]))) if len(si) else 0.0
    c1_b = float(np.nanmax(np.abs(D['B'][si, sj] - Bv[si, sj]))) if len(si) else 0.0
    c2 = float(np.nanmax(np.abs(D['a'][si, sj] - D['a_dec'][si, sj]))) if len(si) else 0.0
    c3_viol = 0
    for j in np.unique(sj):
        src = si[sj == j]
        o = np.argsort(D['a'][src, j], kind='stable')
        a_s = D['a'][src, j][o]; te_s = TE[src, j][o]
        # equal a may carry rounding-level differences; require TE non-decreasing in a
        c3_viol += int(np.sum((np.diff(te_s) < -1e-12) & (np.diff(a_s) > 0)))
    res = dict(C1_max_abs_te_minus_formula=c1_te, C1_max_abs_B_minus_formula=c1_b,
               C2_max_abs_a_minus_decomposition=c2, C3_monotonicity_violations=c3_viol)
    res['all_checks_pass'] = bool(c1_te < TOL and c1_b < TOL and c2 < TOL and c3_viol == 0)
    return res


# ---------------------------------------------------------------- tests
def deciles(v):
    return np.searchsorted(np.unique(np.quantile(v, np.linspace(.1, .9, 9))), v, side='right')


def run_tests(pos, jj, te, B, rt, kap, cT, q0T, piT, phase_stratum, mult):
    """A9.1-A9.4 and the descriptive AUCs over the pairs of S (arrays aligned by pair)."""
    out = {}

    def auc(name, score, strata):
        r, bs = auc_with_ci(np.asarray(score, float), pos, np.asarray(strata), jj, mult)
        r['verdict'] = verdict(r); out[name] = r
        return r, bs

    bd = deciles(B); cq = quintile(cT)
    r, _ = auc('A9.1_followership_within_target', kap, jj)
    r['prediction'] = 'above 0.5'; r['prediction_met'] = bool(r['lo'] > 0.5)
    rTE, bTE = auc('desc_TE_within_target', te, jj)
    rR, bR = auc('desc_activity_within_target', rt, jj)
    d = diff_ci(bTE, bR, rTE['auc'], rR['auc'])
    d['prediction'] = 'below 0'; d['prediction_met'] = bool(d['hi'] < 0)
    out['A9.2_TE_minus_activity_within_target'] = d
    r, _ = auc('A9.3_target_c_within_B_deciles', cT, bd)
    r['prediction'] = 'below 0.5'; r['prediction_met'] = bool(r['hi'] < 0.5)
    rB, bB = auc('desc_TE_within_B_deciles', te, bd)
    rBc, bBc = auc('A9.4_TE_within_B_decile_x_c_quintile', te, bd * 5 + cq)
    d = diff_ci(bBc, bB, rBc['auc'], rB['auc'])
    d['prediction'] = 'above 0'; d['prediction_met'] = bool(d['lo'] > 0)
    out['A9.4_primary_difference'] = d
    includes = bool(rBc['lo'] <= 0.5 <= rBc['hi'])
    if d['lo'] > 0 and includes:
        reading = 'accounted for'
    elif d['lo'] > 0 and rBc['hi'] < 0.5:
        reading = 'partly accounted for'
    elif d['lo'] <= 0 <= d['hi']:
        reading = 'not accounted for'
    else:
        reading = 'other: see intervals'
    out['A9.4_reading'] = reading
    out['A9.4_stratified_interval_includes_0.5'] = includes
    auc('desc_TE_within_phase', te, phase_stratum)
    auc('desc_not_followed_within_target', 1 - kap, jj)
    auc('desc_target_q0_within_B_deciles', q0T, bd)
    auc('desc_target_pi_within_B_deciles', piT, bd)
    out['n_pairs'] = int(len(pos)); out['n_wired'] = int(pos.sum())
    return out


# ---------------------------------------------------------------- recording
def load_recording():
    X = np.asarray(np.load(RUN / 'thesis_v5/idtxl/data_int8.npy')).astype(np.uint8)
    with np.load(RUN / 'te.npz') as z:
        TE1 = z['te_bits'].copy(); K = z['optimal_k'].astype(int)
    with np.load(OUT / 'A3_bound_arrays.npz') as z:
        B3, Wd, ph, strict = z['B'], z['wired'], z['phase'], z['strict']
    return X, TE1, K, B3, Wd, ph, strict


def run_recording(draws):
    X, TE1, K, B3, Wd, ph, strict = load_recording()
    Gn, T = X.shape; par = (np.arange(T) % 2).astype(np.int64)
    ph2, st2 = phase_class(X)
    assert (ph2 == ph).all() and (st2 == strict).all(), 'phase classes differ from A3'
    Xsp = tk.as_sparse_sources(X)
    off = ~np.eye(Gn, dtype=bool)
    locked = ph < 2
    same = (ph[:, None] == ph[None, :]) & locked[:, None] & off
    S = same & strict[:, None] & strict[None, :]
    approx = same & ~S
    t0 = time.time()
    D, tg = decompose(X, Xsp, K, ph, np.flatnonzero(locked), par)
    log(f'recording: formula terms for {locked.sum()} locked targets in {time.time() - t0:.0f} s')
    res = dict(checks=checks(TE1, B3, D, S))
    log(f"recording checks: {res['checks']}")
    ai, aj = np.nonzero(approx)
    resid = np.abs(D['te'][ai, aj] - TE1[ai, aj])
    res['approx_95pct_locked_pairs'] = dict(n=int(len(ai)), median_abs_resid=float(np.median(resid)),
                                            p95_abs_resid=float(np.quantile(resid, .95)), max_abs_resid=float(resid.max()))
    si, sj = np.nonzero(S)
    pos = Wd[si, sj]
    res['scope'] = dict(pairs_in_S=int(len(si)), wired_in_S=int(pos.sum()),
                        wired_same_phase_all=int(Wd[same].sum()), pairs_same_phase_all=int(same.sum()))
    mult = boot_mult(Gn, draws, SEED)
    res['tests'] = run_tests(pos, sj, TE1[si, sj], B3[si, sj], D['rt'][si, sj], D['kap'][si, sj],
                             tg['c'][sj], tg['q0'][sj], tg['pi'][sj], ph[sj], mult)
    return res


# ---------------------------------------------------------------- simulated chip
def load_chip(cdir):
    from sim6502 import Chip
    with np.load(Path(cdir) / 'states.npz') as z:
        states, nodes = z['states'], z['nodes']
    with np.load(Path(cdir) / 'interventions.npz') as z:
        U = z['nodes']
    chip = Chip(ROOT / 'data/visual6502')
    g_, c1_, c2_ = chip.transistor_matrix()
    A = tk.node_reference_from_transistors(g_, c1_, c2_, nodes, supply={chip.gnd, chip.pwr}).toarray() > 0
    return states, U, A[np.ix_(U, U)]


def run_chip(cdir, draws):
    states, U, W = load_chip(cdir)
    T = states.shape[0]; n = len(U); par = (np.arange(T) % 2).astype(np.int64)
    chg = change(states.T[U].astype(np.uint8))
    reps, ug = tk.group_identical(chg); Xg = chg[reps]; Gn = len(reps)
    ks = np.asarray([tk.ais_select_k(x, 6) for x in Xg])
    phg, strg = phase_class(Xg)
    Xsp = tk.as_sparse_sources(Xg)
    offg = ~np.eye(Gn, dtype=bool)
    Sg = (phg[:, None] == phg[None, :]) & (phg[:, None] < 2) & strg[:, None] & strg[None, :] & offg
    targets = np.flatnonzero(Sg.any(0))
    t0 = time.time()
    TEg = np.full((Gn, Gn), np.nan); Bg = np.full((Gn, Gn), np.nan)
    for j in targets:
        TEg[:, j] = tk.te_all_sources(Xsp, Xg[j], ks[j], delay=1)[0]
        Bg[:, j] = bound_all_sources(Xsp, Xg[j], ks[j], par, delay=1)
    D, tg = decompose(Xg, Xsp, ks, phg, targets, par)
    log(f'chip {Path(cdir).name}: {Gn} groups, {len(targets)} targets with strict same-phase sources, {time.time() - t0:.0f} s')
    res = dict(groups=int(Gn), checks=checks(TEg, Bg, D, Sg))
    log(f"chip checks: {res['checks']}")
    off = ~np.eye(n, dtype=bool)
    valid = (ug[:, None] >= 0) & (ug[None, :] >= 0) & (ug[:, None] != ug[None, :]) & off
    ii, jj = np.nonzero(valid)
    gi, gj = ug[ii], ug[jj]
    inS = Sg[gi, gj]
    pos_all = W[ii, jj]; neg_all = ~(W | W.T)[ii, jj]
    keep = inS & (pos_all | neg_all)
    ii, jj, gi, gj, pos = ii[keep], jj[keep], gi[keep], gj[keep], pos_all[keep]
    res['scope'] = dict(node_pairs_in_S=int(len(ii)), wired_in_S=int(pos.sum()),
                        wired_all=int(pos_all.sum()))
    mult = boot_mult(n, draws, SEED)
    res['tests'] = run_tests(pos, jj, TEg[gi, gj], Bg[gi, gj], D['rt'][gi, gj], D['kap'][gi, gj],
                             tg['c'][gj], tg['q0'][gj], tg['pi'][gj], phg[gj], mult)
    return res


# ---------------------------------------------------------------- self-test (synthetic only)
def synthetic(Gn, T, rng):
    par = np.arange(T) % 2
    ph = rng.integers(0, 2, Gn); X = np.zeros((Gn, T), np.uint8)
    kind = rng.choice(3, Gn, p=[.7, .2, .1])          # strict, 95%-locked, mixed
    for g in range(Gn):
        r = rng.uniform(.01, .2)
        if kind[g] == 2:
            X[g] = rng.random(T) < r / 2
        else:
            X[g] = (par == ph[g]) & (rng.random(T) < r)
            if kind[g] == 1:
                X[g] |= (par != ph[g]) & (rng.random(T) < r * .02)
    Wd = np.zeros((Gn, Gn), bool)
    for _ in range(Gn * 2):                            # followers: target copies source in phase
        s, t = rng.integers(0, Gn, 2)
        if s != t and ph[s] == ph[t] and kind[t] == 0 and kind[s] == 0:
            X[t] |= (X[s] == 1) & (rng.random(T) < rng.uniform(.3, .9))
            Wd[s, t] = True
    K = rng.integers(1, 13, Gn)
    return X, Wd, K


def selftest():
    rng = np.random.default_rng(7)
    log('self-test 1: identities on synthetic phase-locked data')
    X, Wd, K = synthetic(120, 30000, rng)
    ph, strict = phase_class(X); T = X.shape[1]; par = (np.arange(T) % 2).astype(np.int64)
    Xsp = tk.as_sparse_sources(X)
    TE = tk.te_matrix(X, K, delay=1)[0]
    Bv = np.full_like(TE, np.nan)
    for j in range(len(X)):
        Bv[:, j] = bound_all_sources(Xsp, X[j], K[j], par)
    off = ~np.eye(len(X), dtype=bool)
    S = (ph[:, None] == ph[None, :]) & (ph[:, None] < 2) & strict[:, None] & strict[None, :] & off
    D, tg = decompose(X, Xsp, K, ph, np.flatnonzero(ph < 2), par)
    ck = checks(TE, Bv, D, S)
    log(f"  checks pass: {ck['all_checks_pass']} (max |TE - formula| {ck['C1_max_abs_te_minus_formula']:.1e} bits)")
    si, sj = np.nonzero(S)
    out = run_tests(Wd[si, sj], sj, TE[si, sj], Bv[si, sj], D['rt'][si, sj], D['kap'][si, sj],
                    tg['c'][sj], tg['q0'][sj], tg['pi'][sj], ph[sj], boot_mult(len(X), 20, 1))
    log(f"  test machinery ran: {len(out)} outputs (synthetic labels; values not meaningful)")
    assert ck['all_checks_pass']

    log('self-test 2: timing at the recording\'s size (519 groups x 133,926 samples, synthetic)')
    X, Wd, K = synthetic(519, 133926, rng); ph, strict = phase_class(X)
    T = X.shape[1]; par = (np.arange(T) % 2).astype(np.int64); Xsp = tk.as_sparse_sources(X)
    t0 = time.time(); D, tg = decompose(X, Xsp, K, ph, np.flatnonzero(ph < 2), par); t_dec = time.time() - t0
    off = ~np.eye(len(X), dtype=bool)
    S = (ph[:, None] == ph[None, :]) & (ph[:, None] < 2) & strict[:, None] & strict[None, :] & off
    si, sj = np.nonzero(S)
    t0 = time.time()
    run_tests(Wd[si, sj], sj, D['te'][si, sj], D['B'][si, sj], D['rt'][si, sj], D['kap'][si, sj],
              tg['c'][sj], tg['q0'][sj], tg['pi'][sj], ph[sj], boot_mult(len(X), DRAWS, 1))
    t_test = time.time() - t0
    log(f'  recording-size decomposition {t_dec:.0f} s; tests with {DRAWS} draws on {len(si):,} pairs {t_test:.0f} s')

    log('self-test 3: timing of one chip program\'s TE, bound and AIS (610 groups x 60,000, synthetic)')
    X, Wd, K = synthetic(610, 60000, rng); ph, strict = phase_class(X)
    T = X.shape[1]; par = (np.arange(T) % 2).astype(np.int64); Xsp = tk.as_sparse_sources(X)
    t0 = time.time()
    ks = np.asarray([tk.ais_select_k(x, 6) for x in X])
    for j in range(len(X)):
        tk.te_all_sources(Xsp, X[j], ks[j], delay=1); bound_all_sources(Xsp, X[j], ks[j], par)
    decompose(X, Xsp, ks, ph, np.flatnonzero(ph < 2), par)   # the identity is defined for locked targets only
    log(f'  chip-size AIS + TE + bound + decomposition {time.time() - t0:.0f} s (tests as above)')
    log('SELF-TEST PASSED')


def check_load():
    X, TE1, K, B3, Wd, ph, strict = load_recording()
    log(f'recording: X {X.shape}, TE {TE1.shape}, K {K.shape}, B {B3.shape}, wired {Wd.shape}')
    for cdir in (CHIP_DIR, OUT / 'A7_chip'):
        states, U, W = load_chip(cdir)
        log(f'chip {cdir.name}: states {states.shape}, intervened nodes {U.shape}, wiring {W.shape}')
    log('LOAD CHECK PASSED')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--selftest', action='store_true'); ap.add_argument('--check-load', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if a.check_load:
        return check_load()
    hashes = start('A9')
    res = dict(hashes=hashes, draws=DRAWS, seed=SEED)
    res['recording'] = run_recording(DRAWS); save('A9_results', res)
    log('recording done')
    res['chip_program1'] = run_chip(CHIP_DIR, DRAWS); save('A9_results', res)
    log('chip program 1 done')
    res['chip_program2'] = run_chip(OUT / 'A7_chip', DRAWS); save('A9_results', res)
    log('chip program 2 done')
    log('done')


if __name__ == '__main__':
    main()
