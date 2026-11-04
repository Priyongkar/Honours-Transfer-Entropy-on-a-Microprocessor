"""A12, A13 (recording), A14 (recording), A18, A19 on the reduced recording
(design: PREDICTIONS_A12.md).

Usage:
  python a12_recording.py --selftest   synthetic data only: estimator checks and timing
  python a12_recording.py              full run
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / 'src'))  # repository library
from te6502.paths import ANALYSES, DATA_ROOT, RESULTS, HELPERS_FILE, load_analysis  # noqa: E402,F401
import argparse, hashlib, time
from pathlib import Path
import numpy as np
from scipy import sparse
from te6502.test_helpers import HERE, RUN, OUT, tk, log, save, phase_class, boot_mult, auc_with_ci, verdict, diff_ci, quintile
a9 = load_analysis('04_followership_mechanism_A9_to_A11/followership_decomposition_A9.py')

PRED = ANALYSES / '05_consequences_of_the_account_A12_to_A19/design/predictions_A12_to_A19.md'
SEED, DRAWS = 2026, 1000


def start(tag):
    h = {}
    for p in (PRED, Path(__file__).resolve(), HELPERS_FILE, ANALYSES / '04_followership_mechanism_A9_to_A11/followership_decomposition_A9.py'):
        h[p.name] = hashlib.sha256(Path(p).read_bytes()).hexdigest()
        log(f'{tag}: {p.name} sha256 {h[p.name]}')
    return h


def lock_index(X):
    n = X.sum(1).astype(float)
    even = X[:, 0::2].sum(1) / np.maximum(n, 1)
    return np.maximum(even, 1 - even)


def history_codes(y, k, st, T):
    code = tk._encode_history(np.asarray(y).astype(np.int64), int(k), st, T)
    _, code = np.unique(code, return_inverse=True)
    return code.reshape(-1)


def cmi_from_counts(Cs, Ntot, n):
    """Plug-in I(source; y | h) in bits from per-source-value count arrays Cs[v] of shape (G, S, 2)."""
    ns = Ntot.sum(axis=1)
    out = np.zeros(Cs[0].shape[0])
    for C in Cs:
        n_sx = C.sum(axis=2, keepdims=True)
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = (C * ns[None, :, None]) / (n_sx * Ntot[None])
            out += np.where(C > 0, C * np.log2(ratio), 0.0).sum(axis=(1, 2))
    return np.maximum(out / n, 0.0)


def te_block2_and_lag1(Xsp, P2, y, k):
    """I((X_{t-1}, X_{t-2}); Y_t | Y^(k)) and I(X_{t-1}; Y_t | Y^(k)) on the same samples t = max(k,2)..T-1.
    Xsp: (G, T) csr; P2: (G, T-1) csr with P2[:, i] = X[:, i] * X[:, i+1]."""
    y = np.asarray(y).astype(np.int64); T = len(y); st = max(int(k), 2); t = np.arange(st, T); n = T - st
    h = history_codes(y, k, st, T); S = int(h.max()) + 1
    j = h * 2 + y[st:T]
    ones = np.ones(n)
    Ma = sparse.csr_matrix((ones, (t - 1, j)), shape=(T, 2 * S))
    Mb = sparse.csr_matrix((ones, (t - 2, j)), shape=(T, 2 * S))
    Mc = sparse.csr_matrix((ones, (t - 2, j)), shape=(T - 1, 2 * S))
    G = Xsp.shape[0]
    C1x = np.asarray((Xsp @ Ma).todense()).reshape(G, S, 2)
    Cx1 = np.asarray((Xsp @ Mb).todense()).reshape(G, S, 2)
    C11 = np.asarray((P2 @ Mc).todense()).reshape(G, S, 2)
    Ntot = np.bincount(j, minlength=2 * S).astype(np.float64).reshape(S, 2)
    n3 = C11; n2 = C1x - C11; n1 = Cx1 - C11; n0 = Ntot[None] - C1x - Cx1 + C11
    block2 = cmi_from_counts([n0, n1, n2, n3], Ntot, n)
    lag1 = cmi_from_counts([Ntot[None] - C1x, C1x], Ntot, n)
    return block2, lag1


def phase_split_te(Xsp, y, k, p):
    """Lag-1 TE for every source, split into the parts summed over samples of parity p and of the other parity.
    Returns (te_total_bits, sum_in, sum_off) with sum_in + sum_off = te_total * n."""
    y = np.asarray(y).astype(np.int64); T = len(y); st = max(int(k), 1); t = np.arange(st, T); n = T - st
    h = history_codes(y, k, st, T); S = int(h.max()) + 1
    j = h * 2 + y[st:T]; inph = (t % 2) == p
    G = Xsp.shape[0]
    def counts(mask):
        M = sparse.csr_matrix((np.ones(mask.sum()), (t[mask] - 1, j[mask])), shape=(T, 2 * S))
        C1 = np.asarray((Xsp @ M).todense()).reshape(G, S, 2)
        Nt = np.bincount(j[mask], minlength=2 * S).astype(np.float64).reshape(S, 2)
        return C1, Nt
    C1_in, N_in = counts(inph); C1_off, N_off = counts(~inph)
    C1 = C1_in + C1_off; Ntot = N_in + N_off; C0 = Ntot[None] - C1
    ns = Ntot.sum(axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        l1 = np.log2((C1 / C1.sum(2, keepdims=True)) / (Ntot / ns[:, None])[None])
        l0 = np.log2((C0 / C0.sum(2, keepdims=True)) / (Ntot / ns[:, None])[None])
    l1 = np.where(C1 > 0, l1, 0.0); l0 = np.where(C0 > 0, l0, 0.0)
    C0_in = N_in[None] - C1_in; C0_off = N_off[None] - C1_off
    s_in = (C1_in * l1 + C0_in * l0).sum(axis=(1, 2))
    s_off = (C1_off * l1 + C0_off * l0).sum(axis=(1, 2))
    return (s_in + s_off) / n, s_in, s_off


def auc(out, name, score, pos, strata, target, mult):
    r, b = auc_with_ci(np.asarray(score, float), np.asarray(pos, bool), np.asarray(strata), np.asarray(target), mult)
    r['verdict'] = verdict(r); out[name] = r
    return r, b


# ---------------------------------------------------------------- full run
def main_run():
    hashes = start('A12-recording')
    X = np.asarray(np.load(RUN / 'thesis_v5/idtxl/data_int8.npy')).astype(np.uint8)
    with np.load(RUN / 'te.npz') as z:
        TE1 = z['te_bits'].copy(); K = z['optimal_k'].astype(int)
    TE0 = np.load(RUN / 'thesis_v5/te_lag0.npz')['te_bits']
    with np.load(OUT / 'A3_bound_arrays.npz') as z:
        B3, strict = z['B'], z['strict']
    ref = np.load(OUT / 'recording_reference.npz'); Wd, Lat, ph, rate = ref['wired'], ref['latch'], ref['phase'], ref['rate']
    G, T = X.shape; par = (np.arange(T) % 2).astype(np.int64)
    XS = tk.as_sparse_sources(X); off = ~np.eye(G, dtype=bool)
    clock_g = int(np.argmax(rate)); mult = boot_mult(G, DRAWS, SEED)
    res = dict(hashes=hashes, draws=DRAWS)

    # ---- A12
    L = lock_index(X); cls = np.where(L < 0.95, 0, np.where(L < 1, 1, 2))
    keep = off.copy(); keep[clock_g, :] = False; keep[:, clock_g] = False
    ii, jj = np.nonzero(keep & (Wd | ~(Wd | Wd.T))); pos = Wd[ii, jj]
    o = {}; bs = {}
    for sname, Smat in [('lag1', TE1), ('lag0', TE0)]:
        s = Smat[ii, jj]
        for c, cn in enumerate(['mixed', 'locked', 'strict']):
            m = cls[jj] == c
            r, b = auc(o, f'{sname}_{cn}', s[m], pos[m], jj[m], jj[m], mult)
            r['targets_with_wired'] = int(len(np.unique(jj[m & pos]))); bs[f'{sname}_{cn}'] = b
    d = diff_ci(bs['lag1_mixed'], bs['lag1_strict'], o['lag1_mixed']['auc'], o['lag1_strict']['auc'])
    o['A12.1_mixed_minus_strict'] = d
    o['predictions_met'] = {'A12.1': bool(d['lo'] > 0), 'A12.2': bool(o['lag1_strict']['hi'] < 0.5),
                            'A12.3': bool(o['lag1_mixed']['lo'] > 0.5)}
    o['targets_per_class'] = {cn: int((cls == c).sum()) for c, cn in enumerate(['mixed', 'locked', 'strict'])}
    res['A12'] = o; save('A12_results', res); log(f"A12 done: {o['predictions_met']}")

    # ---- A13 (recording)
    rq = quintile(rate)
    N_ = ~(Wd | Wd.T) & ~(Lat | Lat.T) & off
    ii, jj = np.nonzero((Lat | N_) & off)
    base = ((ph[ii] * 3 + ph[jj]) * 5 + rq[ii]) * 5 + rq[jj]
    S_all, P_all, ST_all, TG_all, COV_all = [], [], [], [], []
    for u in range(2, 8):
        Su = np.full((G, G), np.nan)
        for t_ in range(G):
            Su[:, t_] = tk.te_all_sources(XS, X[t_], K[t_], delay=u)[0]
        S_all.append(Su[ii, jj]); P_all.append(Lat[ii, jj]); ST_all.append(base * 10 + u); TG_all.append(jj)
        COV_all.append(K[jj] >= u - 1)
        log(f'A13 recording lag {u} computed')
    s, p_, stt, tg, cov = (np.concatenate(v) for v in (S_all, P_all, ST_all, TG_all, COV_all))
    o = {}
    rc, bc = auc(o, 'covered', s[cov], p_[cov], stt[cov], tg[cov], mult)
    o['uncovered_latch_pairs'] = int(p_[~cov].sum()); o['covered_latch_pairs'] = int(p_[cov].sum())
    if p_[~cov].sum() > 0:
        ru, bu = auc(o, 'uncovered', s[~cov], p_[~cov], stt[~cov], tg[~cov], mult)
        o['A13.2_uncovered_minus_covered'] = diff_ci(bu, bc, ru['auc'], rc['auc'])
    o['descriptive_only'] = bool(p_[~cov].sum() < 10)
    res['A13_recording'] = o; save('A12_results', res); log('A13 recording done')

    # ---- A14 (recording)
    D, tg9 = a9.decompose(X, XS, K, ph, np.flatnonzero(ph < 2), par)
    same = (ph[:, None] == ph[None, :]) & (ph[:, None] < 2) & off
    S = same & strict[:, None] & strict[None, :]
    si, sj = np.nonzero(S); pos = Wd[si, sj]
    bd = a9.deciles(B3[si, sj]); q0q = quintile(tg9['q0'][sj]); piq = quintile(tg9['pi'][sj])
    o = {}
    rB, bB = auc(o, 'TE_within_B_deciles', TE1[si, sj], pos, bd, sj, mult)
    rX, bX = auc(o, 'A14.2_TE_within_B_q0_pi', TE1[si, sj], pos, (bd * 5 + q0q) * 5 + piq, sj, mult)
    d = diff_ci(bX, bB, rX['auc'], rB['auc']); o['A14.1_difference'] = d
    inc = bool(rX['lo'] <= 0.5 <= rX['hi'])
    o['reading'] = 'accounted for' if (d['lo'] > 0 and inc) else ('partly accounted for' if d['lo'] > 0 else 'not accounted for')
    o['predictions_met'] = {'A14.1': bool(d['lo'] > 0), 'A14.2': inc}
    res['A14_recording'] = o; save('A12_results', res); log(f"A14 recording: {o['reading']}")

    # ---- A18
    P2 = XS[:, :-1].multiply(XS[:, 1:]).tocsr()
    BL2 = np.full((G, G), np.nan); L1s = np.full((G, G), np.nan)
    for t_ in range(G):
        b2, l1 = te_block2_and_lag1(XS, P2, X[t_], K[t_]); BL2[:, t_] = b2; L1s[:, t_] = l1
    np.fill_diagonal(BL2, np.nan); np.fill_diagonal(L1s, np.nan)
    o = {}
    ii, jj = np.nonzero(off); pos = Wd[ii, jj]; one = np.zeros(len(ii), int)
    r, _ = auc(o, 'A18.1_block2_group_auc', BL2[ii, jj], pos, one, jj, mult)
    auc(o, 'lag1_group_auc_same_samples', L1s[ii, jj], pos, one, jj, mult)
    r2, _ = auc(o, 'A18.2_block2_within_target_strict_same_phase', BL2[si, sj], Wd[si, sj], sj, sj, mult)
    auc(o, 'lag1_within_target_strict_same_phase', L1s[si, sj], Wd[si, sj], sj, sj, mult)
    o['predictions_met'] = {'A18.1': bool(r['lo'] <= 0.5), 'A18.2': bool(r2['lo'] <= 0.5)}
    res['A18'] = o; save('A12_results', res); log(f"A18 done: {o['predictions_met']}")

    # ---- A19
    SH = np.full((G, G), np.nan); TT = np.full((G, G), np.nan)
    for t_ in np.flatnonzero(ph < 2):
        te, s_in, s_off = phase_split_te(XS, X[t_], K[t_], int(ph[t_]))
        tot = s_in + s_off
        TT[:, t_] = te
        with np.errstate(divide='ignore', invalid='ignore'):
            SH[:, t_] = np.where(te > 1e-6, s_in / tot, np.nan)
    np.fill_diagonal(SH, np.nan)
    o = {}
    stor = np.unique(np.nonzero(Lat)[1]); stor = stor[ph[stor] < 2]
    m = np.zeros((G, G), bool); m[:, stor] = True; m &= off & np.isfinite(SH)
    ii, jj = np.nonzero(m); pos = Lat[ii, jj]
    r, _ = auc(o, 'A19.1_in_phase_share_latch_vs_others', SH[ii, jj], pos, jj, jj, mult)
    lat_sh = SH[Lat & np.isfinite(SH)]; str_sh = SH[S & np.isfinite(SH)]
    o['median_share_latch'] = float(np.median(lat_sh)) if len(lat_sh) else None
    o['median_share_strict_same_phase'] = float(np.median(str_sh)) if len(str_sh) else None
    o['n_latch_scored'] = int(len(lat_sh)); o['n_strict_scored'] = int(len(str_sh))
    o['max_abs_check_te_vs_saved'] = float(np.nanmax(np.abs(TT[Lat] - TE1[Lat]))) if Lat.any() else None
    o['predictions_met'] = {'A19.1': bool(r['lo'] > 0.5),
                            'A19.2': bool((o['median_share_latch'] or 0) > 0.5 and (o['median_share_strict_same_phase'] or 1) < 0.5)}
    res['A19'] = o; save('A12_results', res); log(f"A19 done: {o['predictions_met']}")
    log('done')


# ---------------------------------------------------------------- self-test (synthetic only)
def brute_cmi(src_sym, y, k):
    from collections import Counter
    T = len(y); st = max(k, 2); n = T - st
    c = Counter()
    for t in range(st, T):
        h = tuple(y[t - k:t]); c[(h, src_sym[t], y[t])] += 1
    hx = Counter(); hy = Counter(); hh = Counter()
    for (h, x, yy), v in c.items():
        hx[(h, x)] += v; hy[(h, yy)] += v; hh[h] += v
    return sum(v / n * np.log2(v * hh[h] / (hx[(h, x)] * hy[(h, yy)])) for (h, x, yy), v in c.items())


def selftest():
    rng = np.random.default_rng(5)
    X = (rng.random((6, 3000)) < 0.2).astype(np.uint8)
    X[1, 1:] |= X[0, :-1]; X[2, 2:] |= X[0, :-2]
    XS = tk.as_sparse_sources(X); P2 = XS[:, :-1].multiply(XS[:, 1:]).tocsr()
    worst = 0
    for k in (1, 2, 3):
        b2, l1 = te_block2_and_lag1(XS, P2, X[2], k)
        for g in range(6):
            sym = np.zeros(3000, int); sym[2:] = 2 * X[g, 1:-1] + X[g, :-2]
            worst = max(worst, abs(b2[g] - brute_cmi(sym, X[2], k)))
            sym1 = np.zeros(3000, int); sym1[1:] = X[g, :-1]
            worst = max(worst, abs(l1[g] - brute_cmi(sym1, X[2], k)))
    log(f'block-2 and lag-1 estimator vs brute force: max abs diff {worst:.1e} bits')
    assert worst < 1e-10
    te, s_in, s_off = phase_split_te(XS, X[1], 2, 0)
    ref = tk.te_all_sources(XS, X[1], 2, delay=1)[0]
    d = float(np.max(np.abs(te - ref)))
    log(f'phase-split TE sums to the toolkit TE: max abs diff {d:.1e} bits'); assert d < 1e-10
    Xb = (rng.random((519, 133926)) < 0.07).astype(np.uint8); XS = tk.as_sparse_sources(Xb)
    P2 = XS[:, :-1].multiply(XS[:, 1:]).tocsr()
    t0 = time.time()
    for t_ in range(10):
        te_block2_and_lag1(XS, P2, Xb[t_], 12); phase_split_te(XS, Xb[t_], 12, 0)
    per = (time.time() - t0) / 10
    log(f'timing at recording size: {per:.2f} s per target for A18+A19; about {per * 519 / 60:.0f} min for all targets')
    log('SELF-TEST PASSED')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--selftest', action='store_true'); a = ap.parse_args()
    selftest() if a.selftest else main_run()
