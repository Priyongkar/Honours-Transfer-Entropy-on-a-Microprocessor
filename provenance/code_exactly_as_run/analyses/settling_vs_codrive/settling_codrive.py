"""Settling versus co-drive on the simulated 6502 (design: PREDICTIONS.md).

Compares lag-0 CMI of directly influencing wired node pairs with that of
co-driven, non-influencing, unwired pairs matched on phase class, change
rates and co-drive strength. Reads the saved chip outputs; simulates nothing.

Usage: python settling_codrive.py [--draws 1000] [--quick]
"""
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'notebooks'))
import te6502_toolkit as tk          # noqa: E402
from sim6502 import Chip             # noqa: E402

CHIP_DIR = ROOT / 'results/v5_bitcorrect_every_third_mi_ais_te_k1_12_L1_d1_ais_x1/thesis_v5/simulated_chip'
OUT = HERE / 'output'
K_MAX, SEED, LOCK = 6, 2026, 0.95


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


def quintile(v):
    """Quintile index 0-4 of v (ties share a bin)."""
    edges = np.unique(np.quantile(v, [.2, .4, .6, .8]))
    return np.searchsorted(edges, v, side='right')


class StratifiedAUC:
    """AUC of score, positives vs negatives, compared only within strata,
    pooled over strata by number of comparisons; reweightable by a
    multiplicity per target unit for the target-cluster bootstrap."""

    def __init__(self, score, is_pos, stratum, target):
        keep = np.isin(stratum, np.intersect1d(stratum[is_pos], stratum[~is_pos]))
        s, p, st, tg = score[keep], is_pos[keep], stratum[keep], target[keep]
        neg_order = np.lexsort((s[~p], st[~p]))
        self.neg_s, self.neg_st, self.neg_t = s[~p][neg_order], st[~p][neg_order], tg[~p][neg_order]
        self.pos_s, self.pos_st, self.pos_t = s[p], st[p], tg[p]
        # negatives of each stratum occupy [start, stop) in the sorted array
        uniq, start = np.unique(self.neg_st, return_index=True)
        stop = np.r_[start[1:], len(self.neg_st)]
        k = np.searchsorted(uniq, self.pos_st)
        self.lo = np.empty(len(self.pos_s), np.int64); self.hi = np.empty_like(self.lo)
        for j, (a, b) in enumerate(zip(start, stop)):
            m = k == j
            self.lo[m] = a + np.searchsorted(self.neg_s[a:b], self.pos_s[m], side='left')
            self.hi[m] = a + np.searchsorted(self.neg_s[a:b], self.pos_s[m], side='right')
        self.start_of_pos = start[k]; self.stop_of_pos = stop[k]
        self.n_pos, self.n_neg, self.n_strata = len(self.pos_s), len(self.neg_s), len(uniq)

    def __call__(self, mult=None):
        wn = np.ones(self.n_neg) if mult is None else mult[self.neg_t].astype(float)
        wp = np.ones(self.n_pos) if mult is None else mult[self.pos_t].astype(float)
        cw = np.r_[0.0, np.cumsum(wn)]
        below = cw[self.lo] - cw[self.start_of_pos]
        tied = cw[self.hi] - cw[self.lo]
        total = cw[self.stop_of_pos] - cw[self.start_of_pos]
        den = (wp * total).sum()
        return float((wp * (below + 0.5 * tied)).sum() / den) if den > 0 else np.nan


def plain_auc(score, is_pos):
    return StratifiedAUC(score, is_pos, np.zeros(len(score), np.int64), np.zeros(len(score), np.int64))()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--draws', type=int, default=1000)
    ap.add_argument('--quick', action='store_true', help='first 6,000 half-steps, 50 draws (pipeline test only)')
    ap.add_argument('--same-phase-only', action='store_true',
                    help='post hoc (not in PREDICTIONS.md): only pairs whose members are locked to the same phase')
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    pred = (HERE / 'PREDICTIONS.md').read_bytes()
    pred_sha = hashlib.sha256(pred).hexdigest()
    log(f'PREDICTIONS.md sha256 {pred_sha}')
    draws = 50 if args.quick else args.draws

    with np.load(CHIP_DIR / 'states.npz') as z:
        states, nodes = z['states'], z['nodes']
    with np.load(CHIP_DIR / 'interventions.npz') as z:
        E0, int_idx = z['E0'], z['nodes']
    if args.quick:
        states = states[:6000]
    T, N = states.shape
    log(f'states {T} x {N}; intervened nodes {len(int_idx)}')

    chip = Chip(ROOT / 'data/visual6502')
    g_, c1_, c2_ = chip.transistor_matrix()
    A = tk.node_reference_from_transistors(g_, c1_, c2_, nodes, supply={chip.gnd, chip.pwr}).toarray() > 0

    # Universe: intervened nodes, in E0's row order
    U = np.asarray(int_idx)
    n = len(U)
    W = A[np.ix_(U, U)]
    infl = E0[:, U] > 0                      # infl[c, u]: clamping U[c] changes U[u] within the half-step
    np.fill_diagonal(infl, False)
    off = ~np.eye(n, dtype=bool)

    # common causes c != i, j of i and j (c ranges over U)
    Ci = infl.astype(np.int32)
    common = Ci.T @ Ci                       # counts c with infl[c,i] & infl[c,j]; c == i or j cannot occur (diag False)
    deg = Ci.sum(0)                          # |causes(u)|
    n_i = deg[:, None] - infl.T              # causes(i) without j  (infl.T[i, j] = infl[j, i])
    n_j = deg[None, :] - infl                # causes(j) without i
    union = n_i + n_j - common
    jacc = np.where(union > 0, common / np.maximum(union, 1), 0.0)

    direct = W & infl & off
    idle = W & ~infl & off
    codriven = ~W & ~W.T & ~infl & ~infl.T & (common > 0) & off
    log(f'pairs: direct {direct.sum()}, wired idle {idle.sum()}, co-driven only {codriven.sum()}')

    st = states[:, U].T.astype(np.uint8)                      # [n, T]
    chg = np.zeros_like(st); chg[:, 1:] = st[:, 1:] ^ st[:, :-1]
    parity = (np.arange(T) % 2).astype(np.uint8)
    ev = chg[:, 1:].astype(bool)
    rate = ev.mean(1)
    even_frac = np.where(ev.sum(1) > 0, ev[:, 1::2].sum(1) / np.maximum(ev.sum(1), 1), 0.5)  # samples t = 1, 3, ... -> index 1::2 is t even
    phase = np.where(even_frac >= LOCK, 0, np.where(even_frac <= 1 - LOCK, 1, 2))
    log(f'phase classes (even, odd, mixed): {np.bincount(phase, minlength=3).tolist()}')

    results, rows = dict(predictions_sha256=pred_sha, T=T, quick=args.quick, draws=draws), []
    rng = np.random.default_rng(SEED)
    boot_mult = rng.multinomial(n, np.full(n, 1 / n), size=draws)   # target-cluster multiplicities over nodes

    for rep, X in [('change', chg), ('state', st)]:
        t0 = time.perf_counter()
        reps, ug = tk.group_identical(X)
        G = len(reps); Xg = X[reps]
        ks = np.asarray([tk.ais_select_k(x, K_MAX) for x in Xg])
        scores = {
            'lag0_cmi': tk.te_matrix(Xg, ks, delay=0)[0],
            'lag0_cmi_parity': tk.te_matrix(Xg, ks, delay=0, cond=lambda t: parity[None, :])[0],
            'lag1_te': tk.te_matrix(Xg, ks, delay=1)[0],
        }
        log(f'{rep}: {G} groups, k {np.bincount(ks).tolist()}, scores in {time.perf_counter() - t0:.0f}s')
        valid = (ug[:, None] >= 0) & (ug[None, :] >= 0) & (ug[:, None] != ug[None, :]) & off
        if args.same_phase_only:
            valid &= (phase[:, None] == phase[None, :]) & (phase[:, None] < 2)
        ii, jj = np.nonzero(valid)
        cls = np.full(len(ii), '', dtype=object)
        cls[direct[ii, jj]] = 'direct'; cls[idle[ii, jj]] = 'wired idle'; cls[codriven[ii, jj]] = 'co-driven only'
        jq = np.zeros(len(ii), np.int64)
        m_cd = (cls == 'direct') | (cls == 'co-driven only')
        jq[m_cd] = quintile(jacc[ii, jj][m_cd])
        rq = quintile(rate)
        base_stratum = ((phase[ii] * 3 + phase[jj]) * 5 + rq[ii]) * 5 + rq[jj]
        stratum = base_stratum * 5 + jq
        pairs = pd.DataFrame(dict(cls=cls, jacc=jacc[ii, jj], i=ii, j=jj))
        res = {}
        for sname, S in scores.items():
            s = S[ug[ii], ug[jj]]
            out = {}
            for label, mpos, mneg, strat in [
                    ('P1 direct vs co-driven (stratified)', cls == 'direct', cls == 'co-driven only', stratum),
                    ('P3 direct vs wired idle (stratified, no co-drive stratum)', cls == 'direct', cls == 'wired idle', base_stratum)]:
                m = mpos | mneg
                if mpos.sum() == 0 or mneg.sum() == 0:
                    continue
                est = StratifiedAUC(s[m], mpos[m], strat[m], jj[m])
                point = est(); bs = np.array([est(b) for b in boot_mult])
                lo, hi = np.nanquantile(bs, [.025, .975])
                out[label] = dict(auc=point, ci=[float(lo), float(hi)], n_pos_used=est.n_pos, n_neg_used=est.n_neg,
                                  strata=est.n_strata, unstratified_auc=plain_auc(s[m], mpos[m]))
            cd = cls == 'co-driven only'
            out['P2 spearman(score, co-drive strength | co-driven only)'] = float(
                pd.Series(s[cd]).corr(pd.Series(jacc[ii, jj][cd]), method='spearman'))
            out['class_means_millibits'] = {c: float(1e3 * s[cls == c].mean()) for c in ['direct', 'wired idle', 'co-driven only']}
            out['unrelated_mean_millibits'] = float(1e3 * s[cls == ''].mean())
            res[sname] = out
            for label, v in out.items():
                if isinstance(v, dict) and 'auc' in v:
                    rows.append(dict(representation=rep, score=sname, contrast=label, AUC=v['auc'], lo=v['ci'][0], hi=v['ci'][1],
                                     unstratified=v['unstratified_auc'], positives=v['n_pos_used'], negatives=v['n_neg_used'], strata=v['strata']))
        res['pair_counts'] = pairs.cls.value_counts().to_dict()
        res['groups'] = int(G)
        results[rep] = res
        log(f'{rep}: done; P1 lag-0 CMI AUC {res["lag0_cmi"]["P1 direct vs co-driven (stratified)"]["auc"]:.3f} '
            f'{res["lag0_cmi"]["P1 direct vs co-driven (stratified)"]["ci"]}')

    tag = ('_quick' if args.quick else '') + ('_same_phase' if args.same_phase_only else '')
    (OUT / f'results{tag}.json').write_text(json.dumps(results, indent=2, default=float))
    pd.DataFrame(rows).to_csv(OUT / f'auc_table{tag}.csv', index=False)
    log(f'wrote {OUT}')


if __name__ == '__main__':
    main()
