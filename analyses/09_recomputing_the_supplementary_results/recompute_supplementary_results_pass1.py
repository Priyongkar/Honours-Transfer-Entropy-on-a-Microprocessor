"""Recompute the supplementary analyses whose output files are not in the project (4 November 2026).

The thesis cites supplementary files (supp_r3_baselines, supp_parity_*, supp_rq2_decomposition_*,
supp_rq2_enrichment_*, ...) that are not on this machine. This script recomputes them from the main
run's saved matrices, with the main run's definitions, so that every value in the thesis rests on a file
in the evidence package. It forms no new hypothesis; each value is compared with the thesis afterwards.

Blocks: B baselines (tab:r1-baselines), P parity and strata (tab:discrimination supplementary rows,
tab:r1-strata), D RQ2 decomposition (tab:rq2decomp, tab:rq2increments), E enrichment (tab:rq2enrich).
Usage: python recompute_supplementary.py [--selftest]
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / 'src'))
from te6502.paths import DATA_ROOT, RESULTS, load_analysis  # noqa: E402
import argparse, ast, glob, json, sys, time
from pathlib import Path
import numpy as np

ROOT = DATA_ROOT
import te6502.transfer_entropy_estimators as tk  # noqa: E402
from te6502.switch_level_simulator import load_netlist  # noqa: E402

RUN = ROOT / 'results/v5_bitcorrect_every_third_mi_ais_te_k1_12_L1_d1_ais_x1'
OUT = RESULTS / '09_recomputed_supplementary'; OUT.mkdir(exist_ok=True)
VSS, VCC = 558, 657
SEED, DRAWS = 2026, 1000
ALIGN = 1407


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', flush=True)


# ---------------------------------------------------------------- data and references
def load():
    X = np.asarray(np.load(RUN / 'thesis_v5/idtxl/data_int8.npy')).astype(np.uint8)
    te = np.load(RUN / 'te.npz'); K = te['optimal_k'].astype(int); TE1 = te['te_bits']
    TE0 = np.load(RUN / 'thesis_v5/te_lag0.npz')['te_bits']; MI = np.load(RUN / 'mi.npz')['mi_bits']
    p2g = np.load(RUN / 'groups.npz')['physical_to_group'].astype(int)
    NET = np.load(RUN / 'thesis_v5/idtxl/network_multivariate.npz')['adjacency'].astype(bool)
    for M in (TE1, TE0, MI):
        np.fill_diagonal(M, np.nan)
    return X, K, TE1, TE0, MI, p2g, NET


def netlist():
    defs = []
    for line in (ROOT / 'data/visual6502/transdefs.js').read_text().splitlines():
        if line.strip().startswith("['t"):
            defs.append(ast.literal_eval(line.strip().rstrip(',')))
    defs.sort(key=lambda d: int(d[0][1:]))
    gate, c1, c2 = (np.asarray([d[i] for d in defs]) for i in (1, 2, 3))
    box = np.asarray([d[4] for d in defs], float)
    centre = np.column_stack(((box[:, 0] + box[:, 1]) / 2, (box[:, 2] + box[:, 3]) / 2))
    return gate, c1, c2, centre


def references(p2g, G):
    gate, c1, c2, centre = netlist()
    gated_by = {}
    for t, g in enumerate(gate):
        gated_by.setdefault(int(g), []).append(t)
    T = len(gate)
    # transistor-level directed edges (terminal-to-gate), unique pairs
    edges = set()
    for s in range(T):
        for term in (c1[s], c2[s]):
            if term in (VSS, VCC): continue
            for t in gated_by.get(int(term), []):
                if t != s: edges.add((s, t))
    E = np.array(sorted(edges))
    act = p2g >= 0
    Wphys = np.zeros((G, G)); Wd = np.zeros((G, G), bool)
    for s, t in E:
        if act[s] and act[t] and p2g[s] != p2g[t]:
            Wphys[p2g[s], p2g[t]] += 1; Wd[p2g[s], p2g[t]] = True
    size = np.bincount(p2g[act], minlength=G).astype(float)
    allp = np.outer(size, size); np.fill_diagonal(allp, 0)
    # layout distance: minimum distance between transistor centres of the two groups
    from scipy.spatial.distance import cdist
    idx = np.flatnonzero(act); order = idx[np.argsort(p2g[idx])]
    bounds = np.r_[0, np.cumsum(np.bincount(p2g[order], minlength=G))]
    dist = np.zeros((G, G))
    for g in range(G):
        dmin = cdist(centre[order[bounds[g]:bounds[g + 1]]], centre[order]).min(0)
        dist[g] = np.minimum.reduceat(dmin, bounds[:-1])
    return dict(E=E, Wphys=Wphys, Nphys=allp - Wphys, Wd=Wd, size=size, dist=dist, gate=gate, c1=c1, c2=c2)


# ---------------------------------------------------------------- metrics with bootstrap
def group_auc(score, pos, mask, w_t=None, w_s=None):
    """Tie-aware AUC over the pairs in mask, each pair weighted by w_s[src]*w_t[tgt]."""
    ii, jj = np.nonzero(mask)
    s = score[ii, jj]; p = pos[ii, jj].astype(float)
    w = np.ones(len(ii)) if w_t is None else w_t[jj] * (1 if w_s is None else w_s[ii])
    return weighted_auc(s, p * w, (1 - p) * w)


def weighted_auc(s, pw, nw):
    keep = np.isfinite(s); s, pw, nw = s[keep], pw[keep], nw[keep]
    o = np.argsort(s, kind='mergesort'); s, pw, nw = s[o], pw[o], nw[o]
    _, idx = np.unique(s, return_index=True)
    P = np.add.reduceat(pw, idx); N = np.add.reduceat(nw, idx)
    below = np.cumsum(N) - N
    den = P.sum() * N.sum()
    return float((P * (below + N / 2)).sum() / den) if den > 0 else np.nan


def strat_auc(s, pw, nw, strata):
    num = den = 0.0
    for st in np.unique(strata):
        m = strata == st
        a = weighted_auc(s[m], pw[m], nw[m])
        d = pw[m].sum() * nw[m].sum()
        if d > 0 and np.isfinite(a):
            num += a * d; den += d
    return num / den if den > 0 else np.nan


def ap_lift(s, pw, nw):
    c = tk.weighted_curves(s, pw, nw)
    return float(c['ap_lift'])


def within_target_auc(score, pos, mask, w_t=None):
    G = score.shape[0]; vals = np.full(G, np.nan)
    for j in range(G):
        m = mask[:, j]
        if pos[m, j].any() and (~pos[m, j]).any():
            vals[j] = weighted_auc(score[m, j], pos[m, j].astype(float), (~pos[m, j]).astype(float))
    ok = np.isfinite(vals)
    if w_t is None:
        return float(vals[ok].mean())
    return float((vals[ok] * w_t[ok]).sum() / w_t[ok].sum())


def ci(f, mults):
    b = np.array([f(m) for m in mults]); b = b[np.isfinite(b)]
    return [float(np.quantile(b, .025)), float(np.quantile(b, .975))], b


def alloc_precision(score, pos, alloc, w_t=None, w_s=None, wphys=None, nphys=None):
    """Each target gets its alloc[j] top-ranked sources; group precision (pairs) and physical precision."""
    G = score.shape[0]; tp = n = 0.0; ptp = pall = 0.0
    for j in range(G):
        kj = int(alloc[j])
        if kj == 0: continue
        col = np.where(np.isfinite(score[:, j]), score[:, j], -np.inf); col[j] = -np.inf
        top = np.argsort(-col, kind='stable')[:kj]
        wt = 1.0 if w_t is None else w_t[j]
        ws = np.ones(len(top)) if w_s is None else w_s[top]
        tp += wt * (pos[top, j] * ws).sum(); n += wt * ws.sum()
        if wphys is not None:
            ptp += wt * (wphys[top, j] * ws).sum(); pall += wt * ((wphys[top, j] + nphys[top, j]) * ws).sum()
    return (100 * tp / n if n else np.nan), (100 * ptp / pall if pall else np.nan)


# ---------------------------------------------------------------- scores computed here
def parity_scores(X, K, sub=None):
    """Lag-1 TE | parity, lag-0 CMI | parity, MI | parity (zero-lag MI with k = 0)."""
    Xs = X if sub is None else X[:, sub:]
    G, T = Xs.shape; par = (np.arange(T) % 2).astype(np.uint8)[None, :]
    XS = tk.as_sparse_sources(Xs)
    out = {n: np.full((G, G), np.nan) for n in ('lag1_parity', 'lag0_parity', 'mi_parity', 'lag1', 'lag0', 'mi')}
    for t in range(G):
        out['lag1_parity'][:, t] = tk.te_all_sources(XS, Xs[t], K[t], delay=1, cond=par, cond_delay=0)[0]
        out['lag0_parity'][:, t] = tk.te_all_sources(XS, Xs[t], K[t], delay=0, cond=par, cond_delay=0)[0]
        out['mi_parity'][:, t] = tk.te_all_sources(XS, Xs[t], 0, delay=0, cond=par, cond_delay=0)[0]
        if sub is not None:
            out['lag1'][:, t] = tk.te_all_sources(XS, Xs[t], K[t], delay=1)[0]
            out['lag0'][:, t] = tk.te_all_sources(XS, Xs[t], K[t], delay=0)[0]
            out['mi'][:, t] = tk.te_all_sources(XS, Xs[t], 0, delay=0)[0]
    for M in out.values():
        np.fill_diagonal(M, np.nan)
    return out


def te_target_lags(X, y, lags):
    """Lag-1 TE with the target's history given by an arbitrary set of lags (IDTxl's selected target lags)."""
    lags = sorted(int(l) for l in lags); y = np.asarray(y).astype(np.int64); T = len(y)
    start = max(max(lags), 1); n = T - start
    code = np.zeros(n, dtype=np.int64)
    for l in lags:
        code = (code << 1) | y[start - l:T - l]
    _, state = np.unique(code, return_inverse=True); state = state.reshape(-1); S = int(state.max()) + 1
    from scipy import sparse
    j = state * 2 + y[start:T]
    M = sparse.csr_matrix((np.ones(n), (np.arange(start, T) - 1, j)), shape=(T, 2 * S))
    XS = tk.as_sparse_sources(X)
    C1 = np.asarray((XS @ M).todense()).reshape(XS.shape[0], S, 2)
    Ntot = np.bincount(j, minlength=2 * S).astype(float).reshape(S, 2); C0 = Ntot[None] - C1; ns = Ntot.sum(1)
    outv = np.zeros(C1.shape[0])
    for C in (C0, C1):
        nsx = C.sum(2, keepdims=True)
        with np.errstate(divide='ignore', invalid='ignore'):
            outv += np.where(C > 0, C * np.log2((C * ns[None, :, None]) / (nsx * Ntot[None])), 0).sum((1, 2))
    return np.maximum(outv / n, 0)


def phase_classes(X):
    n = X.sum(1).astype(float); even = X[:, 0::2].sum(1) / np.maximum(n, 1)
    return np.where(even >= .95, 0, np.where(even <= .05, 1, 2))


def deciles(v):
    return np.searchsorted(np.unique(np.quantile(v, np.linspace(.1, .9, 9))), v, side='right')


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--selftest', action='store_true'); a = ap.parse_args()
    X, K, TE1, TE0, MI, p2g, NET = load(); G, T = X.shape
    log(f'loaded: {G} groups, {T} samples')
    ref = references(p2g, G); Wd, Wphys, Nphys = ref['Wd'], ref['Wphys'], ref['Nphys']
    size, dist = ref['size'], ref['dist']
    off = ~np.eye(G, dtype=bool)
    rng = np.random.default_rng(SEED)
    draws = 50 if a.selftest else DRAWS
    mt = rng.multinomial(G, np.full(G, 1 / G), size=draws)          # target-cluster
    ms = np.random.default_rng(SEED + 1).multinomial(G, np.full(G, 1 / G), size=draws)  # second axis for two-way
    res = dict(n_groups=G, n_samples=T, wired_group_pairs=int(Wd.sum()), wired_transistor_pairs=int(Wphys.sum()))
    log(f"references: {res['wired_group_pairs']} wired group pairs, {res['wired_transistor_pairs']} wired transistor pairs")

    # ---- check: the main run's headline group AUCs are reproduced by this code
    chk = {}
    for name, S in (('lag1_TE', TE1), ('lag0_CMI', TE0), ('MI', MI)):
        g = group_auc(S, Wd, off); c, _ = ci(lambda m: group_auc(S, Wd, off, m.astype(float)), mt)
        ii, jj = np.nonzero(off)
        ph = weighted_auc(S[ii, jj], Wphys[ii, jj], Nphys[ii, jj])
        chk[name] = dict(group_auc=g, group_ci=c, physical_auc=ph,
                         within_target=within_target_auc(S, Wd, off))
    res['check_against_main_run'] = chk
    log(f'check: {json.dumps({k: round(v["group_auc"], 3) for k, v in chk.items()})} (main run: lag1 0.469, lag0 0.632, MI 0.597)')
    if a.selftest:
        (OUT / 'selftest.json').write_text(json.dumps(res, indent=1)); log('SELF-TEST DONE'); return

    rate = X.mean(1); clock_g = int(np.argmax(rate))
    alloc = NET.sum(0)
    scores = {'source_group_size': np.repeat(size[:, None], G, 1), 'size_product': np.outer(size, size),
              'layout_distance': -dist, 'lag0_CMI': TE0, 'zero_lag_MI': MI, 'lag1_TE': TE1}
    ii, jj = np.nonzero(off)

    def block(mask, keep_alloc_src=None):
        out = {}
        mi_, mj_ = np.nonzero(mask)
        for name, S in scores.items():
            s = S[mi_, mj_]; p = Wd[mi_, mj_].astype(float)
            r = dict(group_auc=weighted_auc(s, p, 1 - p))
            r['group_auc_ci'], _ = ci(lambda m: weighted_auc(s, p * m[mj_], (1 - p) * m[mj_]), mt)
            r['group_ap_lift'] = ap_lift(s, p, 1 - p)
            r['group_ap_lift_ci'], _ = ci(lambda m: ap_lift(s, p * m[mj_], (1 - p) * m[mj_]), mt[:200])
            r['within_target_auc'] = within_target_auc(S, Wd, mask)
            r['within_target_ci'], _ = ci(lambda m: within_target_auc(S, Wd, mask, m.astype(float)), mt[:200])
            r['physical_auc'] = weighted_auc(s, Wphys[mi_, mj_], Nphys[mi_, mj_])
            r['physical_auc_ci'], _ = ci(lambda m: weighted_auc(s, Wphys[mi_, mj_] * m[mj_], Nphys[mi_, mj_] * m[mj_]), mt[:200])
            Sa = S.copy()
            if keep_alloc_src is not None:
                Sa[~keep_alloc_src, :] = np.nan
            r['precision_at_allocation'], r['physical_precision_at_allocation'] = alloc_precision(Sa, Wd, alloc, wphys=Wphys, nphys=Nphys)
            r['precision_at_allocation_ci'], _ = ci(lambda m: alloc_precision(Sa, Wd, alloc, w_t=m.astype(float))[0], mt[:200])
            out[name] = r
            log(f'  {name}: group AUC {r["group_auc"]:.3f}, precision at allocation {r["precision_at_allocation"]:.2f}')
        return out

    # ---- B. baselines
    log('B: baselines, all groups')
    res['B_baselines_all'] = block(off)
    noclk = off.copy(); noclk[clock_g, :] = False; noclk[:, clock_g] = False
    srcok = np.ones(G, bool); srcok[clock_g] = False
    log('B: baselines, clock group excluded')
    res['B_baselines_clock_excluded'] = block(noclk, keep_alloc_src=srcok)
    # panel B: within covariate deciles
    pb = {}
    for cname, cov in (('source_group_size', np.repeat(size[:, None], G, 1)), ('size_product', np.outer(size, size)),
                       ('layout_distance', dist)):
        strata = deciles(cov[ii, jj]); p = Wd[ii, jj].astype(float)
        r = {}
        bs = {}
        for sname, S in (('lag0', TE0), ('lag1', TE1)):
            s = S[ii, jj]
            r[f'{sname}_group'] = strat_auc(s, p, 1 - p, strata)
            r[f'{sname}_group_ci'], bs[sname] = ci(lambda m: strat_auc(s, p * m[jj], (1 - p) * m[jj], strata), mt[:200])
            bs[sname + '_phys'] = np.array([strat_auc(s, Wphys[ii, jj] * m[jj], Nphys[ii, jj] * m[jj], strata) for m in mt[:200]])
            r[f'{sname}_physical'] = strat_auc(s, Wphys[ii, jj], Nphys[ii, jj], strata)
        d = bs['lag0'] - bs['lag1']; dp = bs['lag0_phys'] - bs['lag1_phys']
        r['delta_group'] = r['lag0_group'] - r['lag1_group']; r['delta_group_ci'] = [float(np.quantile(d, .025)), float(np.quantile(d, .975))]
        r['delta_physical'] = r['lag0_physical'] - r['lag1_physical']; r['delta_physical_ci'] = [float(np.quantile(dp, .025)), float(np.quantile(dp, .975))]
        pb[cname] = r; log(f'B panel B {cname}: lag0 {r["lag0_group"]:.3f}, lag1 {r["lag1_group"]:.3f}')
    res['B_within_covariate_deciles'] = pb
    # panel C: covariate logistic model and lag-0 increment
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold
    ph = phase_classes(X)
    pc = np.zeros((len(ii), 3)); pc[np.arange(len(ii)), np.where(ph[ii] == ph[jj], 0, np.where((ph[ii] == 2) | (ph[jj] == 2), 2, 1))] = 1
    base_feats = np.column_stack([np.log(size[ii]), np.log(size[jj]), np.log1p(dist[ii, jj]),
                                  np.log(np.maximum(rate[ii], 1e-6)), np.log(np.maximum(rate[jj], 1e-6)), pc])
    y = Wd[ii, jj].astype(int)
    pc_res = {}
    for fit in ('group', 'physical'):
        w = np.ones(len(ii)) if fit == 'group' else Wphys[ii, jj] + Nphys[ii, jj]
        preds = {}
        for name, F in (('covariates', base_feats), ('covariates_plus_lag0', np.column_stack([base_feats, np.log1p(1e4 * np.nan_to_num(TE0[ii, jj]))]))):
            pr = np.zeros(len(ii))
            for tr, te_ in GroupKFold(5).split(F, y, groups=jj):
                mu, sd = F[tr].mean(0), F[tr].std(0) + 1e-12
                m = LogisticRegression(max_iter=2000).fit((F[tr] - mu) / sd, y[tr], sample_weight=w[tr] / w[tr].mean())
                pr[te_] = m.decision_function((F[te_] - mu) / sd)
            preds[name] = pr
        r = {}
        if fit == 'group':
            met = lambda s, m: weighted_auc(s, y * m[jj], (1 - y) * m[jj])
            metap = lambda s, m: ap_lift(s, y * m[jj], (1 - y) * m[jj])
        else:
            met = lambda s, m: weighted_auc(s, Wphys[ii, jj] * m[jj], Nphys[ii, jj] * m[jj])
            metap = lambda s, m: ap_lift(s, Wphys[ii, jj] * m[jj], Nphys[ii, jj] * m[jj])
        one = np.ones(G)
        for mname, mf in (('auc', met), ('ap_lift', metap)):
            a0 = mf(preds['covariates'], one); a1 = mf(preds['covariates_plus_lag0'], one)
            b = np.array([mf(preds['covariates_plus_lag0'], m) - mf(preds['covariates'], m) for m in mt[:200]])
            c0, _ = ci(lambda m: mf(preds['covariates'], m), mt[:200])
            r[mname] = dict(covariates=a0, covariates_ci=c0, increment=a1 - a0, increment_ci=[float(np.quantile(b, .025)), float(np.quantile(b, .975))])
        pc_res[f'{fit}_fitted'] = r; log(f'B panel C {fit}: AUC {r["auc"]["covariates"]:.4f}, increment {r["auc"]["increment"]:+.4f}')
    res['B_covariate_model'] = pc_res
    (OUT / 'recomputed_supplementary.json').write_text(json.dumps(res, indent=1, default=float))

    # ---- P. parity conditioning and strata
    log('P: parity-conditioned scores, whole recording')
    pw = parity_scores(X, K)
    log('P: clock-aligned recording')
    pa = parity_scores(X, K, sub=ALIGN)
    disc = {}
    for label, d, extra in (('whole', pw, dict(lag1=TE1, lag0=TE0, mi=MI)), ('clock_aligned', pa, {})):
        mats = dict(d); mats.update(extra if extra else {})
        for name in ('lag0_parity', 'mi_parity', 'lag1_parity', 'lag1', 'lag0', 'mi'):
            if name not in mats or np.all(np.isnan(mats[name])): continue
            S = mats[name]; s = S[ii, jj]; p = Wd[ii, jj].astype(float)
            r = dict(group_auc=weighted_auc(s, p, 1 - p), group_ap_lift=ap_lift(s, p, 1 - p),
                     physical_auc=weighted_auc(s, Wphys[ii, jj], Nphys[ii, jj]),
                     physical_ap_lift=ap_lift(s, Wphys[ii, jj], Nphys[ii, jj]),
                     within_target=within_target_auc(S, Wd, off))
            r['group_auc_ci'], _ = ci(lambda m: weighted_auc(s, p * m[jj], (1 - p) * m[jj]), mt[:200])
            disc[f'{label}:{name}'] = r
        if label == 'whole':
            S = -TE1; s = S[ii, jj]; p = Wd[ii, jj].astype(float)
            disc['whole:lag1_reversed'] = dict(group_auc=weighted_auc(s, p, 1 - p), physical_auc=weighted_auc(s, Wphys[ii, jj], Nphys[ii, jj]),
                                               within_target=within_target_auc(S, Wd, off), group_ap_lift=ap_lift(s, p, 1 - p))
            cls = np.where(ph[ii] == ph[jj], 2, np.where((ph[ii] == 2) | (ph[jj] == 2), 1, 0)).astype(float)
            disc['whole:phase_class_only'] = dict(group_auc=weighted_auc(cls, p, 1 - p), physical_auc=weighted_auc(cls, Wphys[ii, jj], Nphys[ii, jj]),
                                                  group_ap_lift=ap_lift(cls, p, 1 - p))
            combo = cls * 10 + np.argsort(np.argsort(np.nan_to_num(TE0[ii, jj]))) / len(ii)
            disc['whole:phase_class_then_lag0'] = dict(group_auc=weighted_auc(combo, p, 1 - p), group_ap_lift=ap_lift(combo, p, 1 - p),
                                                       physical_auc=weighted_auc(combo, Wphys[ii, jj], Nphys[ii, jj]))
    res['P_discrimination_supplementary'] = disc
    # strata (tab:r1-strata)
    strata_res = {}
    same = ph[ii] == ph[jj]; mixed = (ph[ii] == 2) | (ph[jj] == 2); oppo = ~same & ~mixed
    same = same & ~mixed
    notclk = (ii != clock_g) & (jj != clock_g)
    for sname, m in (('same_parity', same), ('opposite_parity', oppo), ('involves_mixed', mixed), ('involves_mixed_no_clock', mixed & notclk)):
        r = {'wired_group_pairs': int(Wd[ii[m], jj[m]].sum())}
        for scn, S in (('lag1', TE1), ('lag0', TE0)):
            s = S[ii[m], jj[m]]; p = Wd[ii[m], jj[m]].astype(float)
            r[f'{scn}_group'] = weighted_auc(s, p, 1 - p)
            r[f'{scn}_group_ci'], _ = ci(lambda mm: weighted_auc(s, p * mm[jj[m]], (1 - p) * mm[jj[m]]), mt[:200])
            r[f'{scn}_physical'] = weighted_auc(s, Wphys[ii[m], jj[m]], Nphys[ii[m], jj[m]])
        strata_res[sname] = r
    ent = np.array([tk.ais_select_k(x, 1) for x in X[:1]])  # placeholder to keep tk import used
    srcprops = {'source_rate': rate}
    for name, v in srcprops.items():
        strata = deciles(v[ii]); r = {}
        for scn, S in (('lag1', TE1), ('lag0', TE0)):
            s = S[ii, jj]; p = Wd[ii, jj].astype(float)
            r[f'{scn}_group'] = strat_auc(s, p, 1 - p, strata)
            r[f'{scn}_physical'] = strat_auc(s, Wphys[ii, jj], Nphys[ii, jj], strata)
        strata_res[f'within_deciles_{name}'] = r
    rr = deciles(rate[ii]) * 10 + deciles(rate[jj]); r = {}
    for scn, S in (('lag1', TE1), ('lag0', TE0)):
        s = S[ii, jj]; p = Wd[ii, jj].astype(float)
        r[f'{scn}_group'] = strat_auc(s, p, 1 - p, rr); r[f'{scn}_physical'] = strat_auc(s, Wphys[ii, jj], Nphys[ii, jj], rr)
    strata_res['within_source_rate_x_target_rate'] = r
    res['P_strata'] = strata_res
    (OUT / 'recomputed_supplementary.json').write_text(json.dumps(res, indent=1, default=float))
    log('P done')

    # ---- D. RQ2 decomposition at IDTxl's allocation
    log('D: lag-1 TE at IDTxl target lags')
    TEemb = np.full((G, G), np.nan); XS = tk.as_sparse_sources(X)
    for f in sorted(glob.glob(str(RUN / 'thesis_v5/idtxl/multivariate/target_*.json'))):
        d = json.load(open(f)); j = int(d['target'])
        lags = [l for (_, l) in d.get('selected_vars_target', [])] or [1]
        TEemb[:, j] = te_target_lags(XS, X[j], lags)
    np.fill_diagonal(TEemb, np.nan)
    rates_model = np.log(np.maximum(rate, 1e-6))[:, None] + 0 * np.log(np.maximum(rate, 1e-6))[None, :]
    nets = {'lag1_global_cut': None, 'lag1_allocation_matched': TE1, 'lag1_at_idtxl_target_lags': TEemb,
            'lag0_allocation_matched': TE0, 'lag0_parity_allocation_matched': pw['lag0_parity'], 'mi_allocation_matched': MI,
            'source_event_rate_allocation_matched': rates_model,
            'source_group_size_allocation_matched': np.repeat(size[:, None], G, 1), 'layout_distance_allocation_matched': -dist}
    dres = {}
    ntot = int(NET.sum())
    s1 = TE1[ii, jj]; order = np.argsort(-np.nan_to_num(s1, nan=-np.inf), kind='stable')[:ntot]
    gc = np.zeros((G, G), bool); gc[ii[order], jj[order]] = True
    def net_prec(A, wt=None, ws=None):
        wt = np.ones(G) if wt is None else wt; ws = np.ones(G) if ws is None else ws
        W = np.outer(ws, wt)
        return 100 * (A & Wd).astype(float).__mul__(W).sum() / (A.astype(float) * W).sum()
    def phys_prec(A):
        return 100 * Wphys[A].sum() / (Wphys[A] + Nphys[A]).sum()
    netsets = {'idtxl_multivariate': NET, 'lag1_global_cut': gc}
    for name, S in nets.items():
        if S is None: continue
        A = np.zeros((G, G), bool)
        for j in range(G):
            kj = int(alloc[j])
            if kj == 0: continue
            col = np.where(np.isfinite(S[:, j]), S[:, j], -np.inf); col[j] = -np.inf
            A[np.argsort(-col, kind='stable')[:kj], j] = True
        netsets[name] = A
    boots = {}
    for name, A in netsets.items():
        p_ = net_prec(A); c, b = ci(lambda m: net_prec(A, m.astype(float)), mt)
        b2 = np.array([net_prec(A, mt[i].astype(float), ms[i].astype(float)) for i in range(len(mt))])
        boots[name] = (b, b2)
        dres[name] = dict(edges=int(A.sum()), wired=int((A & Wd).sum()), group_precision=p_, group_ci=c,
                          physical_precision=phys_prec(A), transistor_pairs=int((Wphys[A] + Nphys[A]).sum()))
        log(f'D {name}: {dres[name]["wired"]}/{dres[name]["edges"]} wired, {p_:.2f} %')
    inc = {}
    for lab, a_, b_ in (('total_idtxl_minus_global_cut', 'idtxl_multivariate', 'lag1_global_cut'),
                        ('allocation', 'lag1_allocation_matched', 'lag1_global_cut'),
                        ('embedding', 'lag1_at_idtxl_target_lags', 'lag1_allocation_matched'),
                        ('conditioning', 'idtxl_multivariate', 'lag1_at_idtxl_target_lags'),
                        ('lag0_minus_idtxl', 'lag0_allocation_matched', 'idtxl_multivariate')):
        d1 = boots[a_][0] - boots[b_][0]; d2 = boots[a_][1] - boots[b_][1]
        inc[lab] = dict(increment=dres[a_]['group_precision'] - dres[b_]['group_precision'],
                        target_ci=[float(np.quantile(d1, .025)), float(np.quantile(d1, .975))], target_positive=float((d1 > 0).mean()),
                        two_way_ci=[float(np.quantile(d2, .025)), float(np.quantile(d2, .975))], two_way_positive=float((d2 > 0).mean()))
    res['D_rq2_networks'] = dres; res['D_rq2_increments'] = inc
    (OUT / 'recomputed_supplementary.json').write_text(json.dumps(res, indent=1, default=float))
    log('D done')

    # ---- E. enrichment of relation classes (transistor-path and group-graph definitions)
    log('E: relation classes')
    gate, c1, c2, E = ref['gate'], ref['c1'], ref['c2'], ref['E']
    from scipy import sparse
    nT = len(gate)
    At = sparse.csr_matrix((np.ones(len(E)), (E[:, 0], E[:, 1])), shape=(nT, nT))
    act = p2g >= 0
    P = sparse.csr_matrix((np.ones(act.sum()), (np.flatnonzero(act), p2g[act])), shape=(nT, G))
    def grp(M):
        R = (P.T @ (M > 0).astype(float) @ P).toarray() > 0; np.fill_diagonal(R, False); return R
    wired_t = grp(At); rev_t = grp(At.T)
    two_t = grp(At @ At); sp_t = grp(At.T @ At); sc_t = grp(At @ At.T)
    Ag = Wd.astype(float)
    two_g = (Ag @ Ag) > 0; sp_g = (Ag.T @ Ag) > 0; sc_g = (Ag @ Ag.T) > 0
    for M in (two_g, sp_g, sc_g): np.fill_diagonal(M, False)
    def classes(w, r, two, sp, sc):
        c = np.full((G, G), 5)
        for k, M in ((4, sc), (3, sp), (2, two), (1, r & ~w), (0, w)):
            c[M] = k
        np.fill_diagonal(c, -1); return c
    defs = {'transistor_paths': classes(wired_t, rev_t, two_t, sp_t, sc_t), 'group_graph': classes(Wd, Wd.T, two_g, sp_g, sc_g)}
    names = ['wired', 'reverse_only', 'two_step', 'shared_parent', 'shared_child', 'unrelated']
    from scipy.stats import binomtest
    A_idtxl = NET; A_l1 = netsets['lag1_allocation_matched']
    enr = {}
    for dname, C in defs.items():
        out = {}
        def rate_(A, k, wt=None, ws=None):
            wt = np.ones(G) if wt is None else wt; ws = np.ones(G) if ws is None else ws
            W = np.outer(ws, wt); m = C == k
            return (A[m] * W[m]).sum() / max(W[m].sum(), 1e-12)
        for k, nm in enumerate(names):
            m = C == k
            r = dict(class_size=int(m.sum()), declared_idtxl=int((A_idtxl & m).sum()), declared_lag1=int((A_l1 & m).sum()),
                     discordant_idtxl_only=int((A_idtxl & ~A_l1 & m).sum()), discordant_lag1_only=int((A_l1 & ~A_idtxl & m).sum()))
            if nm != 'unrelated' and m.sum():
                ratio = lambda wt=None, ws=None: (rate_(A_idtxl, k, wt, ws) / rate_(A_idtxl, 5, wt, ws)) / max(rate_(A_l1, k, wt, ws) / rate_(A_l1, 5, wt, ws), 1e-12)
                r['ratio'] = ratio()
                bt = np.array([ratio(mt[i].astype(float)) for i in range(200)]); bw = np.array([ratio(mt[i].astype(float), ms[i].astype(float)) for i in range(200)])
                r['target_ci'] = [float(np.nanquantile(bt, .025)), float(np.nanquantile(bt, .975))]
                r['two_way_ci'] = [float(np.nanquantile(bw, .025)), float(np.nanquantile(bw, .975))]
            n1, n2 = r['discordant_idtxl_only'], r['discordant_lag1_only']
            r['exact_p'] = float(binomtest(n1, n1 + n2, 0.5).pvalue) if n1 + n2 > 0 else 1.0
            out[nm] = r
        ps = sorted((out[nm]['exact_p'], nm) for nm in names[1:5])
        for i, (pv, nm) in enumerate(ps):
            out[nm]['holm_p'] = float(min(1, max((len(ps) * 2 - j) * ps[j][0] for j in range(i + 1))))
        enr[dname] = out
        log(f"E {dname}: shared parent {out['shared_parent']['declared_idtxl']}/{out['shared_parent']['declared_lag1']}")
    res['E_enrichment'] = enr
    (OUT / 'recomputed_supplementary.json').write_text(json.dumps(res, indent=1, default=float))
    log('ALL DONE')


if __name__ == '__main__':
    main()
