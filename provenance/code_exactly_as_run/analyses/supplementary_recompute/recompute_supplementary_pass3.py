"""Third recomputation pass (4 November 2026): the verdict-bearing intervals of the decision register.

Recomputes, from the main run's saved matrices and with the main run's definitions (notebook cells 47, 50,
52 and 54 of Honours_Every_Third_MI_AIS_TE_v5.ipynb), the values printed in
  1. tab:m-schemes   H1 to H4 under the target-cluster, two-way and vertex bootstraps, and the Holm p-values;
  2. tab:app-h3-clock (H3 without the clock group's pairs) and tab:r1-h4 (H4 within parity strata);
  3. tab:r1-h2 robustness rows: clock group excluded, common history k = 1, node level.

Definitions (thesis methods, Sections "Uncertainty" and "Quantities used by the hypotheses"):
  * AUC: tie-aware, group level (each ordered group pair once, wired = any evaluable edge) or physical level
    (each group pair weighted by its wired / unwired transistor pairs).
  * Orientation accuracy: among unidirectionally wired group pairs (wired i->j, not j->i), the weighted
    fraction of untied pairs for which the orientation with the larger score is the wired one. Group level:
    weight 1; physical: number of physical edges. Rate rule: the member with the lower event rate is the target.
  * Signed influence delta: saved in thesis_v5/signed_effect.npz (cell 47). H4 statistic: percentage of pairs
    with delta < 0, pull-down-only minus unwired; group level unweighted, physical level weighted by the
    pull-down transistor pairs (pull-down class) and unwired transistor pairs (unwired class).
  * Bootstraps (1,000 draws): a pair (i, j) gets weight a[i] * b[j] * (its level weight).
      target_main : b = bincount(rng.integers(0, G, G)), rng = default_rng(2026), a = 1  (tk.bootstrap_targets,
                    the generator of the main run; reproduces the main run's CSV intervals exactly)
      target      : b = multinomial draws, seed 2026, a = 1   (generator of recompute_supplementary.py)
      two_way     : b = multinomial seed 2026, a = independent multinomial seed 2027   (Owen 2007)
      vertex      : a = b = the same multinomial seed 2026 vector                     (Snijders & Borgatti 1999)
  * Holm p: min{1, 2(c + 1)/(B + 1)}, c = draws at or beyond the null (on the side opposite the estimate),
    deciding level, target-cluster draws; Holm over H1..H6; H3: max over its four histories (intersection-union).

Usage: python recompute_supplementary_pass3.py [--selftest]
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import recompute_supplementary as R            # noqa: E402  (load(), netlist(), phase_classes(), weighted_auc)
tk = R.tk
RUN, OUT = R.RUN, HERE / 'output'
CACHE = OUT / 'cache_pass3'; CACHE.mkdir(parents=True, exist_ok=True)
B = 1000; SEED = 2026
LOGF = None


def log(m):
    s = f'[{time.strftime("%H:%M:%S")}] {m}'
    print(s, flush=True)
    if LOGF is not None:
        LOGF.write(s + '\n'); LOGF.flush()


# ---------------------------------------------------------------- structural reference (notebook cell 40)
def reference(p2g, G):
    gate, c1, c2, _ = R.netlist()
    M = len(gate); lookup = {}
    for t, g in enumerate(gate):
        lookup.setdefault(int(g), []).append(t)
    src, dst, pd = [], [], []
    for s in range(M):
        dests = set()
        for term in (c1[s], c2[s]):
            if term not in (R.VSS, R.VCC):
                dests.update(lookup.get(int(term), []))
        dests.discard(s)
        for t in dests:
            src.append(s); dst.append(t); pd.append(R.VSS in (c1[s], c2[s]))
    src, dst, pd = map(np.asarray, (src, dst, pd))
    act = np.flatnonzero(p2g >= 0)
    P = sparse.csr_matrix((np.ones(len(act)), (act, p2g[act])), shape=(M, G))

    def agg(mask):
        A = sparse.csr_matrix((np.ones(mask.sum()), (src[mask], dst[mask])), shape=(M, M))
        C = np.asarray((P.T @ A @ P).todense()); np.fill_diagonal(C, 0); return C
    POS = agg(np.ones(len(src), bool)); PD = agg(pd); PASS = agg(~pd)
    size = np.bincount(p2g[act], minlength=G).astype(float)
    PAIRS = np.outer(size, size); np.fill_diagonal(PAIRS, 0)
    return dict(POS=POS, PD=PD, PASS=PASS, NEG=PAIRS - POS, size=size, gate=gate, c1=c1, c2=c2)


# ---------------------------------------------------------------- draws
def make_draws(G):
    rng = np.random.default_rng(SEED)
    tmain = np.array([np.bincount(rng.integers(0, G, G), minlength=G) for _ in range(B)], float)
    t = np.random.default_rng(SEED).multinomial(G, np.full(G, 1 / G), size=B).astype(float)
    s = np.random.default_rng(SEED + 1).multinomial(G, np.full(G, 1 / G), size=B).astype(float)
    one = np.ones((B, G))
    return {'target_main': (one, tmain), 'target': (one, t), 'two_way': (s, t), 'vertex': (t, t)}


def q(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return [float(np.quantile(v, .025)), float(np.quantile(v, .975))]


def boot_p(est, draws, null):
    """Conservative two-sided bootstrap p: 2(c+1)/(B+1), c = draws at or beyond the null."""
    d = np.asarray(draws, float); d = d[np.isfinite(d)]
    c = int((d <= null).sum()) if est > null else int((d >= null).sum())
    return min(1.0, 2 * (c + 1) / (len(d) + 1)), c


# ---------------------------------------------------------------- AUC engine (pre-sorted; weights change, order not)
class AUC:
    def __init__(self, S, pos, neg, mask):
        ii, jj = np.nonzero(mask)
        s = S[ii, jj]; p = pos[ii, jj].astype(float); n = neg[ii, jj].astype(float)
        keep = np.isfinite(s) & ((p + n) > 0)
        o = np.argsort(-s[keep], kind='stable')
        self.s, self.p, self.n = s[keep][o], p[keep][o], n[keep][o]
        self.i, self.j = ii[keep][o], jj[keep][o]
        self.ends = np.r_[np.flatnonzero(np.diff(self.s) != 0), len(self.s) - 1]

    def __call__(self, a=None, b=None):
        w = np.ones(len(self.s))
        if a is not None: w = w * a[self.i]
        if b is not None: w = w * b[self.j]
        tp = np.cumsum(self.p * w)[self.ends]; fp = np.cumsum(self.n * w)[self.ends]
        if tp[-1] <= 0 or fp[-1] <= 0: return np.nan
        tpr = np.r_[0, tp / tp[-1]]; fpr = np.r_[0, fp / fp[-1]]
        return float(np.sum(np.diff(fpr) * (tpr[1:] + tpr[:-1]) / 2))


def auc_pair(SA, SB, pos, neg, mask, draws):
    """AUCs of two scores and their difference (A - B), point and per scheme draws."""
    ea, eb = AUC(SA, pos, neg, mask), AUC(SB, pos, neg, mask)
    out = dict(a=ea(), b=eb()); out['diff'] = out['a'] - out['b']; bd = {}
    for sch, (A, Bm) in draws.items():
        va = np.array([ea(A[d], Bm[d]) for d in range(B)]); vb = np.array([eb(A[d], Bm[d]) for d in range(B)])
        bd[sch] = dict(a=va, b=vb, diff=va - vb)
    return out, bd


# ---------------------------------------------------------------- orientation accuracy (cell 52)
def oriented(fwd, rev, w):
    correct = fwd > rev; tie = np.isclose(fwd, rev, rtol=0, atol=1e-15)
    return (w * (correct & ~tie)).sum() / max((w * ~tie).sum(), 1e-300)


def orientation_block(S, rate, I, J, W, draws):
    """Accuracy of S's asymmetry, of the rate rule, and the paired difference, point + draws per scheme."""
    f, r = S[I, J], S[J, I]; fr, rr = rate[I], rate[J]       # rule: source has the higher rate
    pt = dict(acc=oriented(f, r, W), rule=oriented(fr, rr, W)); pt['diff'] = pt['acc'] - pt['rule']
    bd = {}
    for sch, (A, Bm) in draws.items():
        acc = np.empty(B); rule = np.empty(B)
        for d in range(B):
            w = W * A[d][I] * Bm[d][J]
            acc[d] = oriented(f, r, w); rule[d] = oriented(fr, rr, w)
        bd[sch] = dict(acc=acc, rule=rule, diff=acc - rule)
    return pt, bd


# ---------------------------------------------------------------- H4 percentage negative
def pct_neg(neg, w):
    return 100 * (w * neg).sum() / w.sum() if w.sum() > 0 else np.nan


def summarise(pt, bd, keys):
    out = {}
    for k in keys:
        r = dict(estimate=float(pt[k]))
        for sch, v in bd.items():
            r[sch] = q(v[k])
        for sch in ('two_way', 'vertex', 'target'):
            r.setdefault(sch, None)
        out[k] = r
    return out


# ---------------------------------------------------------------- main
def main():
    global LOGF
    ap = argparse.ArgumentParser(); ap.add_argument('--selftest', action='store_true'); a = ap.parse_args()
    LOGF = open(OUT / ('recompute_pass3_selftest.log' if a.selftest else 'recompute_pass3.log'), 'w')
    X, K, TE1, TE0, MI, p2g, NET = R.load(); G, T = X.shape
    log(f'loaded {G} groups x {T} samples')
    ref = reference(p2g, G); POS, PD, PASS, NEG = ref['POS'], ref['PD'], ref['PASS'], ref['NEG']
    sc = np.load(RUN / 'analysis_pairwise/structural_counts.npz')
    assert np.array_equal(sc['positive_count'], POS), 'reference differs from the main run'
    WIRED = (POS > 0).astype(float); UNW = (POS == 0).astype(float)
    off = ~np.eye(G, dtype=bool)
    rate = X.mean(1); clock = int(np.argmax(rate))
    DELTA = np.load(RUN / 'thesis_v5/signed_effect.npz')['delta']
    log(f'reference matches structural_counts.npz: {int(POS.sum())} edges, {int((POS > 0).sum())} wired group pairs; clock group {clock}')
    res = dict(meta=dict(draws=B, seed_target=SEED, seed_second_axis=SEED + 1, groups=G, samples=T, clock_group=clock,
                         schemes=dict(target_main='tk.bootstrap_targets generator (main run), target side only',
                                      target='multinomial seed 2026, target side only',
                                      two_way='targets multinomial seed 2026, sources independent multinomial seed 2027',
                                      vertex='one multinomial seed 2026 vector applied to sources and targets')))
    # sanity: the estimator reproduces the saved lag-1 TE for three targets
    XS = tk.as_sparse_sources(X)
    for t in (0, 6, 300):
        v = tk.te_all_sources(XS, X[t], K[t])[0]; v[t] = np.nan
        assert np.nanmax(np.abs(v - TE1[:, t])) < 1e-12
    draws = make_draws(G)
    if a.selftest:                                    # self-test: point values and the main run's own draws only
        draws = {'target_main': draws['target_main']}
    check = {}

    # ---------- common-history TE (cell 52: tk.te_matrix(RQ_X, np.full(G, k))) and lag-0 at k = 1
    common = {}
    for k in (1, 2, 4, 8):
        f = CACHE / f'te_lag1_common_k{k}.npy'
        if f.exists():
            common[k] = np.load(f)
        else:
            t0 = time.time(); common[k] = tk.te_matrix(X, np.full(G, k))[0]; np.save(f, common[k])
            log(f'common-history lag-1 TE k={k} computed in {time.time() - t0:.0f}s')
    f = CACHE / 'te_lag0_common_k1.npy'
    if f.exists():
        TE0k1 = np.load(f)
    else:
        t0 = time.time(); TE0k1 = tk.te_matrix(X, np.ones(G, int), delay=0)[0]; np.save(f, TE0k1)
        log(f'common-history lag-0 CMI k=1 computed in {time.time() - t0:.0f}s')

    # ---------- H1 and H2 (AUCs, lag0 - lag1 contrast)
    log('H1/H2: AUCs and contrasts')
    H12 = {}
    for level, (pw, nw) in (('group', (WIRED, UNW)), ('physical', (POS, NEG))):
        pt, bd = auc_pair(TE0, TE1, pw, nw, off, draws)
        H12[level] = dict(lag1_auc=summarise(dict(v=pt['b']), {s: dict(v=bd[s]['b']) for s in bd}, ['v'])['v'],
                          lag0_auc=summarise(dict(v=pt['a']), {s: dict(v=bd[s]['a']) for s in bd}, ['v'])['v'],
                          contrast=summarise(dict(v=pt['diff']), {s: dict(v=bd[s]['diff']) for s in bd}, ['v'])['v'])
        H12[level]['contrast']['median_target_main'] = float(np.median(bd['target_main']['diff']))
        H12[level]['contrast']['fraction_positive'] = {s: float((bd[s]['diff'] > 0).mean()) for s in bd}
        p1, c1 = boot_p(pt['b'], bd['target_main']['b'], 0.5)
        p2, c2 = boot_p(pt['diff'], bd['target_main']['diff'], 0.0)
        H12[level]['raw_p_lag1_auc_vs_0.5'] = dict(p=p1, c=c1)
        H12[level]['raw_p_contrast_vs_0'] = dict(p=p2, c=c2)
        for sch in [x for x in draws if x != 'target_main']:
            H12[level][f'raw_p_lag1_auc_vs_0.5_{sch}'] = boot_p(pt['b'], bd[sch]['b'], 0.5)[0]
            H12[level][f'raw_p_contrast_vs_0_{sch}'] = boot_p(pt['diff'], bd[sch]['diff'], 0.0)[0]
        log(f"  {level}: lag1 {pt['b']:.4f} {H12[level]['lag1_auc']['target_main']}  contrast {pt['diff']:+.4f} "
            f"tm {H12[level]['contrast']['target_main']} 2w {H12[level]['contrast']['two_way']} vx {H12[level]['contrast']['vertex']}")
    res['H1_H2'] = H12
    check['H1_H2'] = dict(
        csv_group_lag1_auc='0.469 [0.445, 0.493]', csv_physical_lag1_auc='0.500 [0.411, 0.596]',
        csv_group_lag0_auc='0.632 [0.611, 0.651]', csv_physical_lag0_auc='0.593 [0.504, 0.683]',
        csv_group_contrast_ci=[0.14172079227395068, 0.1842956717001737], csv_group_contrast_median=0.16203561573895892,
        csv_physical_contrast_ci=[0.06464994515733234, 0.12806564631667805], csv_physical_contrast_median=0.0938906464064711)

    # ---------- H3 orientation
    log('H3: orientation accuracy')
    uni = (POS > 0) & (POS.T == 0) & off
    I, J = np.nonzero(uni)
    scores = {'k=1': common[1], 'k=2': common[2], 'k=4': common[4], 'k=8': common[8], 'AIS-selected k': TE1}
    H3 = {}
    for block, keep in (('all', np.ones(len(I), bool)), ('clock_excluded', I != clock)):
        Ib, Jb = I[keep], J[keep]
        H3[block] = dict(unidirectional_group_pairs=int(keep.sum()))
        for level in ('group', 'physical'):
            W = np.ones(len(Ib)) if level == 'group' else POS[Ib, Jb].astype(float)
            for name, S in scores.items():
                pt, bd = orientation_block(S, rate, Ib, Jb, W, draws)
                r = summarise(pt, bd, ['acc', 'rule', 'diff'])
                p, c = boot_p(pt['diff'], bd['target_main']['diff'], 0.0)
                r['raw_p_diff_vs_0'] = dict(p=p, c=c)
                r['raw_p_one_sided_diff_le_0'] = min(1.0, (int((bd['target_main']['diff'] <= 0).sum()) + 1) / (B + 1))
                H3[block][f'{level}|{name}'] = r
                log(f"  {block} {level} {name}: acc {pt['acc']:.3f} rule {pt['rule']:.3f} diff {pt['diff']:+.3f} "
                    f"tm {r['diff']['target_main']} 2w {r['diff']['two_way']} vx {r['diff']['vertex']}")
    res['H3'] = H3
    check['H3'] = dict(csv_physical_accuracy={'AIS-selected k': '0.755 [0.710, 0.797]', 'k=1': '0.695 [0.642, 0.744]',
                                              'k=2': '0.665 [0.612, 0.714]', 'k=4': '0.707 [0.658, 0.751]',
                                              'k=8': '0.744 [0.696, 0.786]'},
                       csv_physical_rate_prespecified='0.268 [0.213, 0.321]')
    res['self_check'] = check
    (OUT / ('selftest_pass3.json' if a.selftest else 'recomputed_supplementary_pass3.json')).write_text(json.dumps(res, indent=1, default=float))

    # ---------- H4 signed influence
    log('H4: percentage negative delta')
    cls = {'pull-down only': (PD > 0) & (PASS == 0) & off, 'pass only': (PASS > 0) & (PD == 0) & off,
           'unwired': (POS == 0) & off}
    wts = {'pull-down only': PD, 'pass only': PASS, 'unwired': NEG}
    ph = R.phase_classes(X)
    gi, gj = np.nonzero(off)
    strat = np.full((G, G), -1)
    sp = ph[gi] == ph[gj]; mx = (ph[gi] == 2) | (ph[gj] == 2)
    strat[gi, gj] = np.where(mx, 2, np.where(sp, 0, 1))         # 0 same parity, 1 opposite, 2 involves mixed
    sname = ['same parity', 'opposite parity', 'involves a mixed group']
    items = {}
    for c, m in cls.items():
        ok = m & np.isfinite(DELTA) & (wts[c] > 0)
        ii, jj = np.nonzero(ok)
        items[c] = dict(i=ii, j=jj, neg=(DELTA[ii, jj] < 0).astype(float), wphys=wts[c][ii, jj].astype(float), st=strat[ii, jj])
    H4 = dict(counts={}, percent_negative={})
    for c, d in items.items():
        H4['counts'][c] = {'all': int(len(d['i']))} | {sname[s]: int((d['st'] == s).sum()) for s in range(3)}
        H4['percent_negative'][c] = {'all_group': pct_neg(d['neg'], np.ones(len(d['i']))),
                                     'all_physical': pct_neg(d['neg'], d['wphys'])} | \
                                    {sname[s]: pct_neg(d['neg'][d['st'] == s], np.ones((d['st'] == s).sum())) for s in range(3)}

    def h4_stats(a_=None, b_=None):
        out = {}
        w = {}
        for c in ('pull-down only', 'unwired'):
            d = items[c]; ww = np.ones(len(d['i']))
            if a_ is not None: ww = ww * a_[d['i']]
            if b_ is not None: ww = ww * b_[d['j']]
            w[c] = ww
        pdi, uni_ = items['pull-down only'], items['unwired']
        out['group'] = pct_neg(pdi['neg'], w['pull-down only']) - pct_neg(uni_['neg'], w['unwired'])
        out['physical'] = pct_neg(pdi['neg'], w['pull-down only'] * pdi['wphys']) - pct_neg(uni_['neg'], w['unwired'] * uni_['wphys'])
        num = den = 0.0
        for s in range(3):
            mp, mu = pdi['st'] == s, uni_['st'] == s
            dsd = pct_neg(pdi['neg'][mp], w['pull-down only'][mp]) - pct_neg(uni_['neg'][mu], w['unwired'][mu])
            out[f'stratum:{sname[s]}'] = dsd
            ws = w['pull-down only'][mp].sum()
            if ws > 0 and np.isfinite(dsd):
                num += ws * dsd; den += ws
        out['composition_adjusted'] = num / den if den > 0 else np.nan
        return out
    pt = h4_stats()
    bd = {sch: {} for sch in draws}
    for sch, (A, Bm) in draws.items():
        rows = [h4_stats(A[d], Bm[d]) for d in range(B)]
        for k in pt:
            bd[sch][k] = np.array([r[k] for r in rows], float)
    H4['differences'] = summarise(pt, bd, list(pt))
    for k in ('group', 'physical', 'composition_adjusted'):
        p, c = boot_p(pt[k], bd['target_main'][k], 0.0)
        H4['differences'][k]['raw_p_vs_0'] = dict(p=p, c=c)
    log(f"  H4 group {pt['group']:+.2f} {H4['differences']['group']}")
    log(f"  H4 physical {pt['physical']:+.2f} {H4['differences']['physical']}")
    log(f"  H4 composition-adjusted {pt['composition_adjusted']:+.2f} {H4['differences']['composition_adjusted']}")
    res['H4'] = H4
    if a.selftest:
        (OUT / 'selftest_pass3.json').write_text(json.dumps(res, indent=1, default=float)); log('SELF-TEST DONE'); return
    (OUT / ('selftest_pass3.json' if a.selftest else 'recomputed_supplementary_pass3.json')).write_text(json.dumps(res, indent=1, default=float))

    # ---------- Holm over H1..H6 at the deciding level (target-cluster, main-run generator)
    prev = json.load(open(OUT / 'recomputed_supplementary.json'))['D_rq2_increments']['total_idtxl_minus_global_cut']
    h6 = json.load(open(OUT / 'h6_run_of_record_R1.json'))['differences']['R1_minus_lag1']['group_precision']
    raw = {'H1': H12['physical']['raw_p_lag1_auc_vs_0.5']['p'],
           'H2': H12['group']['raw_p_contrast_vs_0']['p'],
           'H3': max(H3['all'][f'group|{k}']['raw_p_diff_vs_0']['p'] for k in ('k=1', 'k=2', 'k=4', 'k=8')),
           'H4': H4['differences']['group']['raw_p_vs_0']['p'],
           # H5, H6 from the saved draws' fraction positive (c = draws <= 0 = (1 - fraction positive) * B)
           'H5': min(1.0, 2 * (round((1 - prev['target_positive']) * B) + 1) / (B + 1)),
           'H6': min(1.0, 2 * (round((1 - h6['target_positive']) * B) + 1) / (B + 1))}
    order = sorted(raw, key=lambda h: raw[h]); m = len(order); holm = {}; run = 0.0
    for r_, h in enumerate(order):
        run = max(run, min(1.0, (m - r_) * raw[h])); holm[h] = run
    res['holm'] = dict(raw=raw, holm=holm,
                       note='H3 raw p: maximum over k = 1, 2, 4, 8 of the two-sided p for the group difference '
                            '(intersection-union); H5 and H6 from the saved target-cluster fractions of positive draws.')
    log(f'Holm: raw {raw} -> {holm}')
    (OUT / ('selftest_pass3.json' if a.selftest else 'recomputed_supplementary_pass3.json')).write_text(json.dumps(res, indent=1, default=float))

    # ---------- H2 robustness rows of tab:r1-h2
    log('H2 robustness: clock group excluded, common history k = 1')
    rob = {}
    noclk = off.copy(); noclk[clock, :] = False; noclk[:, clock] = False
    for label, SA, SB, mask in (('clock_group_excluded', TE0, TE1, noclk), ('common_history_k1', TE0k1, common[1], off)):
        rob[label] = {}
        for level, (pw, nw) in (('group', (WIRED, UNW)), ('physical', (POS, NEG))):
            pt_, bd_ = auc_pair(SA, SB, pw, nw, mask, {'target_main': draws['target_main']})
            rob[label][level] = dict(lag0_auc=dict(estimate=pt_['a'], target_main=q(bd_['target_main']['a'])),
                                     lag1_auc=dict(estimate=pt_['b'], target_main=q(bd_['target_main']['b'])),
                                     contrast=dict(estimate=pt_['diff'], target_main=q(bd_['target_main']['diff'])))
            log(f"  {label} {level}: lag0 {pt_['a']:.3f} {q(bd_['target_main']['a'])} lag1 {pt_['b']:.3f} "
                f"contrast {pt_['diff']:+.3f} {q(bd_['target_main']['diff'])}")
        # within-target AUC of lag 0 (panel A)
        rob[label]['lag0_within_target_auc'] = R.within_target_auc(SA, POS > 0, mask)
    res['H2_robustness'] = rob
    (OUT / ('selftest_pass3.json' if a.selftest else 'recomputed_supplementary_pass3.json')).write_text(json.dumps(res, indent=1, default=float))

    # ---------- node level (tab:app-nodelevel; tab:r1-h2 node row)
    log('H2 robustness: node level')
    res['H2_node_level'] = node_level(ref, p2g, G, TE0, TE1, MI)
    (OUT / ('selftest_pass3.json' if a.selftest else 'recomputed_supplementary_pass3.json')).write_text(json.dumps(res, indent=1, default=float))
    log('ALL DONE')


def node_level(ref, p2g, G, TE0, TE1, MI):
    """Gate nodes of active transistors; node u -> node v when a transistor gated by u has v as a channel node
    (supplies excluded), both nodes active, different groups. A node's group: the group of the transistors it gates
    (all identical under change coding; nodes whose gated transistors span groups are reported and dropped)."""
    gate, c1, c2 = ref['gate'], ref['c1'], ref['c2']
    act = p2g >= 0
    node_groups = {}
    for t in np.flatnonzero(act):
        node_groups.setdefault(int(gate[t]), set()).add(int(p2g[t]))
    multi = sum(len(v) > 1 for v in node_groups.values())
    ng = {n: next(iter(v)) for n, v in node_groups.items() if len(v) == 1}
    edges = set()
    for t in range(len(gate)):
        u = int(gate[t])
        if u not in ng: continue
        for c in (int(c1[t]), int(c2[t])):
            if c in (R.VSS, R.VCC) or c == u or c not in ng: continue
            if ng[u] != ng[c]: edges.add((u, c))
    nodes = sorted(ng); info = dict(active_gate_nodes=len(nodes), nodes_spanning_groups=multi, node_edges=len(edges),
                                    printed_node_edges=1645, printed_active_gate_nodes=1023)
    log(f"  node level: {info}")
    if len(edges) != 1645 or len(nodes) != 1023:
        info['status'] = 'not recomputed: this node definition does not reproduce the printed counts (1,645 edges, 1,023 nodes)'
        log('  ' + info['status'])
        return info
    idx = {n: i for i, n in enumerate(nodes)}; N = len(nodes); g = np.array([ng[n] for n in nodes])
    E = np.zeros((N, N), bool)
    for u, v in edges: E[idx[u], idx[v]] = True
    mask = g[:, None] != g[None, :]
    rng = np.random.default_rng(SEED)
    tm = np.array([np.bincount(rng.integers(0, N, N), minlength=N) for _ in range(B)], float)
    out = dict(info)
    res_ = {}
    for name, S in (('lag1_TE', TE1), ('lag0_CMI', TE0), ('MI', MI)):
        Sn = S[g[:, None], g[None, :]]
        e = AUC(Sn, E, ~E, mask); res_[name] = (e(), np.array([e(None, tm[d]) for d in range(B)]))
        out[name] = dict(node_auc=res_[name][0], ci=q(res_[name][1]))
    d = res_['lag0_CMI'][1] - res_['lag1_TE'][1]
    out['contrast'] = dict(estimate=res_['lag0_CMI'][0] - res_['lag1_TE'][0], ci=q(d))
    log(f"  node level contrast {out['contrast']}")
    return out


if __name__ == '__main__':
    main()
