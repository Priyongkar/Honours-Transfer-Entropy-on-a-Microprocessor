"""Second pass for the items the first pass defined differently from the thesis or did not compute:
phase-class baselines (opposite-parity last, others tied); the physical-fitted covariate model (rows per
transistor-pair label); the cross-validated rates model; source self-information and entropy deciles;
Holm correction across all eight classes jointly; zero-lag MI and MI|parity from sample 5,000 onwards.
Writes results/09_recomputed_supplementary/recomputed_supplementary_extras.json."""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'src'))
from te6502.paths import load_analysis, RESULTS  # noqa: E402
R = load_analysis('09_recomputing_the_supplementary_results/recompute_supplementary_results_pass1.py')                                    # noqa: E402
tk = R.tk


def main():
    X, K, TE1, TE0, MI, p2g, NET = R.load(); G, T = X.shape
    ref = R.references(p2g, G); Wd, Wphys, Nphys, size, dist = ref['Wd'], ref['Wphys'], ref['Nphys'], ref['size'], ref['dist']
    off = ~np.eye(G, dtype=bool); ii, jj = np.nonzero(off); p = Wd[ii, jj].astype(float)
    mt = np.random.default_rng(R.SEED).multinomial(G, np.full(G, 1 / G), size=200)
    ph = R.phase_classes(X); rate = X.mean(1)
    out = {}
    # phase class only: opposite parity last, all others tied
    opp = (ph[ii] != ph[jj]) & (ph[ii] < 2) & (ph[jj] < 2)
    cls = np.where(opp, 0.0, 1.0)
    PC = np.full((G, G), np.nan); PC[ii, jj] = cls
    r = dict(group_auc=R.weighted_auc(cls, p, 1 - p), group_ap_lift=R.ap_lift(cls, p, 1 - p),
             physical_auc=R.weighted_auc(cls, Wphys[ii, jj], Nphys[ii, jj]), within_target=R.within_target_auc(PC, Wd, off))
    r['group_auc_ci'], _ = R.ci(lambda m: R.weighted_auc(cls, p * m[jj], (1 - p) * m[jj]), mt)
    out['phase_class_only'] = r
    rk = np.argsort(np.argsort(np.nan_to_num(TE0[ii, jj], nan=-1))) / len(ii)
    combo = cls * 2 + rk; CB = np.full((G, G), np.nan); CB[ii, jj] = combo
    r = dict(group_auc=R.weighted_auc(combo, p, 1 - p), group_ap_lift=R.ap_lift(combo, p, 1 - p),
             physical_auc=R.weighted_auc(combo, Wphys[ii, jj], Nphys[ii, jj]), within_target=R.within_target_auc(CB, Wd, off))
    out['phase_class_then_lag0'] = r
    # covariate model, physical-fitted: one row per label with its transistor-pair weight
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold
    pc = np.zeros((len(ii), 3)); pc[np.arange(len(ii)), np.where(ph[ii] == ph[jj], 0, np.where((ph[ii] == 2) | (ph[jj] == 2), 2, 1))] = 1
    F0 = np.column_stack([np.log(size[ii]), np.log(size[jj]), np.log1p(dist[ii, jj]),
                          np.log(np.maximum(rate[ii], 1e-6)), np.log(np.maximum(rate[jj], 1e-6)), pc])
    F1 = np.column_stack([F0, np.log1p(1e4 * np.nan_to_num(TE0[ii, jj]))])
    wp, wn = Wphys[ii, jj], Nphys[ii, jj]
    res = {}
    for name, F in (('covariates', F0), ('covariates_plus_lag0', F1)):
        pr = np.zeros(len(ii))
        for tr, te_ in GroupKFold(5).split(F, p, groups=jj):
            Xf = np.vstack([F[tr], F[tr]]); yf = np.r_[np.ones(len(tr)), np.zeros(len(tr))]; wf = np.r_[wp[tr], wn[tr]]
            keep = wf > 0; mu, sd = Xf[keep].mean(0), Xf[keep].std(0) + 1e-12
            m = LogisticRegression(max_iter=3000).fit((Xf[keep] - mu) / sd, yf[keep], sample_weight=wf[keep] / wf[keep].mean())
            pr[te_] = m.decision_function((F[te_] - mu) / sd)
        res[name] = pr
    one = np.ones(G)
    met = lambda s, m: R.weighted_auc(s, wp * m[jj], wn * m[jj]); mtap = lambda s, m: R.ap_lift(s, wp * m[jj], wn * m[jj])
    cr = {}
    for mn, f in (('physical_auc', met), ('physical_ap_lift', mtap)):
        a0, a1 = f(res['covariates'], one), f(res['covariates_plus_lag0'], one)
        b = np.array([f(res['covariates_plus_lag0'], m) - f(res['covariates'], m) for m in mt])
        c0, _ = R.ci(lambda m: f(res['covariates'], m), mt)
        cr[mn] = dict(covariates=a0, covariates_ci=c0, increment=a1 - a0, increment_ci=[float(np.quantile(b, .025)), float(np.quantile(b, .975))])
    out['covariate_model_physical_fitted'] = cr
    # rates model: cross-validated logistic of the two log rates, allocation-matched precision
    Fr = np.column_stack([np.log(np.maximum(rate[ii], 1e-6)), np.log(np.maximum(rate[jj], 1e-6))])
    pr = np.zeros(len(ii))
    for tr, te_ in GroupKFold(5).split(Fr, p, groups=jj):
        m = LogisticRegression(max_iter=2000).fit(Fr[tr], p[tr]); pr[te_] = m.decision_function(Fr[te_])
    RM = np.full((G, G), np.nan); RM[ii, jj] = pr
    alloc = NET.sum(0)
    prec, pprec = R.alloc_precision(RM, Wd, alloc, wphys=Wphys, nphys=Nphys)
    c, _ = R.ci(lambda m: R.alloc_precision(RM, Wd, alloc, w_t=m.astype(float))[0], mt)
    out['rates_model'] = dict(group_auc=R.weighted_auc(pr, p, 1 - p), group_ap_lift=R.ap_lift(pr, p, 1 - p),
                              within_target=R.within_target_auc(RM, Wd, off), precision_at_allocation=prec,
                              precision_ci=c, physical_precision=pprec)
    # source self-information and entropy deciles
    def selfinfo(x):
        a, b = x[:-1].astype(int), x[1:].astype(int); j = a * 2 + b
        c = np.bincount(j, minlength=4).astype(float) / len(j); pa = c.reshape(2, 2).sum(1); pb = c.reshape(2, 2).sum(0)
        return float(sum(c[k] * np.log2(c[k] / (pa[k // 2] * pb[k % 2])) for k in range(4) if c[k] > 0))
    def ent(x):
        q = x.mean(); return float(-(q * np.log2(q) + (1 - q) * np.log2(1 - q))) if 0 < q < 1 else 0.0
    props = {'source_self_information': np.array([selfinfo(x) for x in X]), 'source_entropy': np.array([ent(x) for x in X])}
    for nm, v in props.items():
        st = R.deciles(v[ii]); r = {}
        for sc, S in (('lag1', TE1), ('lag0', TE0)):
            s = S[ii, jj]; r[f'{sc}_group'] = R.strat_auc(s, p, 1 - p, st); r[f'{sc}_physical'] = R.strat_auc(s, Wphys[ii, jj], Nphys[ii, jj], st)
        out[f'within_deciles_{nm}'] = r
    # Holm across all eight non-unrelated classes jointly
    base = json.load(open(RESULTS / '09_recomputed_supplementary/recomputed_supplementary.json'))['E_enrichment']
    ps = sorted((base[d][c]['exact_p'], d, c) for d in base for c in ('reverse_only', 'two_step', 'shared_parent', 'shared_child'))
    holm, run = {}, 0.0
    for i, (pv, d, c) in enumerate(ps):
        run = max(run, min(1.0, (len(ps) - i) * pv)); holm[f'{d}:{c}'] = run
    out['holm_joint'] = holm
    # zero-lag MI and MI|parity from sample 5,000 onwards
    Xs = X[:, 5000:]; XS = tk.as_sparse_sources(Xs); par = (np.arange(Xs.shape[1]) % 2).astype(np.uint8)[None, :]
    MIa = np.full((G, G), np.nan); MIp = np.full((G, G), np.nan)
    for t in range(G):
        MIa[:, t] = tk.te_all_sources(XS, Xs[t], 0, delay=0)[0]
        MIp[:, t] = tk.te_all_sources(XS, Xs[t], 0, delay=0, cond=par, cond_delay=0)[0]
    for M in (MIa, MIp): np.fill_diagonal(M, np.nan)
    for nm, S in (('mi_from_5000', MIa), ('mi_parity_from_5000', MIp)):
        s = S[ii, jj]
        out[nm] = dict(group_auc=R.weighted_auc(s, p, 1 - p), physical_auc=R.weighted_auc(s, Wphys[ii, jj], Nphys[ii, jj]),
                       within_target=R.within_target_auc(S, Wd, off), group_ap_lift=R.ap_lift(s, p, 1 - p))
    (RESULTS / '09_recomputed_supplementary/recomputed_supplementary_extras.json').write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
