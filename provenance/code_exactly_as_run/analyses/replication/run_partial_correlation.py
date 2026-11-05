"""R6: zero-lag partial correlation on the clock-aligned Donkey Kong recording (design: PREDICTIONS_R6.md)."""
import hashlib, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SR = HERE.parent / 'supplementary_recompute'
sys.path.insert(0, str(SR)); sys.path.insert(0, str(HERE))
import recompute_supplementary as R                                    # noqa: E402
from run_remedy_heldout import zero_lag                                # noqa: E402

PRED = HERE / 'PREDICTIONS_R6.md'; START = 1408


def log(m):
    import time; print(f'[{time.strftime("%H:%M:%S")}] {m}', flush=True)


def partial_corr(Z):
    C = np.cov(Z); Th = np.linalg.pinv(C, rcond=1e-10)
    d = np.sqrt(np.clip(np.diag(Th), 1e-300, None))
    P = -Th / np.outer(d, d); np.fill_diagonal(P, np.nan)
    return np.abs(P)


def main():
    for p in (PRED, Path(__file__).resolve()):
        log(f'R6: {p.name} sha256 {hashlib.sha256(p.read_bytes()).hexdigest()}')
    X, K, TE1, TE0, MI, p2g, NET = R.load(); G = X.shape[0]
    ref = R.references(p2g, G); Wd, Wphys, Nphys, dist = ref['Wd'], ref['Wphys'], ref['Nphys'], ref['dist']
    off = ~np.eye(G, dtype=bool); ii, jj = np.nonzero(off); p = Wd[ii, jj].astype(float)
    Xa = X[:, START:].astype(float); T = Xa.shape[1]; par = (np.arange(T) % 2)
    Xr = Xa.copy()
    for q in (0, 1):
        Xr[:, par == q] -= Xa[:, par == q].mean(1, keepdims=True)
    S = {'partial': partial_corr(Xa), 'partial_parity': partial_corr(Xr),
         'correlation': np.abs(np.corrcoef(Xa)) * off}
    np.fill_diagonal(S['correlation'], np.nan)
    log('partial correlations done; zero-lag MI next')
    Xu = X[:, START:]
    S['mi'] = zero_lag(Xu, False); S['mi_parity'] = zero_lag(Xu, True); S['distance'] = -dist
    mt = np.random.default_rng(R.SEED).multinomial(G, np.full(G, 1 / G), size=R.DRAWS)
    out, draws = {}, {}
    for nm, M in S.items():
        s = np.nan_to_num(M[ii, jj], nan=-np.inf)
        r = dict(group_auc=R.weighted_auc(s, p, 1 - p), group_ap_lift=R.ap_lift(s, p, 1 - p),
                 physical_auc=R.weighted_auc(s, Wphys[ii, jj], Nphys[ii, jj]),
                 within_target=R.within_target_auc(np.where(np.isnan(M), -np.inf, M), Wd, off))
        b = np.array([R.weighted_auc(s, p * m[jj], (1 - p) * m[jj]) for m in mt])
        r['group_auc_ci'] = [float(np.quantile(b, .025)), float(np.quantile(b, .975))]
        out[nm] = r; draws[nm] = b; log(f"{nm}: group AUC {r['group_auc']:.3f} {r['group_auc_ci']}")
    def diff(a, b):
        d = draws[a] - draws[b]
        return dict(diff=out[a]['group_auc'] - out[b]['group_auc'], lo=float(np.quantile(d, .025)), hi=float(np.quantile(d, .975)))
    out['R6.1_partial_minus_mi'] = diff('partial', 'mi')
    out['R6.2_partial_parity_minus_partial'] = diff('partial_parity', 'partial')
    out['R6.3_distance_minus_partial_parity'] = diff('distance', 'partial_parity')
    out['R6.4_partial_parity_minus_mi_parity'] = diff('partial_parity', 'mi_parity')
    # Added after the first run (descriptive): the paired difference of partial correlation from plain correlation,
    # the physical AP lift, and group precision at the multivariate network's allocation (as Table 5.2).
    out['R6.5_descriptive_partial_minus_correlation'] = diff('partial', 'correlation')
    alloc = NET.sum(0)
    for nm in ('partial', 'partial_parity', 'correlation', 'mi', 'mi_parity'):
        s = np.nan_to_num(S[nm][ii, jj], nan=-np.inf)
        out[nm]['physical_ap_lift'] = R.ap_lift(s, Wphys[ii, jj], Nphys[ii, jj])
        Msc = np.where(np.isnan(S[nm]), -np.inf, S[nm])
        prec, pprec = R.alloc_precision(Msc, Wd, alloc, wphys=Wphys, nphys=Nphys)
        c, _ = R.ci(lambda m: R.alloc_precision(Msc, Wd, alloc, w_t=m.astype(float))[0], mt)
        out[nm]['precision_at_allocation'] = prec; out[nm]['precision_at_allocation_ci'] = c; out[nm]['physical_precision_at_allocation'] = pprec
        log(f"{nm}: precision at allocation {prec:.2f} {c}")
    out['predictions_met'] = {'R6.1': out['R6.1_partial_minus_mi']['lo'] > 0, 'R6.2': out['R6.2_partial_parity_minus_partial']['lo'] > 0,
                              'R6.3': out['R6.3_distance_minus_partial_parity']['lo'] > 0}
    (HERE / 'output/R6_partial_correlation_results.json').write_text(json.dumps(out, indent=1, default=float))
    log(f"done {out['predictions_met']}")


if __name__ == '__main__':
    main()
