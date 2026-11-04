"""H1 robustness checks quoted in the text: each target's strongest source wired (lag-1 TE vs lag-0 CMI);
group AUC of lag-1 TE after subtracting the analytic bias; group AUC ranked by exact analytic p-value.
Supplementary. Writes output/h1_robustness_checks.json."""
import json, sys
from pathlib import Path
import numpy as np
from scipy.stats import chi2
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'src'))
from te6502.paths import load_analysis, RESULTS, DATA_ROOT  # noqa: E402
R = load_analysis('09_recomputing_the_supplementary_results/recompute_supplementary_results_pass1.py')                                    # noqa: E402


def main():
    X, K, TE1, TE0, MI, p2g, NET = R.load(); G, T = X.shape
    ref = R.references(p2g, G); Wd = ref['Wd']
    off = ~np.eye(G, dtype=bool); ii, jj = np.nonzero(off); p = Wd[ii, jj].astype(float)
    has = Wd.sum(0) > 0
    def top1(S):
        hit = [Wd[np.nanargmax(np.where(np.isfinite(S[:, t]), S[:, t], -np.inf)), t] for t in np.flatnonzero(has)]
        return 100 * float(np.mean(hit))
    k = K[jj]; n = T - k - 1
    df = 2.0 ** k                                                        # binary source, l = 1
    bias = df / (2 * n * np.log(2))
    s = TE1[ii, jj]; s_bc = s - bias
    lp = chi2.logsf(2 * n * np.log(2) * np.nan_to_num(s), df)
    out = dict(targets_with_wired_source=int(has.sum()), top1_wired_lag1_te=top1(TE1), top1_wired_lag0_cmi=top1(TE0),
               group_auc_lag1=R.weighted_auc(s, p, 1 - p), group_auc_lag1_bias_subtracted=R.weighted_auc(s_bc, p, 1 - p),
               group_auc_lag1_by_exact_p=R.weighted_auc(-lp, p, 1 - p))
    (RESULTS / '09_recomputed_supplementary/h1_robustness_checks.json').write_text(json.dumps(out, indent=1)); print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
