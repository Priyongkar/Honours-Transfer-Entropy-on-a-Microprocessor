"""Post hoc: sensitivity of the Donkey Kong clock-aligned scores to the start sample (1,407, used by
recompute_supplementary.py, against 1,408, the clock-alignment rule of run_replication.py).
Writes output/posthoc_align_offset.json."""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import recompute_supplementary as R                                    # noqa: E402
X, K, TE1, TE0, MI, p2g, NET = R.load(); G = X.shape[0]
ref = R.references(p2g, G); Wd, Wphys, Nphys = ref['Wd'], ref['Wphys'], ref['Nphys']
off = ~np.eye(G, dtype=bool); ii, jj = np.nonzero(off); p = Wd[ii, jj].astype(float)
mt = np.random.default_rng(R.SEED).multinomial(G, np.full(G, 1 / G), size=200)
ph_all = R.phase_classes(X)
out = {}
for start in (1407, 1408):
    d = R.parity_scores(X, K, sub=start); o = {}
    Xs = X[:, start:]; ph = R.phase_classes(Xs)
    PC = np.full((G, G), np.nan); opp = (ph[ii] != ph[jj]) & (ph[ii] < 2) & (ph[jj] < 2); PC[ii, jj] = np.where(opp, 0., 1.)
    d['phase_class'] = PC
    for nm, S in d.items():
        s = np.nan_to_num(S[ii, jj], nan=-1)
        c, _ = R.ci(lambda m: R.weighted_auc(s, p * m[jj], (1 - p) * m[jj]), mt)
        o[nm] = dict(group_auc=R.weighted_auc(s, p, 1 - p), group_auc_ci=c, group_ap_lift=R.ap_lift(s, p, 1 - p),
                     physical_auc=R.weighted_auc(s, Wphys[ii, jj], Nphys[ii, jj]),
                     physical_ap_lift=R.ap_lift(s, Wphys[ii, jj], Nphys[ii, jj]),
                     within_target=R.within_target_auc(np.where(np.isnan(S), -1, S), Wd, off))
    out[str(start)] = o
    print(start, {k: round(v['group_auc'], 3) for k, v in o.items()}, flush=True)
(HERE / 'output/posthoc_align_offset.json').write_text(json.dumps(out, indent=1, default=float))
