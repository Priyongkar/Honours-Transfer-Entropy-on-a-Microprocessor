"""Post hoc, after R3.4 failed on both held-out programs (PREDICTIONS_R3.md): descriptive only.
(1) Donkey Kong clock-aligned: group AUC and AP lift of MI|parity and of phase class only, and their
difference, from sample 1,408 (the clock-alignment rule) and from sample 5,000 (the samples behind 0.661).
(2) Held-out programs, clock-aligned: AP lift of MI|parity and of phase class only.
Writes results/10_remedy_held_out_R3/posthoc_r3_phase_class.json."""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'src'))
from te6502.paths import load_analysis, RESULTS, DATA_ROOT  # noqa: E402
R = load_analysis('09_recomputing_the_supplementary_results/recompute_supplementary_results_pass1.py')
R2 = load_analysis('07_replication_on_pitfall_and_space_invaders_R2/replication_pairwise_pipeline_R2.py')
zero_lag = load_analysis('10_phase_aware_remedy_on_held_out_programs_R3/remedy_on_held_out_programs_R3.py').zero_lag


def phase_score(X, ii, jj):
    ph = R.phase_classes(X)
    opp = (ph[ii] != ph[jj]) & (ph[ii] < 2) & (ph[jj] < 2)
    return np.where(opp, 0.0, 1.0)


def block(X, Wd, mt):
    G = X.shape[0]; off = ~np.eye(G, dtype=bool); ii, jj = np.nonzero(off); p = Wd[ii, jj].astype(float)
    mi = np.nan_to_num(zero_lag(X, True)[ii, jj], nan=-1); pc = phase_score(X, ii, jj)
    out = {}
    for nm, s in (('mi_parity', mi), ('phase_class', pc)):
        c, _ = R.ci(lambda m: R.weighted_auc(s, p * m[jj], (1 - p) * m[jj]), mt)
        out[nm] = dict(group_auc=R.weighted_auc(s, p, 1 - p), group_auc_ci=c, ap_lift=R.ap_lift(s, p, 1 - p))
    d = np.array([R.weighted_auc(mi, p * m[jj], (1 - p) * m[jj]) - R.weighted_auc(pc, p * m[jj], (1 - p) * m[jj]) for m in mt])
    out['difference'] = dict(est=out['mi_parity']['group_auc'] - out['phase_class']['group_auc'],
                             ci=[float(np.quantile(d, .025)), float(np.quantile(d, .975))])
    return out


def main():
    out = {}
    X, K, TE1, TE0, MI, p2g, NET = R.load(); G = X.shape[0]
    Wd = R.references(p2g, G)['Wd']
    mt = np.random.default_rng(R.SEED).multinomial(G, np.full(G, 1 / G), size=200)
    for start in (1408, 5000):
        out[f'donkeykong_from_{start}'] = block(X[:, start:], Wd, mt)
        print(start, json.dumps(out[f'donkeykong_from_{start}'], default=float))
    D = DATA_ROOT / 'data/regen'
    for g, f in (('Pitfall', 'pitfall'), ('SpaceInvaders', 'spaceinvaders')):
        gate, c1, c2, centre, pullup, names, raw, phase, red, start, n_off, reps, p2g2, X2, K2 = R2.prepare(D / f'{f}.regen.ts.spikes.npy')
        G2 = len(reps)
        Wd2 = R2.references(p2g2, G2, gate, c1, c2, centre, pullup, X2.mean(1))[0]
        mt2 = np.random.default_rng(R.SEED).multinomial(G2, np.full(G2, 1 / G2), size=200)
        out[g] = block(X2[:, start:], Wd2, mt2); print(g, json.dumps(out[g], default=float))
    (RESULTS / '10_remedy_held_out_R3' / 'posthoc_r3_phase_class.json').write_text(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
