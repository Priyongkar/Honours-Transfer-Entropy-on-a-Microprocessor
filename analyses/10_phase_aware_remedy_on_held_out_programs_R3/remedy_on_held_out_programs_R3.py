"""R3: the phase-aware remedy (zero-lag MI given the phase) on the held-out programs (design: PREDICTIONS_R3.md).

Usage: python run_remedy_heldout.py GAME RECORDING.npy      GAME in {Pitfall, SpaceInvaders}
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / 'src'))  # repository library
from te6502.paths import ANALYSES, DATA_ROOT, RESULTS, HELPERS_FILE, load_analysis  # noqa: E402,F401
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
R2 = load_analysis('07_replication_on_pitfall_and_space_invaders_R2/replication_pairwise_pipeline_R2.py')
tk, log, phase_class, boot_mult, auc_with_ci, diff_ci = R2.tk, R2.log, R2.phase_class, R2.boot_mult, R2.auc_with_ci, R2.diff_ci

PRED = HERE / 'design/predictions_R3.md'
OUT = RESULTS / '10_remedy_held_out_R3'; OUT.mkdir(exist_ok=True)


def zero_lag(X, cond):
    """I(X_t ; Y_t), or I(X_t ; Y_t | parity_t), for every source and target, no target history."""
    G, T = X.shape; XS = tk.as_sparse_sources(X); M = np.full((G, G), np.nan)
    par = (np.arange(T) % 2).astype(np.uint8)[None, :]
    for t in range(G):
        if cond:
            M[:, t] = tk.te_all_sources(XS, X[t], 0, delay=0, cond=par, cond_delay=0)[0]
        else:
            M[:, t] = tk.te_all_sources(XS, X[t], 0, delay=0)[0]
    np.fill_diagonal(M, np.nan)
    return M


def analyse(X, K, Wd, Wphys, Nphys, dist, size):
    G, T = X.shape; off = ~np.eye(G, dtype=bool); ii, jj = np.nonzero(off)
    sc = R2.scores(X, K); TE1, TE0 = sc[1], sc[0]
    MIp, MI = zero_lag(X, True), zero_lag(X, False)
    ph, _ = phase_class(X)
    opp = (ph[ii] != ph[jj]) & (ph[ii] < 2) & (ph[jj] < 2)
    PC = np.full((G, G), np.nan); PC[ii, jj] = np.where(opp, 0.0, 1.0)
    SZ = np.repeat(size[:, None], G, 1).astype(float)
    mult = boot_mult(G, R2.DRAWS, R2.SEED); one = np.zeros(len(ii), int); pos = Wd[ii, jj]
    res, bs = {}, {}
    for name, S in (('mi_parity', MIp), ('mi', MI), ('phase_class', PC), ('source_size', SZ),
                    ('lag1_te', TE1), ('lag0_cmi', TE0), ('distance', -dist)):
        res[f'group_auc_{name}'], bs[name] = auc_with_ci(np.nan_to_num(S[ii, jj], nan=-np.inf), pos, one, jj, mult)
    res['R3.1_group_auc_mi_parity'] = res['group_auc_mi_parity']
    a = lambda n: res[f'group_auc_{n}']['auc']
    res['R3.2_mi_parity_minus_lag1'] = diff_ci(bs['mi_parity'], bs['lag1_te'], a('mi_parity'), a('lag1_te'))
    res['R3.3_mi_parity_minus_lag0'] = diff_ci(bs['mi_parity'], bs['lag0_cmi'], a('mi_parity'), a('lag0_cmi'))
    res['R3.4_mi_parity_minus_phase_class'] = diff_ci(bs['mi_parity'], bs['phase_class'], a('mi_parity'), a('phase_class'))
    res['R3.5_distance_minus_mi_parity'] = diff_ci(bs['distance'], bs['mi_parity'], a('distance'), a('mi_parity'))
    res['physical_auc_mi_parity'] = R2.phys_auc_ci(MIp, Wphys, Nphys, mult)
    res['within_target_auc_mi_parity'], _ = auc_with_ci(np.nan_to_num(MIp[ii, jj], nan=-np.inf), pos, jj, jj, mult)
    return res


def met(r):
    return {'R3.1': bool(r['R3.1_group_auc_mi_parity']['lo'] > 0.5),
            'R3.2': bool(r['R3.2_mi_parity_minus_lag1']['lo'] > 0),
            'R3.3': bool(r['R3.3_mi_parity_minus_lag0']['lo'] > 0),
            'R3.4': bool(r['R3.4_mi_parity_minus_phase_class']['lo'] > 0),
            'R3.5': bool(r['R3.5_distance_minus_mi_parity']['lo'] > 0)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('game'); ap.add_argument('recording'); a = ap.parse_args()
    hashes = {p.name: R2.sha(p) for p in (PRED, Path(__file__).resolve(), HELPERS_FILE,
                                          ANALYSES / '07_replication_on_pitfall_and_space_invaders_R2/replication_pairwise_pipeline_R2.py')}
    for k, v in hashes.items():
        log(f'R3-{a.game}: {k} sha256 {v}')
    gate, c1, c2, centre, pullup, names, raw, phase, red, start, n_off, reps, p2g, X, K = R2.prepare(a.recording)
    G = len(reps); rate = X.mean(1)
    Wd, Wphys, Nphys, Lat, dist, size, clock_g = R2.references(p2g, G, gate, c1, c2, centre, pullup, rate)
    res = dict(game=a.game, hashes=hashes, samples=int(red.shape[1]), clock_aligned_start=start, groups=int(G),
               wired_group_pairs=int(Wd.sum()))
    log(f'{a.game}: {res}')
    res['clock_aligned'] = analyse(X[:, start:], K, Wd, Wphys, Nphys, dist, size)
    log(f'{a.game}: clock-aligned done')
    res['whole_recording'] = analyse(X, K, Wd, Wphys, Nphys, dist, size)
    res['predictions_met'] = met(res['clock_aligned'])
    def conv(v):
        if isinstance(v, dict): return {str(k): conv(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)): return [conv(x) for x in v]
        if isinstance(v, np.integer): return int(v)
        if isinstance(v, np.floating): return float(v)
        if isinstance(v, np.bool_): return bool(v)
        return v
    (OUT / f'R3_{a.game}_results.json').write_text(json.dumps(conv(res), indent=1))
    log(f"{a.game}: done {res['predictions_met']}")


if __name__ == '__main__':
    main()
