"""R7 (design: PREDICTIONS_R7.md): locking gradient within one wiring class; matched-bound inversion matched on k."""
import hashlib, json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; NT = HERE.parent / 'new_tests'
sys.path.insert(0, str(NT))
from common import OUT, RUN, log, boot_mult, auc_with_ci, diff_ci, DRAWS, SEED   # noqa: E402
from a12_recording import lock_index                                             # noqa: E402
from a9_decompose import deciles                                                 # noqa: E402

PRED = HERE / 'PREDICTIONS_R7.md'


def main():
    for p in (PRED, Path(__file__).resolve()):
        log(f'R7: {p.name} sha256 {hashlib.sha256(p.read_bytes()).hexdigest()}')
    ref = np.load(OUT / 'recording_reference.npz'); Wd, PD, K, rate = ref['wired'], ref['pulldown'], ref['k'], ref['rate']
    X = np.load(RUN / 'thesis_v5/idtxl/data_int8.npy').astype(np.uint8)
    TE1 = np.load(RUN / 'te.npz')['te_bits']; G = len(K); off = ~np.eye(G, dtype=bool)
    mult = boot_mult(G, DRAWS, SEED)
    clock_g = int(np.argmax(rate))
    L = lock_index(X); cls = np.where(L < 0.95, 0, np.where(L < 1, 1, 2))
    keep = off.copy(); keep[clock_g, :] = False; keep[:, clock_g] = False
    unw = ~(Wd | Wd.T)
    out = {'classes': {}}
    for cname, Wc in (('pulldown', Wd & PD), ('pass_only', Wd & ~PD)):
        ii, jj = np.nonzero(keep & (Wc | unw)); pos = Wc[ii, jj]; s = np.nan_to_num(TE1[ii, jj], nan=-np.inf)
        o = {'wired_pairs': int(Wc[keep].sum())}; bs = {}
        for c, cn in enumerate(['mixed', 'locked', 'strict']):
            m = cls[jj] == c
            r, b = auc_with_ci(s[m], pos[m], jj[m], jj[m], mult); r['targets_with_wired'] = int(len(np.unique(jj[m & pos])))
            o[cn] = r; bs[cn] = b
        o['mixed_minus_strict'] = diff_ci(bs['mixed'], bs['strict'], o['mixed']['auc'], o['strict']['auc'])
        out['classes'][cname] = o
        log(f"{cname}: mixed {o['mixed']['auc']:.3f} locked {o['locked']['auc']:.3f} strict {o['strict']['auc']:.3f}; diff {o['mixed_minus_strict']['diff']:+.3f} [{o['mixed_minus_strict']['lo']:+.3f}, {o['mixed_minus_strict']['hi']:+.3f}]")
    z = np.load(OUT / 'A3_bound_arrays.npz'); TE, B, W3, ph, strict = z['TE'], z['B'], z['wired'], z['phase'], z['strict']
    S = (ph[:, None] == ph[None, :]) & (ph[:, None] < 2) & strict[:, None] & strict[None, :] & off
    si, sj = np.nonzero(S); te = np.nan_to_num(TE[si, sj], nan=-np.inf); pos = W3[si, sj]; bd = deciles(B[si, sj])
    rB, bB = auc_with_ci(te, pos, bd, sj, mult)
    rK, bK = auc_with_ci(te, pos, bd * 20 + K[sj], sj, mult)
    d = diff_ci(bK, bB, rK['auc'], rB['auc'])
    out['matched_bound'] = dict(pairs=int(len(si)), wired=int(pos.sum()), bound_deciles=rB, bound_deciles_x_k=rK, k_matched_minus_bound_matched=d)
    log(f"matched bound: deciles {rB['auc']:.3f} [{rB['lo']:.3f}, {rB['hi']:.3f}]; x k {rK['auc']:.3f} [{rK['lo']:.3f}, {rK['hi']:.3f}]; diff {d['diff']:+.3f} [{d['lo']:+.3f}, {d['hi']:+.3f}]")
    # Added after the first run (descriptive): the same two stratifications on A3's own scope, all same-phase
    # pairs of locked units (770 wired), which the design file described as strict by mistake.
    S2 = (ph[:, None] == ph[None, :]) & (ph[:, None] < 2) & off
    si, sj = np.nonzero(S2); te = np.nan_to_num(TE[si, sj], nan=-np.inf); pos = W3[si, sj]; bd = deciles(B[si, sj])
    rB2, bB2 = auc_with_ci(te, pos, bd, sj, mult); rK2, bK2 = auc_with_ci(te, pos, bd * 20 + K[sj], sj, mult)
    d2 = diff_ci(bK2, bB2, rK2['auc'], rB2['auc'])
    out['matched_bound_same_phase_scope'] = dict(pairs=int(len(si)), wired=int(pos.sum()), bound_deciles=rB2, bound_deciles_x_k=rK2, k_matched_minus_bound_matched=d2)
    log(f"A3 scope (770 wired): deciles {rB2['auc']:.3f} [{rB2['lo']:.3f}, {rB2['hi']:.3f}]; x k {rK2['auc']:.3f} [{rK2['lo']:.3f}, {rK2['hi']:.3f}]; diff {d2['diff']:+.3f} [{d2['lo']:+.3f}, {d2['hi']:+.3f}]")
    out['predictions_met'] = {'R7.1': bool(out['classes']['pulldown']['mixed_minus_strict']['lo'] > 0),
                              'R7.2': bool(out['classes']['pass_only']['mixed_minus_strict']['lo'] > 0),
                              'R7.3': bool(rK['hi'] < 0.5)}
    def conv(v):
        if isinstance(v, dict): return {str(k): conv(x) for k, x in v.items()}
        if isinstance(v, (np.integer,)): return int(v)
        if isinstance(v, (np.floating,)): return float(v)
        if isinstance(v, np.bool_): return bool(v)
        return v
    (HERE / 'output/R7_locking_within_class_results.json').write_text(json.dumps(conv(out), indent=1))
    log(f"done {out['predictions_met']}")


if __name__ == '__main__':
    main()
