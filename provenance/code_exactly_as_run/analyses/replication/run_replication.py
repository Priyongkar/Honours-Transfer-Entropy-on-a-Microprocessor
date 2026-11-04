"""R2: the pairwise pipeline on a regenerated recording (design: PREDICTIONS_R2.md).

Usage:
  python run_replication.py --selftest            Donkey Kong only: the pipeline must reproduce the main
                                                  run's groups, histories and lag-1 TE; no new statistic
  python run_replication.py GAME RECORDING.npy    GAME in {Pitfall, SpaceInvaders}
"""
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'analyses/new_tests'))
from common import (tk, log, phase_class, boot_mult, auc_with_ci, verdict, diff_ci, quintile,   # noqa: E402
                    bound_all_sources, RUN)
from a9_decompose import decompose, run_tests                                                   # noqa: E402
from a12_recording import lock_index                                                            # noqa: E402
from sim6502 import load_netlist                                                                # noqa: E402

PRED = HERE / 'PREDICTIONS_R2.md'
OUT = HERE / 'output'; OUT.mkdir(exist_ok=True)
VSS, VCC = 558, 657
SEED, DRAWS = 2026, 1000


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def netlist():
    import ast
    defs = []
    for line in (ROOT / 'data/visual6502/transdefs.js').read_text().splitlines():
        if line.strip().startswith("['t"):
            defs.append(ast.literal_eval(line.strip().rstrip(',')))
    defs.sort(key=lambda d: int(d[0][1:]))
    gate, c1, c2 = (np.asarray([d[i] for d in defs]) for i in (1, 2, 3))
    box = np.asarray([d[4] for d in defs], float)
    centre = np.column_stack(((box[:, 0] + box[:, 1]) / 2, (box[:, 2] + box[:, 3]) / 2))
    pullup, _, names = load_netlist(str(ROOT / 'data/visual6502'))
    return gate, c1, c2, centre, pullup, names


def reduce(raw):
    """Main-run reduction: crop to the last nonzero row; keep every third row at the offset with most events."""
    last = int(np.flatnonzero(raw.any(1))[-1]); raw = raw[:last + 1]
    phase = int(np.argmax([raw[o::3].sum() for o in range(3)]))
    return raw, phase, np.ascontiguousarray(raw[phase::3].T)


def clock_start(raw, phase, gate, names):
    """First reduced sample after the last raw clock event off the kept offset (the Donkey Kong rule)."""
    ct = np.flatnonzero(np.isin(gate, [names['cp1'], names['cclk']]))
    rows = np.flatnonzero(raw[:, ct].any(1))
    off = rows[(rows % 3) != phase]
    return 0 if len(off) == 0 else int(off[-1] // 3 + 1), int(len(off))


def references(p2g, G, gate, c1, c2, centre, pullup, rate):
    gated_by = {}
    for t, g in enumerate(gate):
        gated_by.setdefault(int(g), []).append(t)
    pd_nodes = {int(c1[t]) for t in range(len(gate)) if c2[t] == VSS} | {int(c2[t]) for t in range(len(gate)) if c1[t] == VSS}
    driven = lambda nd: bool(pullup.get(nd, False)) or nd in pd_nodes
    Wd = np.zeros((G, G), bool); Wphys = np.zeros((G, G))
    for s in range(len(gate)):
        if p2g[s] < 0:
            continue
        for term in (c1[s], c2[s]):
            if term in (VSS, VCC):
                continue
            for t in gated_by.get(int(term), []):
                if p2g[t] >= 0 and p2g[t] != p2g[s]:
                    Wd[p2g[s], p2g[t]] = True
    # physical counts: wired transistor pairs per group pair (unique pairs)
    pairs = set()
    for s in range(len(gate)):
        if p2g[s] < 0: continue
        for term in (c1[s], c2[s]):
            if term in (VSS, VCC): continue
            for t in gated_by.get(int(term), []):
                if p2g[t] >= 0 and p2g[t] != p2g[s]:
                    pairs.add((s, t))
    for s, t in pairs:
        Wphys[p2g[s], p2g[t]] += 1
    size = np.bincount(p2g[p2g >= 0], minlength=G).astype(float)
    allp = np.outer(size, size); np.fill_diagonal(allp, 0)
    clock_g = int(np.argmax(rate))
    clock_gates = {int(gate[t]) for t in np.flatnonzero(p2g == clock_g)}
    groups_of = lambda nd: {int(p2g[t]) for t in gated_by.get(int(nd), []) if p2g[t] >= 0}
    Lat = np.zeros((G, G), bool)
    for t in range(len(gate)):
        if int(gate[t]) in clock_gates and c1[t] not in (VSS, VCC) and c2[t] not in (VSS, VCC):
            a1, a2 = int(c1[t]), int(c2[t]); d1, d2 = driven(a1), driven(a2)
            if d1 == d2: continue
            for ga in groups_of(a1):
                for gb in groups_of(a2):
                    if ga != gb:
                        if d1: Lat[ga, gb] = True
                        else: Lat[gb, ga] = True
    act = np.flatnonzero(p2g >= 0); order = act[np.argsort(p2g[act])]
    bounds = np.r_[0, np.cumsum(np.bincount(p2g[order], minlength=G))]
    from scipy.spatial.distance import cdist
    dist = np.zeros((G, G))
    for g in range(G):
        dmin = cdist(centre[order[bounds[g]:bounds[g + 1]]], centre[order]).min(0)
        dist[g] = np.minimum.reduceat(dmin, bounds[:-1])
    return Wd, Wphys, allp - Wphys, Lat, dist, size, clock_g


def scores(X, K, which=(1, 0)):
    G = X.shape[0]; XS = tk.as_sparse_sources(X); out = {}
    for d in which:
        M = np.full((G, G), np.nan)
        for t in range(G):
            M[:, t] = tk.te_all_sources(XS, X[t], K[t], delay=d)[0]
        np.fill_diagonal(M, np.nan); out[d] = M
    return out


def phys_auc_ci(S, Wphys, Nphys, mult):
    G = len(S); off = ~np.eye(G, dtype=bool)
    s = np.nan_to_num(S[off], nan=-np.inf); jj = np.nonzero(off)[1]
    p, n = Wphys[off], Nphys[off]
    def wauc(pw, nw):
        o = np.argsort(s, kind='mergesort'); ss, pp, nn = s[o], pw[o], nw[o]
        _, idx = np.unique(ss, return_index=True); bnd = np.r_[idx, len(ss)]
        P = np.add.reduceat(pp, idx); N = np.add.reduceat(nn, idx)
        cumN = np.cumsum(N) - N
        return float((P * (cumN + N / 2)).sum() / (P.sum() * N.sum()))
    point = wauc(p, n); bs = np.array([wauc(p * m[jj], n * m[jj]) for m in mult])
    lo, hi = np.quantile(bs, [.025, .975])
    return dict(auc=point, lo=float(lo), hi=float(hi))


def analyse(X, K, Wd, Wphys, Nphys, Lat, dist, size, clock_g, label):
    G, T = X.shape; off = ~np.eye(G, dtype=bool); ii, jj = np.nonzero(off)
    sc = scores(X, K); TE1, TE0 = sc[1], sc[0]
    mult = boot_mult(G, DRAWS, SEED); one = np.zeros(len(ii), int)
    res = {}
    r1, b1 = auc_with_ci(np.nan_to_num(TE1[ii, jj], nan=-np.inf), Wd[ii, jj], one, jj, mult)
    r0, b0 = auc_with_ci(np.nan_to_num(TE0[ii, jj], nan=-np.inf), Wd[ii, jj], one, jj, mult)
    rd, bd = auc_with_ci(-dist[ii, jj], Wd[ii, jj], one, jj, mult)
    res['R2.1_group_auc_lag1'] = r1; res['group_auc_lag0'] = r0; res['group_auc_distance'] = rd
    res['R2.1b_physical_auc_lag1'] = phys_auc_ci(TE1, Wphys, Nphys, mult)
    res['R2.2_lag0_minus_lag1'] = diff_ci(b0, b1, r0['auc'], r1['auc'])
    res['R2.7_distance_minus_lag0'] = diff_ci(bd, b0, rd['auc'], r0['auc'])
    if label != 'clock_aligned':
        return res
    rate = X.mean(1); ph, strict = phase_class(X); rq = quintile(rate)
    res['R2.3'] = dict(groups=int(G), locked_95=int((ph < 2).sum()), strict=int(strict.sum()),
                       share_locked=float((ph < 2).mean()))
    # R2.4 (A9): strictly locked same-phase pairs
    par = (np.arange(T) % 2).astype(np.int64); Xsp = tk.as_sparse_sources(X)
    locked = ph < 2
    same = (ph[:, None] == ph[None, :]) & locked[:, None] & off
    S = same & strict[:, None] & strict[None, :]
    D, tg = decompose(X, Xsp, K, ph, np.flatnonzero(locked), par)
    B = np.full((G, G), np.nan)
    for j in np.flatnonzero(locked):
        B[:, j] = bound_all_sources(Xsp, X[j], K[j], par)
    si, sj = np.nonzero(S)
    res['R2.4_scope'] = dict(pairs=int(len(si)), wired=int(Wd[si, sj].sum()),
                             identity_max_abs=float(np.nanmax(np.abs(D['te'][si, sj] - TE1[si, sj]))) if len(si) else None)
    if Wd[si, sj].sum() >= 10:
        res['R2.4'] = run_tests(Wd[si, sj], sj, TE1[si, sj], B[si, sj], D['rt'][si, sj], D['kap'][si, sj],
                                tg['c'][sj], tg['q0'][sj], tg['pi'][sj], ph[sj], mult)
    # R2.5 (A12)
    L = lock_index(X); cls = np.where(L < 0.95, 0, np.where(L < 1, 1, 2))
    keep = off.copy(); keep[clock_g, :] = False; keep[:, clock_g] = False
    ki, kj = np.nonzero(keep & (Wd | ~(Wd | Wd.T))); kpos = Wd[ki, kj]; o = {}; bs = {}
    for c, cn in enumerate(['mixed', 'locked', 'strict']):
        m = cls[kj] == c
        if kpos[m].sum() >= 5:
            r, b = auc_with_ci(TE1[ki, kj][m], kpos[m], kj[m], kj[m], mult); o[cn] = r; bs[cn] = b
            r['targets_with_wired'] = int(len(np.unique(kj[m & kpos])))
    if 'mixed' in bs and 'strict' in bs:
        o['mixed_minus_strict'] = diff_ci(bs['mixed'], bs['strict'], o['mixed']['auc'], o['strict']['auc'])
    o['targets_per_class'] = {cn: int((cls == c).sum()) for c, cn in enumerate(['mixed', 'locked', 'strict'])}
    res['R2.5'] = o
    # R2.6 (A2): latch pairs, stratified by phase class x rate quintiles
    base = ((ph[ii] * 3 + ph[jj]) * 5 + rq[ii]) * 5 + rq[jj]
    res['latch_pairs'] = int(Lat.sum())
    if Lat.sum() >= 5:
        l1, _ = auc_with_ci(np.nan_to_num(TE1[ii, jj], nan=-np.inf), Lat[ii, jj], base, jj, mult)
        l0, _ = auc_with_ci(np.nan_to_num(TE0[ii, jj], nan=-np.inf), Lat[ii, jj], base, jj, mult[:200])
        ld, _ = auc_with_ci(-dist[ii, jj], Lat[ii, jj], base, jj, mult[:200])
        res['R2.6_latch_lag1'] = l1; res['R2.6b_latch_lag0'] = l0; res['latch_distance'] = ld
    # descriptive: share of strictly locked same-phase wired pairs
    res['share_wired_strict_same_phase'] = float(Wd[S].sum() / max(Wd.sum(), 1))
    return res


def met(res):
    r = res['clock_aligned']; m = {}
    m['R2.1'] = bool(r['R2.1_group_auc_lag1']['hi'] < 0.5)
    m['R2.1b'] = bool(r['R2.1b_physical_auc_lag1']['lo'] <= 0.5)
    m['R2.2'] = bool(r['R2.2_lag0_minus_lag1']['lo'] > 0)
    m['R2.3'] = bool(r['R2.3']['share_locked'] > 0.5)
    if 'R2.4' in r:
        t = r['R2.4']
        m['R2.4'] = bool(t['A9.1_followership_within_target']['lo'] > 0.5 and
                         t['A9.2_TE_minus_activity_within_target']['hi'] < 0)
    if 'mixed_minus_strict' in r['R2.5']:
        m['R2.5'] = bool(r['R2.5']['mixed_minus_strict']['lo'] > 0)
    if 'R2.6_latch_lag1' in r:
        m['R2.6'] = bool(r['R2.6_latch_lag1']['lo'] > 0.5)
        m['R2.6b'] = bool(r['R2.6b_latch_lag0']['hi'] < r['R2.6_latch_lag1']['lo'])
    m['R2.7'] = bool(r['R2.7_distance_minus_lag0']['lo'] > 0)
    return m


def prepare(path):
    gate, c1, c2, centre, pullup, names = netlist()
    raw = np.load(path, mmap_mode='r')
    raw, phase, red = reduce(np.asarray(raw))
    start, n_off = clock_start(raw, phase, gate, names)
    reps, p2g = tk.group_identical(red)
    X = red[reps].astype(np.uint8)
    K = np.asarray([tk.ais_select_k(x, 12, base=1.0) for x in X])
    return gate, c1, c2, centre, pullup, names, raw, phase, red, start, n_off, reps, p2g, X, K


def selftest():
    log('selftest: Donkey Kong regenerated, rows 0..401,775')
    raw_path = ROOT / 'data/regen/donkeykong.regen.ts.spikes.npy'
    a = np.load(raw_path, mmap_mode='r')[:401776]
    tmp = OUT / '_dk_selftest.npy'; np.save(tmp, np.asarray(a))
    gate, c1, c2, centre, pullup, names, raw, phase, red, start, n_off, reps, p2g, X, K = prepare(tmp)
    ref = np.load(RUN / 'groups.npz'); Kref = np.load(RUN / 'te.npz')['optimal_k'].astype(int)
    TEref = np.load(RUN / 'te.npz')['te_bits']
    out = dict(samples=int(red.shape[1]), phase=phase, clock_aligned_start=start,
               groups=int(len(reps)), groups_main_run=int(len(ref['representatives'])),
               p2g_equal=bool((p2g == ref['physical_to_group']).all()), K_equal=bool((K == Kref).all()))
    XS = tk.as_sparse_sources(X); d = 0.0
    for t in (0, 100, 300):
        te = tk.te_all_sources(XS, X[t], K[t], delay=1)[0]
        d = max(d, float(np.nanmax(np.abs(np.delete(te, t) - np.delete(TEref[:, t], t)))))
    out['max_abs_lag1_te_vs_main_run'] = d
    t0 = time.time(); tk.te_all_sources(XS, X[5], 12, delay=1); out['seconds_per_target_k12'] = round(time.time() - t0, 2)
    tmp.unlink()
    log(json.dumps(out)); assert out['p2g_equal'] and out['K_equal'] and d < 1e-9 and start == 1408
    log('SELF-TEST PASSED')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--selftest', action='store_true')
    ap.add_argument('game', nargs='?'); ap.add_argument('recording', nargs='?'); a = ap.parse_args()
    if a.selftest:
        return selftest()
    hashes = {p.name: sha(p) for p in (PRED, Path(__file__).resolve(), ROOT / 'analyses/new_tests/common.py')}
    for k, v in hashes.items():
        log(f'R2-{a.game}: {k} sha256 {v}')
    gate, c1, c2, centre, pullup, names, raw, phase, red, start, n_off, reps, p2g, X, K = prepare(a.recording)
    G = len(reps); rate = X.mean(1)
    Wd, Wphys, Nphys, Lat, dist, size, clock_g = references(p2g, G, gate, c1, c2, centre, pullup, rate)
    res = dict(game=a.game, hashes=hashes, raw_rows=int(raw.shape[0]), kept_offset=phase, samples=int(red.shape[1]),
               clock_aligned_start=start, clock_events_off_offset=n_off, groups=int(G), clock_group=clock_g,
               clock_group_rate=float(rate[clock_g]), wired_group_pairs=int(Wd.sum()), K_counts=np.bincount(K).tolist())
    log(f'{a.game}: {res}')
    res['clock_aligned'] = analyse(X[:, start:], K, Wd, Wphys, Nphys, Lat, dist, size, clock_g, 'clock_aligned')
    log(f'{a.game}: clock-aligned done')
    res['whole_recording'] = analyse(X, K, Wd, Wphys, Nphys, Lat, dist, size, clock_g, 'whole')
    res['predictions_met'] = met(res)
    def conv(v):
        if isinstance(v, dict): return {str(k): conv(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)): return [conv(x) for x in v]
        if isinstance(v, (np.integer,)): return int(v)
        if isinstance(v, (np.floating,)): return float(v)
        if isinstance(v, np.bool_): return bool(v)
        return v
    (OUT / f'R2_{a.game}_results.json').write_text(json.dumps(conv(res), indent=1))
    log(f"{a.game}: done {res['predictions_met']}")


if __name__ == '__main__':
    main()
