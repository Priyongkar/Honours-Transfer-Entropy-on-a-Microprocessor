"""A2, A3, A5 (recording), A6 (recording) on the reduced recording (design: PREDICTIONS.md).
Usage: python recording.py [--quick]"""
import argparse, ast
import numpy as np
from sklearn.metrics import average_precision_score
from common import (RUN, ROOT, OUT, tk, log, start, save, phase_class, boot_mult, auc_with_ci, verdict,
                    quintile, bound_all_sources)
from chip import orientation
from sim6502 import load_netlist

VSS, VCC = 558, 657


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); a = ap.parse_args()
    sha = start('recording')
    X = np.asarray(np.load(RUN / 'thesis_v5/idtxl/data_int8.npy')).astype(np.uint8)
    with np.load(RUN / 'te.npz') as z:
        TE1 = z['te_bits'].copy(); K = z['optimal_k'].astype(int)
    TE0 = np.load(RUN / 'thesis_v5/te_lag0.npz')['te_bits']
    p2g = np.load(RUN / 'groups.npz')['physical_to_group'].astype(int)
    NET = np.load(RUN / 'thesis_v5/idtxl/network_multivariate.npz')['adjacency'].astype(bool)
    G, T = X.shape
    XS = tk.as_sparse_sources(X)
    draws = 50 if a.quick else 1000
    # quick: checks the code on a subset of targets; statistics are not meant to be read
    tsub = np.arange(G) if not a.quick else np.arange(0, G, 10)
    # sanity: recomputed lag-1 TE equals the saved matrix
    chk = [int(t) for t in tsub[:3]]
    for t in chk:
        te, _ = tk.te_all_sources(XS, X[t], K[t], delay=1)
        d = np.nanmax(np.abs(np.delete(te, t) - np.delete(TE1[:, t], t)))
        assert d < 1e-9, f'saved lag-1 TE differs for target {t}: {d}'
    log('saved lag-1 TE reproduced')

    # netlist
    defs = []
    for line in (ROOT / 'data/visual6502/transdefs.js').read_text().splitlines():
        if line.strip().startswith("['t"):
            defs.append(ast.literal_eval(line.strip().rstrip(',')))
    defs.sort(key=lambda d: int(d[0][1:]))
    gate, c1, c2 = (np.asarray([d[i] for d in defs]) for i in (1, 2, 3))
    box = np.asarray([d[4] for d in defs], float)
    centre = np.column_stack(((box[:, 0] + box[:, 1]) / 2, (box[:, 2] + box[:, 3]) / 2))
    pullup = load_netlist(str(ROOT / 'data/visual6502'))[0]
    gated_by = {}
    for t, g in enumerate(gate):
        gated_by.setdefault(int(g), []).append(t)
    pd_nodes = {int(c1[t]) for t in range(len(gate)) if c2[t] == VSS} | {int(c2[t]) for t in range(len(gate)) if c1[t] == VSS}
    driven = lambda nd: bool(pullup.get(nd, False)) or nd in pd_nodes
    Wd = np.zeros((G, G), bool); PDw = np.zeros((G, G), bool)
    for s in range(len(gate)):
        if p2g[s] < 0:
            continue
        for term in (c1[s], c2[s]):
            if term in (VSS, VCC):
                continue
            for t in gated_by.get(int(term), []):
                if p2g[t] >= 0 and p2g[t] != p2g[s]:
                    Wd[p2g[s], p2g[t]] = True
                    if VSS in (c1[s], c2[s]):
                        PDw[p2g[s], p2g[t]] = True
    rate = X.mean(1); clock_g = int(np.argmax(rate))
    clock_gates = {int(gate[t]) for t in np.flatnonzero(p2g == clock_g)}
    groups_of = lambda nd: {int(p2g[t]) for t in gated_by.get(int(nd), []) if p2g[t] >= 0}
    Lat = np.zeros((G, G), bool); Unor = np.zeros((G, G), bool)
    for t in range(len(gate)):
        if int(gate[t]) in clock_gates and c1[t] not in (VSS, VCC) and c2[t] not in (VSS, VCC):
            a1, a2 = int(c1[t]), int(c2[t]); d1, d2 = driven(a1), driven(a2)
            for ga in groups_of(a1):
                for gb in groups_of(a2):
                    if ga == gb:
                        continue
                    if d1 != d2:
                        if d1: Lat[ga, gb] = True
                        else: Lat[gb, ga] = True
                    else:
                        Unor[ga, gb] = Unor[gb, ga] = True
    log(f'clock group {clock_g} (rate {rate[clock_g]:.3f}), clock gate nodes {len(clock_gates)}; oriented latch group pairs {Lat.sum()}, '
        f'unoriented {Unor.sum()}, of oriented wired {int((Lat & (Wd | Wd.T)).sum())}')

    ph, strict = phase_class(X); rq = quintile(rate)
    np.savez_compressed(OUT / 'recording_reference.npz', wired=Wd, pulldown=PDw, latch=Lat, phase=ph, rate=rate, k=K)
    off = ~np.eye(G, dtype=bool)
    tmask = np.zeros(G, bool); tmask[tsub] = True
    ii, jj = np.nonzero(off & tmask[None, :])
    base = ((ph[ii] * 3 + ph[jj]) * 5 + rq[ii]) * 5 + rq[jj]
    mult = boot_mult(G, draws)
    size = np.bincount(p2g[p2g >= 0], minlength=G).astype(float)
    act = np.flatnonzero(p2g >= 0); order = act[np.argsort(p2g[act])]
    bounds = np.r_[0, np.cumsum(np.bincount(p2g[order], minlength=G))]
    from scipy.spatial.distance import cdist
    dist = np.zeros((G, G))
    for g in range(G):
        dmin = cdist(centre[order[bounds[g]:bounds[g + 1]]], centre[order]).min(0)
        dist[g] = np.minimum.reduceat(dmin, bounds[:-1])
    res = dict(predictions_sha256=sha, quick=a.quick, clock_group=clock_g,
               latch_pairs=dict(oriented=int(Lat.sum()), unoriented=int(Unor.sum()),
                                oriented_wired=int((Lat & (Wd | Wd.T)).sum()), thesis_count_clock_gated_pass_neighbours=82))

    # A2
    def a2_block(S1, S0, label):
        scores = {'lag1_TE': S1, 'lag0_CMI': S0, 'IDTxl_lag1_network': NET.astype(float), 'source_group_size': np.repeat(size[:, None], G, 1),
                  'minus_layout_distance': -dist, 'source_event_rate': np.repeat(rate[:, None], G, 1)}
        pos = Lat[ii, jj]; out = {}
        kt = NET.sum(0)
        for name, S in scores.items():
            s = np.nan_to_num(S[ii, jj], nan=-np.inf)
            r, _ = auc_with_ci(s, pos, base, jj, mult if name == 'lag1_TE' else mult[:200])
            r['verdict'] = verdict(r)
            fin = np.isfinite(s)
            r['AP_lift'] = float(average_precision_score(pos[fin], s[fin]) / pos[fin].mean())
            hits = tot = 0
            for t in tsub:
                if kt[t] == 0: continue
                col = np.where(np.arange(G) == t, -np.inf, np.nan_to_num(S[:, t], nan=-np.inf))
                top = np.argsort(-col)[:kt[t]]
                hits += int(Lat[top, t].sum()); tot += len(top)
            r['precision_at_IDTxl_allocation'] = hits / max(tot, 1)
            out[name] = r
        out['prevalence'] = float(pos.mean())
        res[f'A2_{label}'] = out
    a2_block(TE1, TE0, 'whole_recording')
    log('A2 whole recording done')
    Xa = X[:, 1407:]; XaS = tk.as_sparse_sources(Xa)
    TE1a = np.full((G, G), np.nan); TE0a = np.full((G, G), np.nan)
    for t in tsub:
        TE1a[:, t] = tk.te_all_sources(XaS, Xa[t], K[t], delay=1)[0]
        TE0a[:, t] = tk.te_all_sources(XaS, Xa[t], K[t], delay=0)[0]
    np.fill_diagonal(TE1a, np.nan); np.fill_diagonal(TE0a, np.nan)
    a2_block(TE1a, TE0a, 'clock_aligned')
    log('A2 clock-aligned done')
    save('recording_results' + ('_quick' if a.quick else ''), res)

    # A3
    P = (np.arange(T) % 2).astype(np.uint8)
    B = np.full((G, G), np.nan)
    for t in tsub:
        B[:, t] = bound_all_sources(XS, X[t], K[t], P, delay=1)
    np.fill_diagonal(B, np.nan)
    same = (ph[:, None] == ph[None, :]) & (ph[:, None] < 2) & off & tmask[None, :]
    strict_same = same & strict[:, None] & strict[None, :]
    exc_strict = TE1[strict_same] - B[strict_same]
    exc_tol = TE1[same & ~strict_same] - B[same & ~strict_same]
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = TE1 / B
    si, sj = np.nonzero(same)
    bd = quintile(B[si, sj]) if False else np.searchsorted(np.unique(np.quantile(B[si, sj], np.linspace(.1, .9, 9))), B[si, sj], side='right')
    wired_same = (Wd | False)[si, sj]
    r3, _ = auc_with_ci(TE1[si, sj], wired_same, bd, sj, mult[:200])
    res['A3'] = dict(strict_same_pairs=int(strict_same.sum()), max_excess_strict=float(np.nanmax(exc_strict)) if exc_strict.size else None,
                     tolerance_pairs=int(exc_tol.size), tolerance_pairs_with_excess=int((exc_tol > 1e-12).sum()),
                     max_excess_tolerance=float(np.nanmax(exc_tol)) if exc_tol.size else None,
                     median_ratio_same_phase=float(np.nanmedian(ratio[same])),
                     median_ratio_by_class={f'{p}{q}': float(np.nanmedian(ratio[(ph[:, None] == p) & (ph[None, :] == q) & off & tmask[None, :]]))
                                            for p in range(3) for q in range(3)},
                     wired_vs_unwired_within_B_deciles=r3)
    np.savez_compressed(OUT / 'A3_bound_arrays.npz', TE=TE1, B=B, wired=Wd, phase=ph, strict=strict)
    log(f"A3 done: max excess strict {res['A3']['max_excess_strict']}")
    save('recording_results' + ('_quick' if a.quick else ''), res)

    # A5 (recording)
    N_ = ~(Wd | Wd.T) & ~(Lat | Lat.T) & off
    prof = []
    for u in range(1, 9):
        if u == 1:
            Su = TE1
        else:
            Su = np.full((G, G), np.nan)
            for t in tsub:
                Su[:, t] = tk.te_all_sources(XS, X[t], K[t], delay=u)[0]
        s = np.nan_to_num(Su[ii, jj], nan=-np.inf); row = dict(delay=u)
        for name, Pm in [('latch', Lat), ('pulldown_wired', PDw)]:
            pos = Pm[ii, jj]; m = pos | N_[ii, jj]
            if pos.sum():
                r, _ = auc_with_ci(s[m], pos[m], base[m], jj[m], mult[:200]); row[name] = r
        prof.append(row)
    res['A5_recording_delay_profile'] = prof
    log('A5 recording done')

    # A6 (recording)
    H = np.asarray([tk.conditional_entropy_bits(X[g], K[g]) for g in range(G)])
    rows = []
    for d, s_ in zip(*np.nonzero(Lat)):
        if not (tmask[s_] and tmask[d]) and not a.quick:
            continue
        f, b = TE1[d, s_], TE1[s_, d]
        rows.append((1.0 if f > b else 0.5 if f == b else 0.0,
                     1.0 if rate[s_] < rate[d] else 0.5 if rate[s_] == rate[d] else 0.0,
                     1.0 if H[s_] > H[d] else 0.5 if H[s_] == H[d] else 0.0))
    res['A6_recording'] = orientation(rows)
    save('recording_results' + ('_quick' if a.quick else ''), res)
    log('done')


if __name__ == '__main__':
    main()
