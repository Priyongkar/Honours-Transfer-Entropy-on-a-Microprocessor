"""A1, A1b, A1c, A5 (chip), A6 (chip) on the simulated 6502 (design: PREDICTIONS.md).
Usage: python chip.py [--quick] [--chip-dir DIR] [--tag NAME]"""
import argparse, glob, json
from pathlib import Path
import numpy as np
from common import (CHIP_DIR, ROOT, tk, log, start, save, change, phase_class, boot_mult,
                    auc_with_ci, verdict, diff_ci, quintile, StratifiedAUC)
from sim6502 import Chip


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--chip-dir', default=str(CHIP_DIR))
    ap.add_argument('--tag', default='chip')
    a = ap.parse_args()
    sha = start(a.tag)
    if a.tag != 'chip':
        import hashlib
        from common import HERE
        log(f"PREDICTIONS_A7.md sha256 {hashlib.sha256((HERE / 'PREDICTIONS_A7.md').read_bytes()).hexdigest()}")
    cdir = Path(a.chip_dir)
    with np.load(cdir / 'states.npz') as z:
        states, nodes, clk0 = z['states'], z['nodes'], z['clk0']
    with np.load(cdir / 'interventions.npz') as z:
        E0, E1, U = z['E0'], z['E1'], z['nodes']
    draws = 50 if a.quick else 1000
    if a.quick:
        states, clk0 = states[:6000], clk0[:6000]
    T, Nn = states.shape; n = len(U)
    chip = Chip(ROOT / 'data/visual6502')
    g_, c1_, c2_ = chip.transistor_matrix()
    gnd, pwr = chip.gnd, chip.pwr
    col = {int(v): i for i, v in enumerate(nodes)}
    A = tk.node_reference_from_transistors(g_, c1_, c2_, nodes, supply={gnd, pwr}).toarray() > 0
    W = A[np.ix_(U, U)]
    off = ~np.eye(n, dtype=bool)
    i0 = E0[:, U] > 0; i1 = E1[:, U] > 0
    np.fill_diagonal(i0, False); np.fill_diagonal(i1, False)
    D = i1 & ~i0 & off
    I = i0 & off
    Nn_ = ~i0 & ~i1 & ~i0.T & ~i1.T & off
    Ci = i0.astype(np.int32); common = Ci.T @ Ci; deg = Ci.sum(0)
    union = (deg[:, None] - i0.T) + (deg[None, :] - i0) - common
    jacc = np.where(union > 0, common / np.maximum(union, 1), 0.0)
    pos_in_U = {int(u): k for k, u in enumerate(U)}

    # latch pairs (node level) and pull-down edges, in U coordinates
    st_all = states.T
    clock_cols = {c for c in range(Nn) if (st_all[c] == clk0).all() or (st_all[c] == 1 - clk0).all()}
    pd_nodes = {int(c1_[t]) for t in range(len(g_)) if c2_[t] == gnd}
    driven = lambda nd: bool(chip.pullup[nd]) or nd in pd_nodes
    latch = set(); latch_unor = set(); pulldown = set()
    for t in range(len(g_)):
        g, c1, c2 = int(g_[t]), int(c1_[t]), int(c2_[t])
        if col.get(g) in clock_cols and c1 not in (gnd, pwr) and c2 not in (gnd, pwr):
            if c1 in col and c2 in col and col[c1] in pos_in_U and col[c2] in pos_in_U:
                a1, a2 = pos_in_U[col[c1]], pos_in_U[col[c2]]
                d1, d2 = driven(c1), driven(c2)
                if d1 != d2:
                    latch.add((a1, a2) if d1 else (a2, a1))
                else:
                    latch_unor.update({(a1, a2), (a2, a1)})
        if c2 == gnd and g in col and c1 in col and col[g] in pos_in_U and col[c1] in pos_in_U:
            pulldown.add((pos_in_U[col[g]], pos_in_U[col[c1]]))
    log(f'clock nodes {len(clock_cols)}; oriented latch pairs {len(latch)}; unoriented {len(latch_unor) // 2}; pull-down edges {len(pulldown)}')
    L = np.zeros((n, n), bool);
    for p in latch: L[p] = True
    PD = np.zeros((n, n), bool)
    for p in pulldown: PD[p] = True

    st = st_all[U].astype(np.uint8); chg = change(st)
    parity = (np.arange(T) % 2).astype(np.uint8)
    ph, _ = phase_class(chg); rate = chg.mean(1); rq = quintile(rate)
    mult = boot_mult(n, draws)
    res = dict(predictions_sha256=sha, T=T, quick=a.quick,
               class_sizes=dict(delayed_only=int(D.sum()), immediate=int(I.sum()), none=int(Nn_.sum()),
                                latch_oriented=len(latch), pulldown=len(pulldown), clock_nodes=len(clock_cols)))
    log(f'classes D {D.sum()}, I {I.sum()}, N {Nn_.sum()}')
    keep = {}
    for rep, X in [('change', chg), ('state', st)]:
        reps, ug = tk.group_identical(X); G = len(reps); Xg = X[reps]
        ks = np.asarray([tk.ais_select_k(x, 6) for x in Xg])
        S = {'lag1': tk.te_matrix(Xg, ks, delay=1)[0], 'lag0': tk.te_matrix(Xg, ks, delay=0)[0]}
        if rep == 'change':
            S['lag1_parity'] = tk.te_matrix(Xg, ks, delay=1, cond=lambda t: parity[None, :])[0]
        log(f'{rep}: {G} groups; scores done')
        valid = (ug[:, None] >= 0) & (ug[None, :] >= 0) & (ug[:, None] != ug[None, :]) & off
        ii, jj = np.nonzero(valid)
        base = ((ph[ii] * 3 + ph[jj]) * 5 + rq[ii]) * 5 + rq[jj]
        strat = base * 5 + quintile(jacc[ii, jj])
        sc = {k: v[ug[ii], ug[jj]] for k, v in S.items()}
        cD, cI, cN, cL, cPD = D[ii, jj], I[ii, jj], Nn_[ii, jj], L[ii, jj], PD[ii, jj]
        cLr = L.T[ii, jj]
        out = dict(groups=G, scorable=dict(D=int(cD.sum()), I=int(cI.sum()), N=int(cN.sum()), latch=int(cL.sum())))
        def run(name, s, pos, neg, stt):
            m = pos | neg
            r, bs = auc_with_ci(s[m], pos[m], stt[m], jj[m], mult)
            r['verdict'] = verdict(r); out[name] = r; return r, bs
        r1, b1 = run('A1_P1_lag1_D_vs_N', sc['lag1'], cD, cN, strat)
        r0, b0 = run('A1_lag0_D_vs_N', sc['lag0'], cD, cN, strat)
        out['A1_P2_lag1_minus_lag0_D'] = diff_ci(b1, b0, r1['auc'], r0['auc'])
        rI0, bI0 = run('A1_lag0_I_vs_N', sc['lag0'], cI, cN, strat)
        rI1, bI1 = run('A1_lag1_I_vs_N', sc['lag1'], cI, cN, strat)
        out['A1_P3_lag0_minus_lag1_I'] = diff_ci(bI0, bI1, rI0['auc'], rI1['auc'])
        unw = ~(W | W.T)[ii, jj]
        run('C2_lag1_wired_vs_unwired', sc['lag1'], W[ii, jj], unw, base)
        run('C2_lag0_wired_vs_unwired', sc['lag0'], W[ii, jj], unw, base)
        if 'lag1_parity' in sc:
            run('A1_P4_lag1parity_D_vs_N', sc['lag1_parity'], cD, cN, strat)
        if cL.sum() > 0:
            run('A1b_lag1_latch_vs_N', sc['lag1'], cL, cN & ~cL, base)
            run('A1b_lag0_latch_vs_N', sc['lag0'], cL, cN & ~cL, base)
            run('A1b_lag1_reverse_latch_vs_N', sc['lag1'], cLr, cN & ~cLr, base)
            out['A1b_share_latch_in_D'] = float(cD[cL].mean())
        res[rep] = out
        keep[rep] = (Xg, ks, ug, ii, jj, base, cL, cPD, cN, S)
        log(f"{rep}: P1 {out['A1_P1_lag1_D_vs_N']['auc']:.3f} [{out['A1_P1_lag1_D_vs_N']['lo']:.3f}, {out['A1_P1_lag1_D_vs_N']['hi']:.3f}]")

    # A5 (chip): delay profile, change coding
    Xg, ks, ug, ii, jj, base, cL, cPD, cN, S = keep['change']
    prof = []
    for u in range(1, 9):
        Su = S['lag1'] if u == 1 else tk.te_matrix(Xg, ks, delay=u)[0]
        s = Su[ug[ii], ug[jj]]
        row = dict(delay=u)
        for name, pos in [('latch', cL), ('pulldown', cPD)]:
            m = pos | (cN & ~pos)
            if pos.sum():
                r, _ = auc_with_ci(s[m], pos[m], base[m], jj[m], mult[:200])
                row[name] = r
        prof.append(row)
    res['A5_chip_delay_profile'] = prof
    log('A5 chip done')

    # A6 (chip): orientation of oriented latch pairs, lag-1 TE, change coding
    te1 = S['lag1']; H = np.asarray([tk.conditional_entropy_bits(x, k) for x, k in zip(Xg, ks)])
    rows = []
    for d, s_ in latch:
        gd, gs = ug[d], ug[s_]
        if gd < 0 or gs < 0 or gd == gs:
            continue
        f, b = te1[gd, gs], te1[gs, gd]
        rows.append((1.0 if f > b else 0.5 if f == b else 0.0,
                     1.0 if rate[s_] < rate[d] else 0.5 if rate[s_] == rate[d] else 0.0,
                     1.0 if H[gs] > H[gd] else 0.5 if H[gs] == H[gd] else 0.0))
    res['A6_chip'] = orientation(rows)
    log(f"A6 chip: {res['A6_chip']}")

    # A1c: IDTxl chip selections (states)
    name = f'{a.tag}_results' + ('_quick' if a.quick else '')
    save(name, res)
    if cdir == Path(CHIP_DIR):
        try:
            res['A1c'] = idtxl_chip(states if not a.quick else np.load(cdir / 'states.npz')['states'], E0, E1, U, A, cdir)
            log(f"A1c done: {res['A1c'].get('grouping_matches_idtxl_data')}")
        except Exception as e:
            res['A1c'] = dict(error=repr(e)); log(f'A1c failed: {e!r}')
    save(name, res)
    log('done')


def orientation(rows):
    if not rows:
        return dict(pairs=0)
    R = np.asarray(rows); rng = np.random.default_rng(2026)
    bs = np.asarray([R[rng.integers(0, len(R), len(R))].mean(0) for _ in range(1000)])
    q = lambda v: [float(np.quantile(v, .025)), float(np.quantile(v, .975))]
    return dict(pairs=len(R), te_accuracy=float(R[:, 0].mean()), te_ci=q(bs[:, 0]),
                rate_rule=float(R[:, 1].mean()), entropy_rule=float(R[:, 2].mean()),
                te_minus_rate=float(R[:, 0].mean() - R[:, 1].mean()), te_minus_rate_ci=q(bs[:, 0] - bs[:, 1]),
                te_minus_entropy=float(R[:, 0].mean() - R[:, 2].mean()), te_minus_entropy_ci=q(bs[:, 0] - bs[:, 2]))


def idtxl_chip(states, E0, E1, U, A, cdir):
    S = states.T.astype(np.uint8)
    reps, ug = tk.group_identical(S)
    data = np.load(cdir / 'idtxl_data_state.npy', mmap_mode='r')
    match = data.shape[0] == len(reps) and np.array_equal(np.asarray(data[:, :200]).astype(np.uint8), S[reps][:, :200])
    members = [np.flatnonzero(ug == g) for g in range(len(reps))]
    Nn = S.shape[0]
    D0 = np.zeros((Nn, Nn), bool); D1 = np.zeros((Nn, Nn), bool); has = np.zeros(Nn, bool)
    D0[U] = E0 > 0; D1[U] = E1 > 0; has[U] = True
    Dd = D1 & ~D0; np.fill_diagonal(Dd, False); np.fill_diagonal(D0, False)
    te1 = np.load(cdir.parent / 'supplementary/sim_state_te.npz')['te1']
    def edges_prec(edges):
        tot = dict(pairs=0, D=0, I=0, wired=0, wired_pairs=0, all_pairs=0)
        for s, t in edges:
            src, tgt = members[s], members[t]
            blk_w = A[np.ix_(src, tgt)]
            tot['wired'] += int(blk_w.sum()); tot['all_pairs'] += blk_w.size
            sv = src[has[src]]
            if len(sv):
                tot['pairs'] += len(sv) * len(tgt)
                tot['D'] += int(Dd[np.ix_(sv, tgt)].sum()); tot['I'] += int(D0[np.ix_(sv, tgt)].sum())
        p = max(tot['pairs'], 1); q = max(tot['all_pairs'], 1)
        return dict(edges=len(edges), precision_D=tot['D'] / p, precision_I=tot['I'] / p, precision_wired=tot['wired'] / q)
    out = dict(grouping_matches_idtxl_data=bool(match))
    for name in ['idtxl_state_multivariate', 'idtxl_state_multivariate_lag0']:
        files = sorted(glob.glob(str(cdir / name / 'target_*.json')))
        e_all, e_l0, e_l1, comp = [], [], [], []
        for f in files:
            j = json.load(open(f)); t = j['target']
            srcs = j['selected_vars_sources']
            e_all += [(s, t) for s, _ in srcs]; e_l0 += [(s, t) for s, l in srcs if l == 0]; e_l1 += [(s, t) for s, l in srcs if l == 1]
            k = len({s for s, _ in srcs})
            if k:
                col = np.where(np.arange(te1.shape[0]) == t, -np.inf, np.nan_to_num(te1[:, t], nan=-np.inf))
                comp += [(int(s), t) for s in np.argsort(-col)[:k]]
        out[name] = dict(targets=len(files), all=edges_prec(e_all), lag0=edges_prec(e_l0), lag1=edges_prec(e_l1),
                         pairwise_lag1_same_allocation=edges_prec(comp))
    return out


if __name__ == '__main__':
    main()
