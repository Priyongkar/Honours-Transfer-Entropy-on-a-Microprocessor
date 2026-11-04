"""A13 (chip), A14 (chip) and A16 (timed clamps) on both simulated programs
(design: PREDICTIONS_A12.md).

Usage:
  python a13_a16_chip.py --selftest   timing only (warm-up, a few clamps); prints no statistic
  python a13_a16_chip.py              full run, both programs
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / 'src'))  # repository library
from te6502.paths import ANALYSES, DATA_ROOT, RESULTS, HELPERS_FILE, load_analysis  # noqa: E402,F401
import argparse, hashlib, time
from pathlib import Path
import numpy as np
from te6502.test_helpers import (HERE, ROOT, OUT, CHIP_DIR, tk, log, save, change, phase_class, boot_mult,
                    auc_with_ci, verdict, diff_ci, quintile, bound_all_sources)
a9 = load_analysis('04_followership_mechanism_A9_to_A11/followership_decomposition_A9.py')
from te6502.switch_level_simulator import Chip, PROGRAM_FIB_LFSR

PRED = ANALYSES / '05_consequences_of_the_account_A12_to_A19/design/predictions_A12_to_A19.md'
SEED, DRAWS = 2026, 1000


def start(tag):
    h = {}
    for p in (PRED, Path(__file__).resolve(), HELPERS_FILE, ANALYSES / '04_followership_mechanism_A9_to_A11/followership_decomposition_A9.py'):
        h[p.name] = hashlib.sha256(Path(p).read_bytes()).hexdigest()
        log(f'{tag}: {p.name} sha256 {h[p.name]}')
    return h


def programs():
    a7_program = load_analysis('02_latch_transfer_tests_A1_to_A8/second_program_build_simulate_and_test_A7.py')
    code2, _ = a7_program.assemble()
    return [('program1', CHIP_DIR, PROGRAM_FIB_LFSR, 0x0000), ('program2', OUT / 'A7_chip', code2, a7_program.ORIGIN)]


def auc(out, name, score, pos, strata, target, mult):
    r, b = auc_with_ci(np.asarray(score, float), np.asarray(pos, bool), np.asarray(strata), np.asarray(target), mult)
    r['verdict'] = verdict(r); out[name] = r
    return r, b


def build(cdir):
    with np.load(Path(cdir) / 'states.npz') as z:
        states, nodes, clk0 = z['states'], z['nodes'], z['clk0']
    with np.load(Path(cdir) / 'interventions.npz') as z:
        E0, E1, U = z['E0'], z['E1'], z['nodes']
    T, Nn = states.shape; n = len(U)
    chip = Chip(ROOT / 'data/visual6502')
    g_, c1_, c2_ = chip.transistor_matrix(); gnd, pwr = chip.gnd, chip.pwr
    col = {int(v): i for i, v in enumerate(nodes)}
    A = tk.node_reference_from_transistors(g_, c1_, c2_, nodes, supply={gnd, pwr}).toarray() > 0
    W = A[np.ix_(U, U)]
    off = ~np.eye(n, dtype=bool)
    i0 = E0[:, U] > 0; i1 = E1[:, U] > 0; np.fill_diagonal(i0, False); np.fill_diagonal(i1, False)
    Nn_ = ~i0 & ~i1 & ~i0.T & ~i1.T & off
    Ci = i0.astype(np.int32); common = Ci.T @ Ci; deg = Ci.sum(0)
    union = (deg[:, None] - i0.T) + (deg[None, :] - i0) - common
    jacc = np.where(union > 0, common / np.maximum(union, 1), 0.0)
    pos_in_U = {int(u): k for k, u in enumerate(U)}
    st_all = states.T
    pol = {}
    for c in range(Nn):
        if (st_all[c] == clk0).all(): pol[c] = 1
        elif (st_all[c] == 1 - clk0).all(): pol[c] = -1
    pd_nodes = {int(c1_[t]) for t in range(len(g_)) if c2_[t] == gnd}
    driven = lambda nd: bool(chip.pullup[nd]) or nd in pd_nodes
    latch = {}
    for t in range(len(g_)):
        g, c1, c2 = int(g_[t]), int(c1_[t]), int(c2_[t])
        if col.get(g) in pol and c1 not in (gnd, pwr) and c2 not in (gnd, pwr):
            if c1 in col and c2 in col and col[c1] in pos_in_U and col[c2] in pos_in_U:
                a1, a2 = pos_in_U[col[c1]], pos_in_U[col[c2]]
                d1, d2 = driven(c1), driven(c2)
                if d1 != d2:
                    latch[(a1, a2) if d1 else (a2, a1)] = pol[col[g]]
    L = np.zeros((n, n), bool)
    for p_ in latch: L[p_] = True
    chg = change(st_all[U].astype(np.uint8))
    ph, _ = phase_class(chg); rate = chg.mean(1); rq = quintile(rate)
    reps, ug = tk.group_identical(chg); Xg = chg[reps]
    ks = np.asarray([tk.ais_select_k(x, 6) for x in Xg])
    valid = (ug[:, None] >= 0) & (ug[None, :] >= 0) & (ug[:, None] != ug[None, :]) & off
    return dict(states=states, nodes=nodes, clk0=clk0, U=U, n=n, W=W, L=L, latch=latch, Nn=Nn_, i0=i0, i1=i1,
                jacc=jacc, ph=ph, rq=rq, ug=ug, Xg=Xg, ks=ks, valid=valid, T=T)


def a13(d, mult):
    ii, jj = np.nonzero(d['valid'] & (d['L'] | (d['Nn'] & ~d['L'])))
    base = ((d['ph'][ii] * 3 + d['ph'][jj]) * 5 + d['rq'][ii]) * 5 + d['rq'][jj]
    ug, ks = d['ug'], d['ks']
    S_all, P_all, ST_all, TG_all, COV_all = [], [], [], [], []
    for u in range(2, 8):
        Su = tk.te_matrix(d['Xg'], ks, delay=u)[0]
        S_all.append(Su[ug[ii], ug[jj]]); P_all.append(d['L'][ii, jj]); ST_all.append(base * 10 + u)
        TG_all.append(jj); COV_all.append(ks[ug[jj]] >= u - 1)
    s, p_, stt, tg, cov = (np.concatenate(v) for v in (S_all, P_all, ST_all, TG_all, COV_all))
    o = dict(covered_latch_pairs=int(p_[cov].sum()), uncovered_latch_pairs=int(p_[~cov].sum()))
    rc, bc = auc(o, 'A13.1_covered', s[cov], p_[cov], stt[cov], tg[cov], mult)
    pm = {'A13.1': bool(rc['hi'] < 0.5)}
    if p_[~cov].sum() > 0:
        ru, bu = auc(o, 'uncovered', s[~cov], p_[~cov], stt[~cov], tg[~cov], mult)
        dd = diff_ci(bu, bc, ru['auc'], rc['auc']); o['A13.2_uncovered_minus_covered'] = dd
        pm['A13.2'] = bool(dd['lo'] > 0); pm['A13.3'] = bool(ru['hi'] >= 0.5)
    o['descriptive_only'] = bool(p_[~cov].sum() < 10); o['predictions_met'] = pm
    return o


def a14(d, mult):
    Xg, ks = d['Xg'], d['ks']; Gn = len(Xg); T = Xg.shape[1]; par = (np.arange(T) % 2).astype(np.int64)
    phg, strg = phase_class(Xg); Xsp = tk.as_sparse_sources(Xg)
    offg = ~np.eye(Gn, dtype=bool)
    Sg = (phg[:, None] == phg[None, :]) & (phg[:, None] < 2) & strg[:, None] & strg[None, :] & offg
    targets = np.flatnonzero(Sg.any(0))
    TEg = np.full((Gn, Gn), np.nan); Bg = np.full((Gn, Gn), np.nan)
    for j in targets:
        TEg[:, j] = tk.te_all_sources(Xsp, Xg[j], ks[j], delay=1)[0]
        Bg[:, j] = bound_all_sources(Xsp, Xg[j], ks[j], par, delay=1)
    D, tg = a9.decompose(Xg, Xsp, ks, phg, targets, par)
    ug = d['ug']; W = d['W']
    ii, jj = np.nonzero(d['valid']); gi, gj = ug[ii], ug[jj]
    keep = Sg[gi, gj] & (W[ii, jj] | ~(W | W.T)[ii, jj])
    ii, jj, gi, gj = ii[keep], jj[keep], gi[keep], gj[keep]; pos = W[ii, jj]
    bd = a9.deciles(Bg[gi, gj]); q0q = quintile(tg['q0'][gj]); piq = quintile(tg['pi'][gj])
    o = {}
    rB, bB = auc(o, 'TE_within_B_deciles', TEg[gi, gj], pos, bd, jj, mult)
    rX, bX = auc(o, 'TE_within_B_q0_pi', TEg[gi, gj], pos, (bd * 5 + q0q) * 5 + piq, jj, mult)
    dd = diff_ci(bX, bB, rX['auc'], rB['auc']); o['difference'] = dd
    inc = bool(rX['lo'] <= 0.5 <= rX['hi'])
    o['reading'] = 'accounted for' if (dd['lo'] > 0 and inc) else ('partly accounted for' if dd['lo'] > 0 else 'not accounted for')
    return o


def timed_clamps(d, code, origin, n_moments=40, limit_nodes=None):
    """E0/E1 per clock-phase set; returns dict with 'high' and 'low' (clk0 at the snapshot)."""
    clk0 = d['clk0']; T = d['T']; rec = d['nodes']; U = d['U']
    cand = np.arange(500, T - 2)
    hi_c = cand[clk0[cand - 1] == 1]; lo_c = cand[clk0[cand - 1] == 0]
    times = {'high': np.sort(np.random.default_rng(3).choice(hi_c, n_moments, replace=False)),
             'low': np.sort(np.random.default_rng(4).choice(lo_c, n_moments, replace=False))}
    allt = sorted([(t, s) for s in times for t in times[s]])
    chip = Chip(ROOT / 'data/visual6502'); chip.reset(program=code, origin=origin); chip.run(200)
    state_of = lambda: np.fromiter((chip.state[x] for x in rec), dtype=bool, count=len(rec))
    snaps, prev, mism = [], 0, 0
    for t, s in allt:
        for _ in range(t - prev): chip.half_step()
        prev = t; snaps.append((s, chip.snapshot()))
        mism += int((state_of() != d['states'][t - 1].astype(bool)).sum())
    nodes = U if limit_nodes is None else U[:limit_nodes]
    out = {s: dict(E0=np.zeros((len(nodes), len(rec)), np.float32), E1=np.zeros((len(nodes), len(rec)), np.float32)) for s in times}
    t0 = time.perf_counter()
    for k, (s, snap) in enumerate(snaps):
        chip.restore(snap); chip.clamp = {}; b0 = state_of(); chip.half_step(); b1 = state_of()
        for i, v in enumerate(nodes):
            chip.restore(snap); chip.clamp = {}
            node = int(rec[v]); chip.set_clamp(node, not chip.state[node])
            out[s]['E0'][i] += state_of() != b0
            chip.half_step(); out[s]['E1'][i] += state_of() != b1
            chip.clamp = {}
        if (k + 1) % 10 == 0:
            log(f'  timed clamps: moment {k + 1}/{len(snaps)} ({time.perf_counter() - t0:.0f}s)')
    for s in out:
        out[s]['E0'] /= n_moments; out[s]['E1'] /= n_moments
    return out, mism, time.perf_counter() - t0


def a16(d, cl, mult):
    U, n = d['U'], d['n']; off = ~np.eye(n, dtype=bool)
    I0 = {s: (cl[s]['E0'][:, U] > 0) & off for s in cl}; I1 = {s: (cl[s]['E1'][:, U] > 0) & off for s in cl}
    delayed = {s: I1[s] & ~I0[s] for s in cl}
    rows = []
    for (a, b), p in d['latch'].items():
        s = 'low' if p == 1 else 'high'          # the set whose moments find the gate low (latch closed)
        rows.append(bool(delayed[s][a, b]))
    rows = np.asarray(rows, bool)
    rng = np.random.default_rng(SEED)
    bs = [rows[rng.integers(0, len(rows), len(rows))].mean() for _ in range(1000)]
    o = dict(latch_pairs=int(len(rows)), share_delayed_only_closed=float(rows.mean()),
             share_lo=float(np.quantile(bs, .025)), share_hi=float(np.quantile(bs, .975)))
    imm_any = I0['high'] | I0['low']; del_any = (I1['high'] | I1['low']) & ~imm_any
    none = ~I0['high'] & ~I0['low'] & ~I1['high'] & ~I1['low']; none = none & none.T & off
    o['class_sizes'] = dict(delayed_only=int(del_any.sum()), immediate=int(imm_any.sum()), none=int(none.sum()),
                            latch_in_delayed_only=int(sum(del_any[a, b] for a, b in d['latch'])))
    ug, ks, Xg = d['ug'], d['ks'], d['Xg']
    L1 = tk.te_matrix(Xg, ks, delay=1)[0]; L0 = tk.te_matrix(Xg, ks, delay=0)[0]
    ii, jj = np.nonzero(d['valid'] & (del_any | none))
    base = ((d['ph'][ii] * 3 + d['ph'][jj]) * 5 + d['rq'][ii]) * 5 + d['rq'][jj]
    strat = base * 5 + quintile(d['jacc'][ii, jj]); pos = del_any[ii, jj]
    r1, b1 = auc(o, 'A16.2_lag1_delayed_vs_none', L1[ug[ii], ug[jj]], pos, strat, jj, mult)
    r0, b0 = auc(o, 'lag0_delayed_vs_none', L0[ug[ii], ug[jj]], pos, strat, jj, mult)
    dd = diff_ci(b1, b0, r1['auc'], r0['auc']); o['A16.3_lag1_minus_lag0'] = dd
    o['predictions_met'] = {'A16.1': bool(o['share_lo'] > 0.5), 'A16.2': bool(r1['lo'] > 0.5), 'A16.3': bool(dd['lo'] > 0)}
    return o


def main_run():
    hashes = start('A13-A16-chip')
    res = dict(hashes=hashes, draws=DRAWS)
    for name, cdir, code, origin in programs():
        d = build(cdir); mult = boot_mult(d['n'], DRAWS, SEED)
        log(f'{name}: {len(d["latch"])} latch pairs; {len(d["Xg"])} groups')
        r = dict(A13=a13(d, mult)); save('A13_A16_chip_results', {**res, name: r})
        log(f'{name}: A13 done {r["A13"]["predictions_met"]}')
        r['A14'] = a14(d, mult); log(f'{name}: A14 {r["A14"]["reading"]}')
        cl, mism, secs = timed_clamps(d, code, origin)
        r['A16_state_mismatches_vs_saved'] = mism; r['A16_seconds'] = secs
        np.savez_compressed(OUT / f'A16_{name}_timed_clamps.npz', **{f'{s}_{k}': v for s in cl for k, v in cl[s].items()})
        r['A16'] = a16(d, cl, mult); res[name] = r; save('A13_A16_chip_results', res)
        log(f'{name}: A16 done {r["A16"]["predictions_met"]}')
    log('done')


def selftest():
    name, cdir, code, origin = programs()[0]
    t0 = time.time(); d = build(cdir); log(f'build {time.time() - t0:.0f} s')
    cl, mism, secs = timed_clamps(d, code, origin, n_moments=2, limit_nodes=50)
    per_clamp = secs / (4 * 50)
    log(f'state reproduction check: {mism} mismatched node-states over 4 snapshots (must be 0)')
    log(f'timing: {per_clamp * 1000:.1f} ms per clamp; full timed-clamp run about {per_clamp * 80 * d["n"] / 60:.0f} min per program')
    assert mism == 0
    log('SELF-TEST PASSED')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--selftest', action='store_true'); a = ap.parse_args()
    selftest() if a.selftest else main_run()
