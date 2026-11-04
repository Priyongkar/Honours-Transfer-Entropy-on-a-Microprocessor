"""Supplementary (not hash-logged; asked for by examiners after the runs):
(1) latch pairs against harder negatives on the simulated chip, both programs:
    immediate-influence pairs, wired pairs, and all other pairs;
(2) A8 paired difference lag-1 TE minus lag-0 CMI, and the count of latches that
    can be scored. Writes output/supp_harder.json.
Usage: python supp_harder.py"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / 'src'))  # repository library
from te6502.paths import ANALYSES, DATA_ROOT, RESULTS, HELPERS_FILE, load_analysis  # noqa: E402,F401
from pathlib import Path
import numpy as np
from te6502.test_helpers import tk, log, start, save, change, phase_class, auc_with_ci, diff_ci, boot_mult, quintile, CHIP_DIR, OUT, ROOT
from te6502.switch_level_simulator import Chip
nmos = load_analysis('02_latch_transfer_tests_A1_to_A8/nmos_circuit_latch_tests_A8.py')


def chip_block(cdir):
    with np.load(cdir / 'states.npz') as z:
        states, nodes, clk0 = z['states'], z['nodes'], z['clk0']
    with np.load(cdir / 'interventions.npz') as z:
        E0, E1, U = z['E0'], z['E1'], z['nodes']
    T, Nn = states.shape; n = len(U)
    chip = Chip(ROOT / 'data/visual6502')
    g_, c1_, c2_ = chip.transistor_matrix(); gnd, pwr = chip.gnd, chip.pwr
    col = {int(v): i for i, v in enumerate(nodes)}
    A = tk.node_reference_from_transistors(g_, c1_, c2_, nodes, supply={gnd, pwr}).toarray() > 0
    W = A[np.ix_(U, U)]
    off = ~np.eye(n, dtype=bool)
    i0 = E0[:, U] > 0; np.fill_diagonal(i0, False)
    I = i0 & off
    pos_in_U = {int(u): k for k, u in enumerate(U)}
    st_all = states.T
    clock_cols = {c for c in range(Nn) if (st_all[c] == clk0).all() or (st_all[c] == 1 - clk0).all()}
    pd_nodes = {int(c1_[t]) for t in range(len(g_)) if c2_[t] == gnd}
    driven = lambda nd: bool(chip.pullup[nd]) or nd in pd_nodes
    L = np.zeros((n, n), bool)
    for t in range(len(g_)):
        g, c1, c2 = int(g_[t]), int(c1_[t]), int(c2_[t])
        if col.get(g) in clock_cols and c1 not in (gnd, pwr) and c2 not in (gnd, pwr):
            if c1 in col and c2 in col and col[c1] in pos_in_U and col[c2] in pos_in_U:
                a1, a2 = pos_in_U[col[c1]], pos_in_U[col[c2]]
                d1, d2 = driven(c1), driven(c2)
                if d1 != d2:
                    L[(a1, a2) if d1 else (a2, a1)] = True
    chg = change(st_all[U].astype(np.uint8))
    ph, _ = phase_class(chg); rq = quintile(chg.mean(1))
    reps, ug = tk.group_identical(chg); Xg = chg[reps]
    ks = np.asarray([tk.ais_select_k(x, 6) for x in Xg])
    S1 = tk.te_matrix(Xg, ks, delay=1)[0]
    valid = (ug[:, None] >= 0) & (ug[None, :] >= 0) & (ug[:, None] != ug[None, :]) & off
    ii, jj = np.nonzero(valid)
    base = ((ph[ii] * 3 + ph[jj]) * 5 + rq[ii]) * 5 + rq[jj]
    s = S1[ug[ii], ug[jj]]
    cL = L[ii, jj]; cI = I[ii, jj] & ~cL; cW = (W | W.T)[ii, jj] & ~cL
    mult = boot_mult(n, 1000)
    out = dict(latch_scored=int(cL.sum()), latch_in_immediate=int((I[ii, jj] & cL).sum()))
    for name, neg in [('vs_immediate', cI), ('vs_wired', cW), ('vs_all_other', ~cL)]:
        m = cL | neg
        r, _ = auc_with_ci(s[m], cL[m], base[m], jj[m], mult)
        r['n_neg'] = int(neg.sum()); out[name] = r
        log(f'{cdir.name}: latch {name}: {r["auc"]:.3f} [{r["lo"]:.3f}, {r["hi"]:.3f}] (neg {neg.sum()})')
    return out


def nmos_block(T=6000):
    cols = {k: [] for k in ['lag1', 'lag0', 'inp', 'rdr', 'neg', 'strat', 'circ']}
    counts = dict(latches=0, latches_with_nodes=0, inp_scorable=0)
    for c, seed in enumerate(nmos.SEEDS):
        defs, meta = nmos.sim6502.make_nmos_circuit(n_gates=80, n_inputs=6, n_latches=16, seed=seed)
        counts['latches'] += len(meta['latches'])
        S, A, inp, rdr = nmos.circuit(seed, T)
        counts['latches_with_nodes'] += len(inp)
        D = change(S); reps, ug = tk.group_identical(D); Xg = D[reps]
        ks = np.asarray([tk.ais_select_k(x, 4) for x in Xg])
        L1 = tk.te_matrix(Xg, ks, delay=1)[0]; L0 = tk.te_matrix(Xg, ks, delay=0)[0]
        n = len(S); ph, _ = phase_class(D); rq = quintile(D.mean(1))
        Im = np.zeros((n, n), bool); Rm = np.zeros((n, n), bool)
        for p in inp: Im[p] = True
        for p in rdr: Rm[p] = True
        rel = A | A.T | Im | Im.T | Rm | Rm.T
        valid = (ug[:, None] >= 0) & (ug[None, :] >= 0) & (ug[:, None] != ug[None, :])
        ii, jj = np.nonzero(valid)
        counts['inp_scorable'] += int(Im[ii, jj].sum())
        cols['lag1'].append(L1[ug[ii], ug[jj]]); cols['lag0'].append(L0[ug[ii], ug[jj]])
        cols['inp'].append(Im[ii, jj]); cols['rdr'].append(Rm[ii, jj]); cols['neg'].append(~rel[ii, jj])
        cols['strat'].append((((c * 3 + ph[ii]) * 3 + ph[jj]) * 5 + rq[ii]) * 5 + rq[jj])
        cols['circ'].append(np.full(len(ii), c))
    C = {k: np.concatenate(v) for k, v in cols.items()}
    mult = np.random.default_rng(2026).multinomial(10, np.full(10, 0.1), size=1000)
    out = dict(counts=counts)
    for pname in ['inp', 'rdr']:
        m = C[pname] | C['neg']
        r1, b1 = auc_with_ci(C['lag1'][m], C[pname][m], C['strat'][m], C['circ'][m], mult)
        r0, b0 = auc_with_ci(C['lag0'][m], C[pname][m], C['strat'][m], C['circ'][m], mult)
        out[f'{pname}_lag1_minus_lag0'] = diff_ci(b1, b0, r1['auc'], r0['auc'])
        log(f"NMOS {pname}: lag1-lag0 {out[f'{pname}_lag1_minus_lag0']}")
    log(f'NMOS counts {counts}')
    return out


def main():
    start('supp_harder')
    res = dict(status='supplementary; not hash-logged',
               program1=chip_block(Path(CHIP_DIR)), program2=chip_block(OUT / 'A7_chip'), nmos_T6000=nmos_block())
    save('supp_harder', res)
    log('done')


if __name__ == '__main__':
    main()
