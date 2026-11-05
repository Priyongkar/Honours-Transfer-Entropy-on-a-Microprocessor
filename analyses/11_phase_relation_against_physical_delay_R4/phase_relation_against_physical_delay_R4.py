"""R4: physical delay (timed clamps) against phase relation (locking) on the simulated chip
(design: PREDICTIONS_R4.md). Usage: python run_delay_vs_phase.py"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / 'src'))  # repository library
from te6502.paths import ANALYSES, DATA_ROOT, RESULTS, HELPERS_FILE, load_analysis  # noqa: E402,F401
import hashlib, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
A16_FILE = ANALYSES / '05_consequences_of_the_account_A12_to_A19/chip_coverage_remainder_and_timed_clamps_A13_A14_A16.py'
A16 = load_analysis(A16_FILE.relative_to(ANALYSES).as_posix())
from te6502.test_helpers import tk, log, boot_mult, auc_with_ci, diff_ci, quintile, CHIP_DIR    # noqa: E402
SAVED = DATA_ROOT / 'analyses/new_tests/output'   # timed clamps and the program-2 run (large; not in the repository)
RES = RESULTS / '11_delay_vs_phase_R4'; RES.mkdir(exist_ok=True)

PRED = HERE / 'design/predictions_R4.md'
SEED, DRAWS = 2026, 1000
MIN_POS = 10


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def classes(d, name):
    z = np.load(SAVED / f'A16_{name}_timed_clamps.npz')
    U, n = d['U'], d['n']; off = ~np.eye(n, dtype=bool)
    I0 = {s: (z[f'{s}_E0'][:, U] > 0) & off for s in ('high', 'low')}
    I1 = {s: (z[f'{s}_E1'][:, U] > 0) & off for s in ('high', 'low')}
    imm = I0['high'] | I0['low']
    one = (I1['high'] | I1['low']) & ~imm
    none = ~I0['high'] & ~I0['low'] & ~I1['high'] & ~I1['low']; none = none & none.T & off
    ph = d['ph']; lk = ph < 2
    rel = {'same': lk[:, None] & lk[None, :] & (ph[:, None] == ph[None, :]),
           'opposite': lk[:, None] & lk[None, :] & (ph[:, None] != ph[None, :]),
           'mixed': ~(lk[:, None] & lk[None, :])}
    return imm, one, none, rel


def analyse(d, name):
    imm, one, none, rel = classes(d, name)
    ug, ks, Xg, rq = d['ug'], d['ks'], d['Xg'], d['rq']
    S = {u: tk.te_matrix(Xg, ks, delay=u)[0] for u in (1, 0, 2)}
    mult = boot_mult(d['n'], DRAWS, SEED); valid = d['valid']
    out, bs = {}, {}
    for r, R in rel.items():
        for cname, P in (('immediate', imm), ('one_step', one)):
            ii, jj = np.nonzero(valid & R & (P | none))
            pos = P[ii, jj]; strata = rq[ii] * 5 + rq[jj]
            cell = dict(positives=int(pos.sum()), negatives=int((~pos).sum()))
            if pos.sum() >= 1 and (~pos).sum() >= 1:
                for u, lab in ((1, 'lag1_te'), (0, 'lag0_cmi'), (2, 'lag2_te')):
                    rr, b = auc_with_ci(np.nan_to_num(S[u][ug[ii], ug[jj]], nan=-np.inf), pos, strata, jj, mult)
                    cell[lab] = rr; bs[(r, cname, lab)] = b
            out[f'{r}_{cname}'] = cell
            log(f"{name}: {r}/{cname}: {cell['positives']} pos, {cell['negatives']} neg; "
                f"lag1 {cell.get('lag1_te', {}).get('auc', float('nan')):.3f}")
    def ok(r, c): return out[f'{r}_{c}']['positives'] >= MIN_POS and (r, c, 'lag1_te') in bs
    pm = {}
    if ok('opposite', 'immediate'):
        a = out['opposite_immediate']
        pm['R4.1'] = bool(a['lag1_te']['lo'] > 0.5)
        out['R4.2_lag1_minus_lag0'] = diff_ci(bs[('opposite', 'immediate', 'lag1_te')], bs[('opposite', 'immediate', 'lag0_cmi')],
                                             a['lag1_te']['auc'], a['lag0_cmi']['auc'])
        pm['R4.2'] = bool(out['R4.2_lag1_minus_lag0']['lo'] > 0)
    if ok('opposite', 'immediate') and ok('same', 'one_step'):
        out['R4.3_opp_imm_minus_same_one'] = diff_ci(bs[('opposite', 'immediate', 'lag1_te')], bs[('same', 'one_step', 'lag1_te')],
                                                    out['opposite_immediate']['lag1_te']['auc'], out['same_one_step']['lag1_te']['auc'])
        pm['R4.3'] = bool(out['R4.3_opp_imm_minus_same_one']['lo'] > 0)
    if ok('mixed', 'one_step') and ok('mixed', 'immediate'):
        out['R4.4_mixed_one_minus_imm'] = diff_ci(bs[('mixed', 'one_step', 'lag1_te')], bs[('mixed', 'immediate', 'lag1_te')],
                                                 out['mixed_one_step']['lag1_te']['auc'], out['mixed_immediate']['lag1_te']['auc'])
        pm['R4.4'] = bool(out['R4.4_mixed_one_minus_imm']['lo'] > 0)
    out['predictions_met'] = pm
    return out


def main():
    hashes = {p.name: sha(p) for p in (PRED, Path(__file__).resolve(), HELPERS_FILE, A16_FILE)}
    for k, v in hashes.items():
        log(f'R4: {k} sha256 {v}')
    res = dict(hashes=hashes, draws=DRAWS)
    for name, cdir in (('program1', CHIP_DIR), ('program2', SAVED / 'A7_chip')):  # as A16.programs(), saved runs only
        d = A16.build(cdir)
        res[name] = analyse(d, name)
        log(f"{name}: {res[name]['predictions_met']}")
    def conv(v):
        if isinstance(v, dict): return {str(k): conv(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)): return [conv(x) for x in v]
        if isinstance(v, np.integer): return int(v)
        if isinstance(v, np.floating): return float(v)
        if isinstance(v, np.bool_): return bool(v)
        return v
    (RES / 'R4_delay_vs_phase_results.json').write_text(json.dumps(conv(res), indent=1))
    log('done')


if __name__ == '__main__':
    main()
