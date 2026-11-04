"""A17: the stratified H6 design (IDTxl, source lags 0 and 1) on the clock-aligned samples (1,408 onwards),
R1's 100 targets, every other setting as R1 (design: PREDICTIONS_A12.md). Reuses idtxl_runs.py unchanged.

Usage:
  python a17_idtxl_aligned.py --selftest   builds the data slice and the job list; runs no IDTxl
  python a17_idtxl_aligned.py [--workers 3]
"""
import argparse, hashlib, json, os
from pathlib import Path
import numpy as np
from common import HERE, RUN, OUT, tk, log, save, boot_mult
from idtxl_runs import BASE, DATA, run_jobs, summaries, network_stats

PRED = HERE / 'PREDICTIONS_A12.md'
ALIGNED = OUT / 'data_clock_aligned.npy'
FIRST = 1407            # index of sample 1,408


def edges_per_target(summ, G, targets, ref):
    for s in summ:
        s['selected_vars_sources'] = [v for v in s['selected_vars_sources'] if v[0] < G]
    keep = tk.network_fdr_omnibus(summ)
    tix = {t: i for i, t in enumerate(targets)}; h = np.zeros(len(targets)); c = np.zeros(len(targets))
    for s in summ:
        if keep.get(s['target']):
            for u, lag in s['selected_vars_sources']:
                h[tix[s['target']]] += ref[u, s['target']]; c[tix[s['target']]] += 1
    return h, c


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--selftest', action='store_true'); ap.add_argument('--workers', type=int, default=3)
    a = ap.parse_args()
    os.environ.setdefault('JAVA_HOME', os.path.expanduser('~/opt/anaconda3/envs/honours6502/lib/jvm'))
    targets = json.load(open(OUT / 'idtxl_results.json'))['targets']
    if not ALIGNED.exists():
        X = np.load(DATA); np.save(ALIGNED, np.ascontiguousarray(X[:, FIRST:]))
    shape = np.load(ALIGNED, mmap_mode='r').shape
    out_dir = OUT / 'A17_targets'; out_dir.mkdir(exist_ok=True)
    jobs = [dict(data_path=str(ALIGNED), target=int(t), method='multivariate', max_heap='2g',
                 settings=dict(BASE, min_lag_sources=0), out_path=str(out_dir / f'target_{t:04d}.json')) for t in targets]
    if a.selftest:
        log(f'aligned data {shape}; {len(jobs)} jobs; settings {jobs[0]["settings"]}; JAVA_HOME set: {bool(os.environ.get("JAVA_HOME"))}')
        log('SELF-TEST PASSED'); return
    hashes = {p.name: hashlib.sha256(Path(p).read_bytes()).hexdigest()
              for p in (PRED, Path(__file__).resolve(), HERE / 'common.py', HERE / 'idtxl_runs.py')}
    for k, v in hashes.items():
        log(f'A17: {k} sha256 {v}')
    log(f'aligned data {shape}')
    run_jobs(tk.idtxl_target_worker, jobs, a.workers)
    ref = np.load(OUT / 'recording_reference.npz'); wired, latch, phase = ref['wired'], ref['latch'], ref['phase']
    net_main = np.load(RUN / 'thesis_v5/idtxl/network_multivariate.npz')['adjacency'].astype(bool)
    mult = boot_mult(len(targets), 1000)
    stats = network_stats(summaries(out_dir), wired, latch, phase, net_main, list(targets), mult)
    G = len(wired)
    h17, c17 = edges_per_target(summaries(out_dir), G, targets, wired)
    hR1, cR1 = edges_per_target(summaries(OUT / 'R1_targets'), G, targets, wired)
    d = [(m @ h17) / max(m @ c17, 1) - (m @ hR1) / max(m @ cR1, 1) for m in mult]
    diff_r1 = dict(diff=h17.sum() / max(c17.sum(), 1) - hR1.sum() / max(cR1.sum(), 1),
                   lo=float(np.nanquantile(d, .025)), hi=float(np.nanquantile(d, .975)))
    res = dict(hashes=hashes, targets=targets, A17=stats, A17_minus_R1=diff_r1,
               predictions_met={'A17.1': bool(stats['diff_vs_main_netlist']['lo'] > 0),
                                'A17.2': bool(diff_r1['lo'] <= 0 <= diff_r1['hi'])})
    save('A17_results', res); log(f"done: {res['predictions_met']}")


if __name__ == '__main__':
    main()
