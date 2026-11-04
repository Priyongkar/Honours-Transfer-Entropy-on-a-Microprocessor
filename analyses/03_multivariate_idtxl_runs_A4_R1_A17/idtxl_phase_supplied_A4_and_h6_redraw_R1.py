"""R1 (stratified H6 design, re-drawn) and A4 (IDTxl with the phase forced into the
conditioning set) on the reduced recording (design: PREDICTIONS.md).
Usage: python idtxl_runs.py [--quick] [--workers 3]"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / 'src'))  # repository library
from te6502.paths import ANALYSES, DATA_ROOT, RESULTS, HELPERS_FILE, load_analysis  # noqa: E402,F401
import argparse, glob, json, os, time
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from te6502.test_helpers import RUN, OUT, tk, log, start, save, boot_mult

DATA = RUN / 'thesis_v5/idtxl/data_int8.npy'
BASE = dict(max_lag_target=12, max_lag_sources=1, alpha=0.05, n_perm=200, analytic_bits_fix=True)


def draw_targets(K, rate, n=100, seed=2026):
    active = np.flatnonzero((rate > 0) & (rate < 1))
    bands = [(1, 1), (2, 4), (5, 8), (9, 11), (12, 12)]
    members = [active[(K[active] >= lo) & (K[active] <= hi)] for lo, hi in bands]
    sizes = np.asarray([len(m) for m in members], float)
    quota = sizes / sizes.sum() * n; alloc = np.floor(quota).astype(int)
    for i in np.argsort(-(quota - alloc))[:n - alloc.sum()]:
        alloc[i] += 1
    rng = np.random.default_rng(seed)
    picks = [rng.choice(m, size=a, replace=False) for m, a in zip(members, alloc) if a]
    return np.sort(np.concatenate(picks)).astype(int), dict(zip([f'{lo}-{hi}' for lo, hi in bands], alloc.tolist()))


def a4_worker(job):
    """As tk.idtxl_target_worker, plus forced conditionals."""
    tk.ensure_jvm(max_heap=job['max_heap']); tk.patch_idtxl()
    from idtxl.data import Data
    from idtxl.multivariate_te import MultivariateTE
    values = np.load(job['data_path'], mmap_mode='r')
    data = tk._quiet(Data, np.asarray(values, dtype=np.int32), dim_order='ps', normalise=False)
    settings = tk.idtxl_settings(**job['settings']); settings['add_conditionals'] = job['add_conditionals']
    t0 = time.perf_counter()
    res = tk._quiet(MultivariateTE().analyse_single_target, settings=settings, data=data,
                    target=int(job['target']), sources=job['sources'])
    s = tk._single_summary(res, job['target'], 'multivariate', time.perf_counter() - t0,
                           n_candidates=len(job['sources']), n_obs=values.shape[1] - 12)
    with open(job['out_path'] + '.tmp', 'w') as fh:
        json.dump(s, fh)
    os.replace(job['out_path'] + '.tmp', job['out_path'])
    return s


def run_jobs(fn, jobs, workers):
    todo = [j for j in jobs if not os.path.exists(j['out_path'])]
    log(f'{len(jobs) - len(todo)} cached, {len(todo)} to run on {workers} workers')
    with ProcessPoolExecutor(workers, mp_context=mp.get_context('spawn')) as ex:
        futs = {ex.submit(fn, j): j['target'] for j in todo}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                f.result()
            except Exception as e:
                log(f'target {futs[f]} failed: {e!r}')
            if i % 10 == 0 or i == len(todo):
                log(f'  {i}/{len(todo)} done')


def summaries(folder):
    return [json.load(open(f)) for f in sorted(glob.glob(str(folder / 'target_*.json')))]


def network_stats(summ, wired, latch, phase, net_main, targets, mult):
    G = len(wired)   # the forced parity conditional (index G) is not an edge
    for s in summ:
        s['selected_vars_sources'] = [v for v in s['selected_vars_sources'] if v[0] < G]
    keep = tk.network_fdr_omnibus(summ)
    E = []
    for s in summ:
        if keep.get(s['target']):
            E += [(u, s['target'], lag) for u, lag in s['selected_vars_sources']]
    def prec(edges, ref):
        return float(np.mean([ref[u, t] for u, t, _ in edges])) if edges else np.nan
    by_lag = {f'lag{l}': dict(edges=sum(1 for e in E if e[2] == l), wired=sum(1 for e in E if e[2] == l and wired[e[0], e[1]]),
                               precision_netlist=prec([e for e in E if e[2] == l], wired)) for l in (0, 1)}
    main_edges = [(u, t, 1) for t in targets for u in np.flatnonzero(net_main[:, t])]
    tix = {t: i for i, t in enumerate(targets)}
    def per_target(edges, ref):
        h = np.zeros(len(targets)); c = np.zeros(len(targets))
        for u, t, _ in edges:
            h[tix[t]] += ref[u, t]; c[tix[t]] += 1
        return h, c
    out = dict(targets=len(summ), targets_kept=int(sum(keep.values())), edges=len(E),
               wired=int(sum(wired[u, t] for u, t, _ in E)), precision_netlist=prec(E, wired),
               precision_latch=prec(E, latch), by_lag=by_lag,
               mean_parent_set=float(np.mean([len(s['selected_vars_sources']) for s in summ])),
               share_opposite_phase_sources=float(np.mean([(phase[u] < 2) and (phase[t] < 2) and phase[u] != phase[t] for u, t, _ in E])) if E else np.nan,
               main_run_same_targets=dict(edges=len(main_edges), precision_netlist=prec(main_edges, wired), precision_latch=prec(main_edges, latch)))
    for name, ref in [('netlist', wired), ('latch', latch)]:
        h1, c1 = per_target(E, ref); h0, c0 = per_target(main_edges, ref)
        d = [(m @ h1) / max(m @ c1, 1) - (m @ h0) / max(m @ c0, 1) for m in mult]
        out[f'diff_vs_main_{name}'] = dict(diff=out[f'precision_{name}'] - out['main_run_same_targets'][f'precision_{name}'],
                                           lo=float(np.nanquantile(d, .025)), hi=float(np.nanquantile(d, .975)))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); ap.add_argument('--workers', type=int, default=3)
    a = ap.parse_args()
    sha = start('idtxl_runs')
    os.environ.setdefault('JAVA_HOME', os.path.expanduser('~/opt/anaconda3/envs/honours6502/lib/jvm'))
    ref = np.load(OUT / 'recording_reference.npz')
    wired, latch, phase, rate, K = ref['wired'], ref['latch'], ref['phase'], ref['rate'], ref['k']
    targets, alloc = draw_targets(K, rate)
    if a.quick:
        targets = targets[[0, len(targets) // 2]]
    log(f'targets {len(targets)}; allocation {alloc}')
    net_main = np.load(RUN / 'thesis_v5/idtxl/network_multivariate.npz')['adjacency'].astype(bool)
    G = len(K); sub = 'quick_' if a.quick else ''
    r1_dir = OUT / f'{sub}R1_targets'; a4_dir = OUT / f'{sub}A4_targets'; r1_dir.mkdir(exist_ok=True); a4_dir.mkdir(exist_ok=True)
    # R1
    jobs = [dict(data_path=str(DATA), target=int(t), method='multivariate', max_heap='2g',
                 settings=dict(BASE, min_lag_sources=0), out_path=str(r1_dir / f'target_{t:04d}.json')) for t in targets]
    run_jobs(tk.idtxl_target_worker, jobs, a.workers)
    # A4
    pdata = OUT / 'data_with_parity.npy'
    if not pdata.exists():
        X = np.load(DATA); P = (np.arange(X.shape[1]) % 2).astype(np.int8)
        np.save(pdata, np.vstack([X, P[None, :]]))
    jobs = [dict(data_path=str(pdata), target=int(t), max_heap='2g', settings=dict(BASE, min_lag_sources=1),
                 add_conditionals=[(G, 1)], sources=[u for u in range(G) if u != t],
                 out_path=str(a4_dir / f'target_{t:04d}.json')) for t in targets]
    run_jobs(a4_worker, jobs, a.workers)
    mult = boot_mult(len(targets), 50 if a.quick else 1000)
    res = dict(predictions_sha256=sha, quick=a.quick, targets=targets, allocation=alloc,
               thesis_H6=dict(edges=236, wired=41, precision_percent=17.37),
               R1=network_stats(summaries(r1_dir), wired, latch, phase, net_main, list(targets), mult),
               A4=network_stats(summaries(a4_dir), wired, latch, phase, net_main, list(targets), mult))
    save('idtxl_results' + ('_quick' if a.quick else ''), res)
    log('done')


if __name__ == '__main__':
    main()
