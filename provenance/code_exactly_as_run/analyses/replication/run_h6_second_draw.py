"""R5: a second, independent stratified draw of the H6 design (design: PREDICTIONS_R5.md).
Usage: python run_h6_second_draw.py [--workers 3]   (run with the honours6502 environment)"""
import argparse, hashlib, json, os, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
NT = HERE.parent / 'new_tests'
sys.path.insert(0, str(NT))
import idtxl_runs as IR                                        # noqa: E402
from common import RUN, OUT, tk, log                             # noqa: E402

PRED = HERE / 'PREDICTIONS_R5.md'
SEED = 2027


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--workers', type=int, default=3); a = ap.parse_args()
    for p in (PRED, Path(__file__).resolve(), NT / 'idtxl_runs.py', NT / 'common.py'):
        log(f'R5: {p.name} sha256 {sha(p)}')
    os.environ.setdefault('JAVA_HOME', os.path.expanduser('~/opt/anaconda3/envs/honours6502/lib/jvm'))
    ref = np.load(OUT / 'recording_reference.npz'); K, rate = ref['k'], ref['rate']
    r1 = set(json.load(open(OUT / 'idtxl_results.json'))['targets'])
    rate_ex = rate.copy(); rate_ex[list(r1)] = 0.0            # excluded targets look inactive to the drawer
    targets, alloc = IR.draw_targets(K, rate_ex, n=100, seed=SEED)
    assert not (set(targets.tolist()) & r1) and len(targets) == 100
    log(f'R5 targets {len(targets)} (none in R1); allocation {alloc}')
    out_dir = OUT / 'R5_targets'; out_dir.mkdir(exist_ok=True)
    (HERE / 'output/R5_targets.json').write_text(json.dumps(dict(seed=SEED, targets=targets.tolist(), allocation=alloc)))
    jobs = [dict(data_path=str(IR.DATA), target=int(t), method='multivariate', max_heap='2g',
                 settings=dict(IR.BASE, min_lag_sources=0), out_path=str(out_dir / f'target_{t:04d}.json')) for t in targets]
    IR.run_jobs(tk.idtxl_target_worker, jobs, a.workers)
    log('R5 search done')


if __name__ == '__main__':
    main()
