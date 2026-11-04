"""Supplementary check on A16 (after the run): each latch pair in its open-latch set (immediate?) and closed-latch set."""
import sys, numpy as np
sys.path.insert(0, '.')
from a13_a16_chip import build, programs, OUT
for name, cdir, code, origin in programs():
    d = build(cdir); U, n = d['U'], d['n']; off = ~np.eye(n, dtype=bool)
    z = np.load(OUT / f'A16_{name}_timed_clamps.npz')
    I0 = {s: (z[f'{s}_E0'][:, U] > 0) & off for s in ('high', 'low')}
    I1 = {s: (z[f'{s}_E1'][:, U] > 0) & off for s in ('high', 'low')}
    op = cl = 0; frac_closed = []
    for (a, b), p in d['latch'].items():
        sc = 'low' if p == 1 else 'high'; so = 'high' if sc == 'low' else 'low'
        op += bool(I0[so][a, b]); cl += bool(I1[sc][a, b] and not I0[sc][a, b])
        frac_closed.append(z[f'{sc}_E1'][a, U[b]])
    print(name, 'latch', len(d['latch']), 'delayed-only closed', cl, 'immediate open', op,
          'median E1 freq closed', float(np.median(frac_closed)),
          'pairs delayed-only in closed set overall share', float((I1['low'] & ~I0['low']).mean()), float((I1['high'] & ~I0['high']).mean()))
