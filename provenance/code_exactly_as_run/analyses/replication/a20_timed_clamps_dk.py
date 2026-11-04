"""A20: timed clamps on the recorded Donkey Kong run (design: PREDICTIONS_A20.md).

Usage: python a20_timed_clamps_dk.py [--selftest]
  --selftest: checks snapshot/restore determinism and that a forced wire holds; no A20 statistic.
"""
import argparse, builtins, copy, hashlib, json, sys, time
from array import array
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'vendor/Sim2600'))
PRED = HERE / 'PREDICTIONS_A20.md'
ADD = HERE / 'PREDICTIONS_A20_ADDENDUM.md'
OUT = HERE / 'output'; OUT.mkdir(exist_ok=True)
RUN = ROOT / 'results/v5_bitcorrect_every_third_mi_ais_te_k1_12_L1_d1_ais_x1'
VSS, VCC = 558, 657
_print = builtins.print


def quiet(f, *a):
    builtins.print = lambda *x, **k: None
    try: return f(*a)
    finally: builtins.print = _print


def log(m):
    _print(f'[{time.strftime("%H:%M:%S")}] {m}', flush=True)


def console():
    from sim2600 import params, sim6502, simTIA
    from sim2600.sim2600Console import Sim2600Console
    return quiet(Sim2600Console, params.ROMS_DONKEY_KONG, sim6502.MySim6502, simTIA.MySimTIA)


ARR = ('_wireState', '_wirePulled', '_transistorState', '_forced')


def snapshot(c):
    s = {}
    for name, obj in (('cpu', c.sim6507), ('tia', c.simTIA)):
        s[name] = {a: getattr(obj, a).copy() for a in ARR}
        s[name + '_scalars'] = {k: v for k, v in obj.__dict__.items() if isinstance(v, (int, float, bool))}
    s['pia'] = copy.deepcopy(c.emuPIA.__dict__)
    s['con'] = {k: copy.deepcopy(v) for k, v in c.__dict__.items() if k not in ('sim6507', 'simTIA', 'emuPIA')}
    return s


def restore(c, s):
    for name, obj in (('cpu', c.sim6507), ('tia', c.simTIA)):
        for a in ARR:
            getattr(obj, a)[:] = s[name][a]
        obj.__dict__.update(s[name + '_scalars'])
    c.emuPIA.__dict__.update(copy.deepcopy(s['pia']))
    for k, v in s['con'].items():
        setattr(c, k, copy.deepcopy(v))


def step(c, n=1):
    for _ in range(n):
        quiet(c.advanceOneHalfClock)


def latch_transistors(cpu):
    p2g = np.load(RUN / 'groups.npz')['physical_to_group'].astype(int)
    clock_g = 6
    fets = cpu._transistorList
    gate = np.array([f.gateWireIndex for f in fets]); c1 = np.array([f.side1WireIndex for f in fets])
    c2 = np.array([f.side2WireIndex for f in fets])
    pullup = {i for i, w in enumerate(cpu._wirePulled) if w == 1}            # PULLED_HIGH at build
    pd_nodes = {int(c1[t]) for t in range(len(gate)) if c2[t] == VSS} | {int(c2[t]) for t in range(len(gate)) if c1[t] == VSS}
    driven = lambda nd: nd in pullup or nd in pd_nodes
    gated_by = {}
    for t, g in enumerate(gate):
        gated_by.setdefault(int(g), []).append(t)
    clock_gates = {int(gate[t]) for t in np.flatnonzero(p2g == clock_g)}
    groups_of = lambda nd: {int(p2g[t]) for t in gated_by.get(int(nd), []) if p2g[t] >= 0}
    lat = []
    for t in range(len(gate)):
        if int(gate[t]) in clock_gates and c1[t] not in (VSS, VCC) and c2[t] not in (VSS, VCC):
            a1, a2 = int(c1[t]), int(c2[t]); d1, d2 = driven(a1), driven(a2)
            if d1 == d2: continue
            a, b = (a1, a2) if d1 else (a2, a1)
            pairs = sorted({(ga, gb) for ga in groups_of(a) for gb in groups_of(b) if ga != gb})
            if pairs:
                lat.append(dict(t=t, driven=a, storage=b, pairs=pairs))
    return lat


def high(obj, w):
    return bool(obj._wireState[w] & (1 | 8 | 16))


def clamp_once(c, s, a, b):
    """E0, E1 for forcing driven node a to the opposite of its value at snapshot s."""
    cpu = c.sim6507
    restore(c, s); b0 = high(cpu, b); va = high(cpu, a)
    step(c, 3); b3 = high(cpu, b)
    restore(c, s)
    cpu._forced[a] = 1 if va else 2
    quiet(cpu.recalcWire, a)
    e0 = high(cpu, b) != b0
    step(c, 3); e1 = high(cpu, b) != b3
    cpu._forced[a] = 0
    return e0, e1


def selftest():
    c = console(); step(c, 5000); s = snapshot(c)
    step(c, 3); ref = c.sim6507.getTransistorState().copy()
    restore(c, s); step(c, 3); again = c.sim6507.getTransistorState().copy()
    det = bool((ref == again).all())
    restore(c, s); cpu = c.sim6507
    lat = latch_transistors(cpu)
    a = lat[0]['driven']; va = high(cpu, a)
    cpu._forced[a] = 1 if va else 2; quiet(cpu.recalcWire, a); held = high(cpu, a) != va
    step(c, 3); held3 = high(cpu, a) != va; cpu._forced[a] = 0
    pairs = sorted({p for l in lat for p in l['pairs']})
    t0 = time.time(); clamp_once(c, s, lat[0]['driven'], lat[0]['storage']); dt = time.time() - t0
    out = dict(restore_deterministic=det, force_holds_now=held, force_holds_after_3=held3,
               latch_transistors=len(lat), latch_group_pairs=len(pairs), seconds_per_clamp=round(dt, 3))
    log(json.dumps(out)); assert det and held and held3
    log('SELF-TEST PASSED')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--selftest', action='store_true'); a = ap.parse_args()
    if a.selftest:
        return selftest()
    hashes = {p.name: hashlib.sha256(Path(p).read_bytes()).hexdigest()
              for p in (PRED, ADD, Path(__file__).resolve(), ROOT / 'vendor/Sim2600/sim2600/mycircuitsimulator.pyx')}
    for k, v in hashes.items():
        log(f'A20: {k} sha256 {v}')
    c = console(); cpu = c.sim6507
    lat = latch_transistors(cpu)
    rows = np.arange(4224, 401776, 3)
    moments = np.sort(np.random.default_rng(2026).choice(rows, 80, replace=False))
    log(f'{len(lat)} latch transistors; moments {moments[:5].tolist()} ... {moments[-3:].tolist()}')
    rec = {l['t']: dict(closed=[], open=[]) for l in lat}
    nonconv = []
    # row r of the recording = state after half-clock r+1 differs from after r; the moment is after half-clock r+1
    done = 0; t0 = time.time()
    step(c, 1)                                            # half-clock 1 (row 0's reference state)
    pos = 1
    for r in moments:
        step(c, int(r + 1 - pos)); pos = int(r + 1)
        s = snapshot(c)
        for l in lat:
            open_ = bool(s['cpu']['_transistorState'][l['t']])
            try:
                e0, e1 = clamp_once(c, s, l['driven'], l['storage'])
                rec[l['t']]['open' if open_ else 'closed'].append((int(e0), int(e1)))
            except RuntimeError:
                nonconv.append(dict(t=l['t'], row=int(r), open=open_))
                c = console(); cpu = c.sim6507
                restore(c, s)
        restore(c, s); done += 1
        if done % 10 == 0:
            log(f'{done}/80 moments, {time.time() - t0:.0f} s')
    per_t = {}
    for l in lat:
        cl, op = rec[l['t']]['closed'], rec[l['t']]['open']
        if len(cl) < 5:
            per_t[l['t']] = dict(driven=l['driven'], storage=l['storage'], pairs=l['pairs'], n_closed=len(cl), n_open=len(op),
                                 delayed_only_closed=False, immediate_open=any(e0 == 1 for e0, _ in op), unclassified=True)
            continue
        d_only = bool(cl) and all(e0 == 0 for e0, _ in cl) and any(e1 == 1 for _, e1 in cl)
        imm = any(e0 == 1 for e0, _ in op)
        per_t[l['t']] = dict(driven=l['driven'], storage=l['storage'], pairs=l['pairs'], n_closed=len(cl), n_open=len(op),
                             delayed_only_closed=d_only, immediate_open=imm,
                             share_closed_e0_0_e1_1=float(np.mean([e0 == 0 and e1 == 1 for e0, e1 in cl])) if cl else None,
                             share_open_e0_1=float(np.mean([e0 == 1 for e0, _ in op])) if op else None)
    pairs = sorted({tuple(p) for l in lat for p in l['pairs']})
    pd = {p: [v for v in per_t.values() if tuple(p) in [tuple(q) for q in v['pairs']]] for p in pairs}
    n1 = sum(any(v['delayed_only_closed'] for v in vs) for vs in pd.values())
    n2 = sum(any(v['immediate_open'] for v in vs) for vs in pd.values())
    n3 = sum(any(v['delayed_only_closed'] and v['immediate_open'] for v in vs) for vs in pd.values())
    res = dict(hashes=hashes, moments=moments.tolist(), latch_transistors=len(lat), latch_group_pairs=len(pairs),
               A20_1_delayed_only_closed_pairs=n1, A20_2_immediate_open_pairs=n2, A20_3_confirmed_pairs=n3,
               predictions_met={'A20.1': n1 > 20, 'A20.2': n2 > 20, 'A20.3': n3 > 20},
               per_transistor={str(k): v for k, v in per_t.items()}, nonconvergent=nonconv, seconds=round(time.time() - t0, 1))
    (OUT / 'A20_results.json').write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o)))
    log(f'non-convergent clamps: {len(nonconv)}')
    log(f"done: pairs {len(pairs)}; delayed-only closed {n1}; immediate open {n2}; confirmed {n3}; {res['predictions_met']}")


if __name__ == '__main__':
    main()
