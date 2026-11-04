"""Convert the Pitfall recording (Wang et al. format) into the Donkey Kong format, and verify it.

Source: ~/6502-information-flow/data/raw/Pitfall/Pitfall/window_*/sequence_step_128_rec_30.npy
  [3510 transistors, 3840 sub-steps] int16 node states per window; eight contiguous windows.
  The file's 'step' is 30 sub-steps; the 6502's clock nodes (cp1, cclk) switch every 90 sub-steps,
  so one clock phase is three steps of 30, as one phase is three rows of the Donkey Kong recording.
Target (identical to data/raw/donkeykong.5000.ts.spikes.npy):
  [rows, 3510] uint8, columns = transistors in netlist order, an entry = a change of state in either
  direction, three rows per clock phase with the events on the first row of each block, so that the
  main run's reduction (row_stride 3, phase with the most events, crop to the last nonzero row)
  returns one sample per phase.

Construction: the chip is silent at sub-steps 48 to 60 of every 90-sub-step period, the last quiet
interval before the clock edge (profile checked below), so the settled state of phase h is the state
at sub-step 60 + 90h; the event of phase h is the XOR of consecutive settled states, a change of
state in either direction during the phase. Change coding is therefore constructed from states
here, not inferred. A transistor that changes twice within one phase shows no event; that count is
reported.

Writes data/raw/pitfall.ts.spikes.npy, data/raw/pitfall_states_by_phase.npy and
analyses/pitfall/output/conversion_checks.json. No transfer-entropy statistic is computed.
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / 'src'))  # repository library
from te6502.paths import ANALYSES, DATA_ROOT, RESULTS, HELPERS_FILE, load_analysis  # noqa: E402,F401
import json, hashlib, sys
from pathlib import Path
import numpy as np

ROOT = DATA_ROOT
from te6502.transfer_entropy_estimators import group_identical          # same grouping as the main run
from te6502.switch_level_simulator import load_netlist

SRC = Path.home() / '6502-information-flow/data/raw/Pitfall/Pitfall'
WINDOWS = [f'window_{a}_{a + 128}' for a in range(0, 1024, 128)]
SUB = 30
PHASE, OFFSET = 90, 60
OUT = RESULTS / '06_regenerated_recordings'; OUT.mkdir(parents=True, exist_ok=True)


def main():
    X = [np.asarray(np.load(SRC / w / 'sequence_step_128_rec_30.npy')).astype(np.uint8) for w in WINDOWS]
    checks = {'source_windows': WINDOWS,
              'source_sha256': {w: hashlib.sha256((SRC / w / 'sequence_step_128_rec_30.npy').read_bytes()).hexdigest() for w in WINDOWS},
              'join_differences': [int((X[i][:, -1] != X[i + 1][:, 0]).sum()) for i in range(len(X) - 1)]}
    full = np.concatenate(X, axis=1)                                   # [3510, 30720]
    ch = np.diff(full.astype(np.int8), axis=1) != 0
    
    prof = np.bincount(np.arange(1, full.shape[1]) % PHASE, weights=ch.sum(0), minlength=PHASE)
    checks['changes_by_substep_mod_90'] = prof.astype(int).tolist()
    checks['quiet_before_offset'] = bool(prof[48:OFFSET + 1].sum() == 0)
    S = full[:, OFFSET::PHASE]                                         # settled state per phase [3510, 341]
    blocks = ch[:, OFFSET:OFFSET + PHASE * (S.shape[1] - 1)].reshape(ch.shape[0], S.shape[1] - 1, PHASE)
    nchg = blocks.sum(2)
    checks['transistor_phases_with_even_changes_>0'] = int(((nchg > 0) & (nchg % 2 == 0)).sum())
    checks['transistor_phases_with_any_change'] = int((nchg > 0).sum())
    E = (S[:, 1:] != S[:, :-1]).astype(np.uint8)                       # change per phase [3510, 340]
    H = E.shape[1]
    raw = np.zeros((3 * H, E.shape[0]), dtype=np.uint8)
    raw[0::3] = E.T                                                    # DK layout: rows x transistors
    last = int(np.flatnonzero(raw.any(1))[-1])
    raw = raw[:last + 1]
    np.save(ROOT / 'data/raw/pitfall.ts.spikes.npy', raw)
    np.save(ROOT / 'data/raw/pitfall_states_by_phase.npy', S[:, 1:].T.astype(np.uint8))

    # 1. the main run's reduction recovers E exactly
    phase = int(np.argmax([raw[o::3].sum() for o in range(3)]))
    red = raw[phase::3].T
    checks['reduction'] = dict(raw_shape=list(raw.shape), phase_chosen=phase, samples=int(red.shape[1]),
                               equals_constructed_events=bool(np.array_equal(red, E[:, :red.shape[1]])))

    # 2. column mapping: transistors sharing a gate node fall in one group
    pullup, tdefs, names = load_netlist(ROOT / 'data/visual6502')
    gate = np.array([t[1] for t in tdefs]); assert len(gate) == red.shape[0] == 3510
    reps, ug = group_identical(red)
    act = ug >= 0
    same_gate = split = 0
    for g in np.unique(gate[act]):
        m = np.flatnonzero((gate == g) & act)
        if len(m) > 1:
            same_gate += len(m) * (len(m) - 1) // 2
            if len(np.unique(ug[m])) > 1: split += 1
    rng = np.random.default_rng(2026); a = rng.choice(np.flatnonzero(act), 20000); b = rng.choice(np.flatnonzero(act), 20000)
    checks['mapping'] = dict(active_transistors=int(act.sum()), groups=int(len(reps)),
                             shared_gate_pairs=int(same_gate), shared_gate_nodes_split=int(split),
                             random_pairs_same_group=float(np.mean(ug[a] == ug[b])))

    # 3. the clock: transistors gated by cp1 and cclk in one group, changing in nearly every sample
    clk_nodes = [names[n] for n in ('cp1', 'cclk') if n in names]
    ct = np.flatnonzero(np.isin(gate, clk_nodes) & act)
    cg = np.unique(ug[ct])
    checks['clock'] = dict(clock_nodes=clk_nodes, clock_transistors=int(len(ct)), groups=cg.tolist(),
                           event_rate=[float(red[reps[g]].mean()) for g in cg],
                           group_size=[int((ug == g).sum()) for g in cg])

    # 4. phase structure (descriptive; already seen before any design file, disclosed as such)
    R = red[reps]; par = np.arange(R.shape[1]) % 2; tot = R.sum(1)
    share = np.maximum(R[:, par == 0].sum(1), R[:, par == 1].sum(1)) / np.maximum(tot, 1)
    checks['phase'] = dict(groups=int(len(reps)), locked_95=int((share >= 0.95).sum()), strict=int((share == 1).sum()),
                           mixed=int((share < 0.95).sum()), median_event_rate=float(np.median(R.mean(1))))
    np.savez_compressed(OUT / 'pitfall_groups.npz', representatives=reps, physical_to_group=ug)
    checks['output_sha256'] = hashlib.sha256((ROOT / 'data/raw/pitfall.ts.spikes.npy').read_bytes()).hexdigest()
    json.dump(checks, open(OUT / 'conversion_checks.json', 'w'), indent=1)
    print(json.dumps({k: v for k, v in checks.items() if k not in ('source_sha256',)}, indent=1))


if __name__ == '__main__':
    main()
