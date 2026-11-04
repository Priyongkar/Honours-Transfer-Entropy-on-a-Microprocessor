"""A7: second 6502 program (design: PREDICTIONS_A7.md). Assembles the program,
validates the switch-level run against py65, records states and clamp effects,
then runs chip.py on the new data.
Usage: python a7_program.py [--test]"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / 'src'))  # repository library
from te6502.paths import ANALYSES, DATA_ROOT, RESULTS, HELPERS_FILE, load_analysis  # noqa: E402,F401
import argparse, subprocess, sys, time
import numpy as np
from py65.assembler import Assembler
from py65.devices.mpu6502 import MPU
from te6502.test_helpers import HERE, ROOT, OUT, log
from te6502.switch_level_simulator import Chip

ORIGIN = 0x0600
SRC = """
start:  LDX #$FF
        TXS
        LDA #$5A
        STA $10
        LDA #$00
        STA $11
loop:   LDY #$00
gen:    LDA $10
        ASL A
        ASL A
        CLC
        ADC $10
        CLC
        ADC #$11
        STA $10
        STA $0300,Y
        INY
        CPY #$10
        BNE gen
        LDA #$0F
        STA $12
pass:   LDX #$00
inner:  LDA $0300,X
        CMP $0301,X
        BCC noswap
        BEQ noswap
        TAY
        LDA $0301,X
        STA $0300,X
        TYA
        STA $0301,X
noswap: INX
        CPX #$0F
        BNE inner
        DEC $12
        BNE pass
        JSR chk
        LDX $11
        STA $0200,X
        LDA $030F
        STA $0280,X
        INC $11
        LDA $11
        CMP #$40
        BNE again
        LDA #$00
        STA $11
again:  JMP loop
chk:    LDA #$00
        STA $14
        LDA #$03
        STA $15
        LDA #$00
        STA $13
        LDY #$00
cl:     SEC
        SBC ($14),Y
        PHA
        EOR $13
        STA $13
        PLA
        INY
        CPY #$10
        BNE cl
        EOR $13
        RTS
"""


def assemble():
    lines = []
    for raw in SRC.strip().splitlines():
        label, _, stmt = raw.partition(':') if ':' in raw.split()[0] else ('', '', raw)
        lines.append((label.strip(), stmt.strip()))
    asm = Assembler(MPU())
    labels = {}
    for pas in (1, 2):
        pc = ORIGIN; code = bytearray()
        for label, stmt in lines:
            if label:
                labels[label] = pc
            s = stmt
            for name, addr in labels.items():
                s = s.replace(name, f'${addr:04X}')
            for name in [l for l, _ in lines if l and l not in labels]:
                s = s.replace(name, f'${pc:04X}')          # forward reference, pass 1 only
            b = asm.assemble(s, pc)
            if b is None:
                raise ValueError(f'cannot assemble {stmt!r}')
            code += bytes(b); pc += len(b)
    return bytes(code), labels


def reference(code, n_instr=400000):
    mpu = MPU(); mpu.memory[ORIGIN:ORIGIN + len(code)] = list(code); mpu.pc = ORIGIN
    vals, last = [], 0
    for _ in range(n_instr):
        mpu.step()
        if mpu.memory[0x11] != last:              # iteration counter advanced
            i = last; vals.append((mpu.memory[0x0200 + i], mpu.memory[0x0280 + i])); last = mpu.memory[0x11]
            if last == 0:
                break
    return vals


def simulate(code, warmup, half_steps):
    chip = Chip(ROOT / 'data/visual6502')
    chip.reset(program=code, origin=ORIGIN); chip.run(warmup)
    states, ann = chip.run(half_steps)
    return chip, states, ann


def validate(mem, ref):
    n = mem[0x11]
    ok = n > 0 and all((mem[0x0200 + i], mem[0x0280 + i]) == ref[i] for i in range(n))
    m2 = bytearray(mem); m2[0x0200] ^= 1
    flipped_fails = not all((m2[0x0200 + i], m2[0x0280 + i]) == ref[i] for i in range(n))
    return n, ok, flipped_fails


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--test', action='store_true'); a = ap.parse_args()
    import hashlib
    log(f"A7: PREDICTIONS_A7.md sha256 {hashlib.sha256((ANALYSES / '02_latch_transfer_tests_A1_to_A8/design/predictions_A7_second_program.md').read_bytes()).hexdigest()}")
    code, labels = assemble(); ref = reference(code)
    log(f'assembled {len(code)} bytes; py65 reference iterations {len(ref)}')
    T = 24000 if a.test else 60000
    t0 = time.perf_counter(); chip, states, ann = simulate(code, 200, T)
    n, ok, flipped = validate(chip.memory, ref)
    log(f'simulated {T} half-steps in {time.perf_counter() - t0:.0f}s; iterations {n}; matches py65 {ok}; altered byte fails {flipped}')
    if not (ok and flipped):
        log('VALIDATION FAILED; stopping'); sys.exit(1)
    if a.test:
        return
    d = OUT / 'A7_chip'; d.mkdir(exist_ok=True)
    np.savez_compressed(d / 'states.npz', states=states, memory=np.frombuffer(bytes(chip.memory), np.uint8), **ann)
    # clamp interventions, as notebook cell 86 (40 moments, seed 2)
    S = states; varying = np.flatnonzero(S.min(0) != S.max(0)); rec = ann['nodes']
    times = np.sort(np.random.default_rng(2).choice(np.arange(500, T - 2), 40, replace=False))
    chip = Chip(ROOT / 'data/visual6502'); chip.reset(program=code, origin=ORIGIN); chip.run(200)
    snaps, prev = [], 0
    for t in times:
        for _ in range(t - prev): chip.half_step()
        prev = t; snaps.append(chip.snapshot())
    state_of = lambda: np.fromiter((chip.state[x] for x in rec), dtype=bool, count=len(rec))
    E0 = np.zeros((len(varying), len(rec)), np.float32); E1 = np.zeros_like(E0)
    t0 = time.perf_counter()
    for k, snap in enumerate(snaps):
        chip.restore(snap); chip.clamp = {}; b0 = state_of(); chip.half_step(); b1 = state_of()
        for i, v in enumerate(varying):
            chip.restore(snap); chip.clamp = {}
            node = int(rec[v]); chip.set_clamp(node, not chip.state[node])
            E0[i] += state_of() != b0
            chip.half_step(); E1[i] += state_of() != b1
            chip.clamp = {}
        log(f'  moment {k + 1}/40 ({time.perf_counter() - t0:.0f}s)')
    np.savez_compressed(d / 'interventions.npz', E0=E0 / 40, E1=E1 / 40, nodes=varying, n_times=40)
    log('interventions saved; running chip analysis')
    r = subprocess.run([sys.executable, '-u', str(ANALYSES / '02_latch_transfer_tests_A1_to_A8/chip_latch_and_intervention_tests_A1_A1b_A5_A6.py'), '--chip-dir', str(d), '--tag', 'A7'],
                       cwd=ANALYSES / '02_latch_transfer_tests_A1_to_A8', stdout=open(OUT / 'A7_chip.log', 'w'), stderr=subprocess.STDOUT)
    log(f'chip analysis exit {r.returncode}')


if __name__ == '__main__':
    main()
