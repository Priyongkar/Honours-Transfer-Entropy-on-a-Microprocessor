"""sim6502: a Python port of the Visual6502 switch-level simulator (chipsim.js,
Greg James, Brian Silverman & Barry Silverman, 2010; netlist CC BY-NC-SA 3.0).

The algorithm is the published Visual6502 one: transistors are switches; when
a node changes, the connected group of nodes (through 'on' transistors) is
resolved to a single value (ground > power > pull-up/pull-down > any high
member), and transistors gated by changed nodes are switched, repeating until
quiescent. Additions for the thesis:
  * per-half-step recording of every node state;
  * node clamping (in-silico lesions / interventions);
  * fast state snapshots (save/restore) for short-horizon perturbations;
  * instruction annotation (SYNC, address bus, opcode) for alignment.
"""
from __future__ import annotations
import re
import numpy as np

PROGRAM_FIB_LFSR = bytes([   # variables live at $80-$82, outside the code
    0xA2, 0x00, 0xA9, 0x01, 0x85, 0x80, 0xA9, 0x00, 0x85, 0x81,       # LDX #0; LDA #1; STA $80; LDA #0; STA $81
    0xA5, 0x80, 0x9D, 0x00, 0x02, 0x18, 0x65, 0x81, 0xA4, 0x80,       # loop: LDA $80; STA $0200,X; CLC; ADC $81; LDY $80
    0x84, 0x81, 0x85, 0x80, 0xE8, 0xE0, 0x20, 0xD0, 0xED,             # STY $81; STA $80; INX; CPX #$20; BNE loop
    0xA9, 0x01, 0x85, 0x82, 0xA0, 0x00,                               # LDA #1; STA $82; LDY #0
    0xA5, 0x82, 0x0A, 0x90, 0x02, 0x49, 0x1D, 0x85, 0x82,             # lf: LDA $82; ASL A; BCC +2; EOR #$1D; STA $82
    0x99, 0x00, 0x03, 0xC8, 0x4C, 0x23, 0x00,                         # STA $0300,Y; INY; JMP lf
])


def reference_fib_lfsr(n_lfsr=64):
    """Expected memory contents written by PROGRAM_FIB_LFSR (Python reference)."""
    a, b, fib = 1, 0, []
    for _ in range(32):
        fib.append(a); a, b = (a + b) & 0xFF, a
    x, lfsr = 1, []
    for _ in range(n_lfsr):
        carry = x & 0x80; x = (x << 1) & 0xFF
        if carry: x ^= 0x1D
        lfsr.append(x)
    return fib, lfsr


def load_netlist(folder):
    folder = str(folder)
    seg = open(folder + '/segdefs.js').read()
    trans = open(folder + '/transdefs.js').read()
    names = open(folder + '/nodenames.js').read()
    pullup = {}
    for m in re.finditer(r"\[\s*(\d+)\s*,\s*'([+-])'", seg):
        n = int(m.group(1)); pullup[n] = pullup.get(n, False) or m.group(2) == '+'
    tdefs = [(int(a), int(b), int(c), int(d)) for a, b, c, d in
             re.findall(r"\['t(\d+)',\s*(\d+),\s*(\d+),\s*(\d+)", trans)]
    tdefs.sort()
    nodenames = {k: int(v) for k, v in re.findall(r"^\s*([A-Za-z_][\w#]*)\s*:\s*(\d+)", names, re.M)}
    return pullup, tdefs, nodenames


class Chip:
    def __init__(self, netlist_folder=None, definitions=None):
        """Build from a Visual6502 netlist folder, or from definitions =
        (pullup dict, [(transistor_id, gate, c1, c2)], names dict with 'vss'/'vcc')."""
        pullup, tdefs, names = load_netlist(netlist_folder) if definitions is None else definitions
        self.names = dict(names)
        self.gnd, self.pwr = self.names['vss'], self.names['vcc']
        n_nodes = max(max(pullup), max(max(t[1:]) for t in tdefs)) + 1
        self.n_nodes = n_nodes
        self.pullup = [bool(pullup.get(i, False)) for i in range(n_nodes)]
        self.pulldown = [False] * n_nodes
        self.state = [False] * n_nodes
        self.clamp = {}                      # node -> forced value (intervention)
        self.t_gate, self.t_c1, self.t_c2 = [], [], []
        self.t_on = []
        self.gates = [[] for _ in range(n_nodes)]   # transistors gated by node
        self.c1c2s = [[] for _ in range(n_nodes)]   # transistors with a channel terminal at node
        for idx, g, c1, c2 in tdefs:
            if c1 == self.gnd: c1, c2 = c2, self.gnd
            if c1 == self.pwr: c1, c2 = c2, self.pwr
            t = len(self.t_gate)
            self.t_gate.append(g); self.t_c1.append(c1); self.t_c2.append(c2); self.t_on.append(False)
            self.gates[g].append(t); self.c1c2s[c1].append(t); self.c1c2s[c2].append(t)
        self.transistor_ids = [d[0] for d in tdefs]
        self.used_nodes = sorted({n for t in tdefs for n in t[1:]} - {self.gnd, self.pwr})
        self.memory = bytearray(65536)
        self.half_steps = 0

    # ---------------- core switch-level algorithm (chipsim.js) ----------------
    def _group(self, start):
        group = {start}; stack = [start]
        gnd, pwr = self.gnd, self.pwr
        t_on, t_c1, t_c2, c1c2s = self.t_on, self.t_c1, self.t_c2, self.c1c2s
        while stack:
            i = stack.pop()
            if i == gnd or i == pwr:
                continue
            for t in c1c2s[i]:
                if not t_on[t]:
                    continue
                other = t_c2[t] if t_c1[t] == i else t_c1[t]
                if other not in group:
                    group.add(other); stack.append(other)
        return group

    def _value(self, group):
        if self.clamp:
            for n in group:
                if n in self.clamp:
                    return self.clamp[n]
        if self.gnd in group: return False
        if self.pwr in group: return True
        pullup, pulldown, state = self.pullup, self.pulldown, self.state
        for n in group:
            if pullup[n]: return True
            if pulldown[n]: return False
        for n in group:
            if state[n]: return True
        return False

    def recalc(self, nodes):
        gnd, pwr = self.gnd, self.pwr
        work = [n for n in dict.fromkeys(nodes) if n != gnd and n != pwr]
        for _ in range(100):
            if not work:
                return
            pending = {}
            for node in work:
                group = self._group(node)
                value = self._value(group)
                for n in group:
                    if self.state[n] == value:
                        continue
                    self.state[n] = value
                    for t in self.gates[n]:
                        if value and not self.t_on[t]:
                            self.t_on[t] = True; pending[self.t_c1[t]] = 1
                        elif not value and self.t_on[t]:
                            self.t_on[t] = False; pending[self.t_c1[t]] = 1; pending[self.t_c2[t]] = 1
            work = [n for n in pending if n != gnd and n != pwr]

    # ---------------- pins, buses, clocking ----------------
    def set_pin(self, name, high):
        n = self.names[name] if isinstance(name, str) else name
        self.pullup[n] = bool(high); self.pulldown[n] = not high
        self.recalc([n])

    def high(self, name):
        return self.state[self.names[name]]

    def read_bits(self, prefix, n):
        return sum(int(self.state[self.names[f'{prefix}{i}']]) << i for i in range(n))

    def write_data_bus(self, x):
        nodes = []
        for i in range(8):
            n = self.names[f'db{i}']; bit = (x >> i) & 1
            self.pullup[n] = bool(bit); self.pulldown[n] = not bit; nodes.append(n)
        self.recalc(nodes)

    def half_step(self):
        if self.high('clk0'):
            self.set_pin('clk0', False)
            if self.high('rw'):
                self.write_data_bus(self.memory[self.read_bits('ab', 16)])
        else:
            self.set_pin('clk0', True)
            if not self.high('rw'):
                self.memory[self.read_bits('ab', 16)] = self.read_bits('db', 8)
        self.half_steps += 1

    def reset(self, program=PROGRAM_FIB_LFSR, origin=0x0000):
        self.memory[:] = bytes(65536)
        self.memory[origin:origin + len(program)] = program
        self.memory[0xFFFC] = origin & 0xFF; self.memory[0xFFFD] = origin >> 8
        self.state = [False] * self.n_nodes; self.state[self.pwr] = True
        self.t_on = [False] * len(self.t_on)
        for name, v in [('res', False), ('clk0', False), ('rdy', True), ('so', False), ('irq', True), ('nmi', True)]:
            n = self.names[name]; self.pullup[n] = v; self.pulldown[n] = not v
        self.recalc(range(self.n_nodes))
        for _ in range(8):
            self.set_pin('clk0', True); self.set_pin('clk0', False)
        self.set_pin('res', True)
        self.half_steps = 0

    # ---------------- snapshots and interventions ----------------
    def snapshot(self):
        return (self.state.copy(), self.t_on.copy(), self.pullup.copy(), self.pulldown.copy(),
                bytes(self.memory), self.half_steps)

    def restore(self, snap):
        s, t, pu, pd, mem, hs = snap
        self.state = s.copy(); self.t_on = t.copy(); self.pullup = pu.copy(); self.pulldown = pd.copy()
        self.memory[:] = mem; self.half_steps = hs

    def set_clamp(self, node, value):
        """Force a node (and its connected group) to a value: an in-silico lesion."""
        self.clamp[int(node)] = bool(value)
        self.recalc([int(node)])

    def clear_clamps(self):
        nodes = list(self.clamp); self.clamp = {}
        self.recalc(nodes)

    # ---------------- recording ----------------
    def run(self, half_steps, record_nodes=None, progress=0):
        """Run and record node states after every half-step.
        Returns (states [T, len(record_nodes)] uint8, annotation dict)."""
        rec = np.asarray(self.used_nodes if record_nodes is None else record_nodes, dtype=np.int64)
        out = np.zeros((half_steps, len(rec)), dtype=np.uint8)
        clk = np.zeros(half_steps, np.uint8); sync = np.zeros(half_steps, np.uint8)
        ab = np.zeros(half_steps, np.int32); db = np.zeros(half_steps, np.int16); rw = np.zeros(half_steps, np.uint8)
        state = self.state
        for t in range(half_steps):
            self.half_step(); state = self.state
            out[t] = [state[n] for n in rec]
            clk[t] = self.high('clk0'); sync[t] = self.high('sync'); rw[t] = self.high('rw')
            ab[t] = self.read_bits('ab', 16); db[t] = self.read_bits('db', 8)
            if progress and (t + 1) % progress == 0:
                print(f'    simulated {t + 1}/{half_steps} half-steps', flush=True)
        return out, dict(nodes=rec, clk0=clk, sync=sync, rw=rw, address=ab, data=db)

    def transistor_matrix(self):
        """[transistor] -> (gate, c1, c2) in original transdefs order (supply-swapped)."""
        return np.asarray(self.t_gate), np.asarray(self.t_c1), np.asarray(self.t_c2)


# ---------------------------------------------------------------------------
# Synthetic NMOS circuits simulated with the same switch-level engine
# ---------------------------------------------------------------------------

def make_nmos_circuit(n_gates=60, n_inputs=6, n_latches=12, seed=0, gate_mix=(('INV', .3), ('NOR2', .3), ('NOR3', .15), ('NAND2', .25))):
    """Random NMOS logic with depletion pull-ups, two-phase clock and dynamic
    latches (pass transistors gated by phi1/phi2), built like the 6502:
    gate outputs have pull-ups; NOR = parallel pull-downs; NAND = series
    pull-downs through an internal node; latches store charge on unpulled nodes.
    Returns definitions for Chip(definitions=...) plus metadata."""
    rng = np.random.default_rng(seed)
    VSS, VCC, PHI1, PHI2 = 0, 1, 2, 3
    nxt = [4]
    def new():
        n = nxt[0]; nxt[0] += 1; return n
    inputs = [new() for _ in range(n_inputs)]
    latches = [new() for _ in range(n_latches)]      # storage nodes (no pull-up)
    pullup = {VSS: False, VCC: False}; tdefs = []; gates = []
    def tr(g, c1, c2):
        tdefs.append((len(tdefs), g, c1, c2))
    pool = list(inputs) + list(latches)
    kinds, probs = zip(*gate_mix); probs = np.asarray(probs) / sum(probs)
    for _ in range(n_gates):
        kind = str(rng.choice(kinds, p=probs)); out = new(); pullup[out] = True
        k = {'INV': 1, 'NOR2': 2, 'NOR3': 3, 'NAND2': 2}[kind]
        ins = [int(v) for v in rng.choice(pool, size=min(k, len(pool)), replace=False)]
        if kind == 'NAND2':
            mid = new(); pullup[mid] = False; tr(ins[0], out, mid); tr(ins[1], mid, VSS)
        else:
            for g in ins:
                tr(g, out, VSS)
        gates.append(dict(kind=kind, output=out, inputs=ins)); pool.append(out)
    outputs = [g['output'] for g in gates]
    latch_src = []
    for i, s_ in enumerate(latches):   # latch input: a gate output (feedback allowed), clocked by phi1 or phi2
        src = int(rng.choice(outputs)); clk = PHI1 if i % 2 == 0 else PHI2
        tr(clk, s_, src); latch_src.append((s_, src, clk)); pullup.setdefault(s_, False)
    for n in range(4, nxt[0]):
        pullup.setdefault(n, False)
    names = {'vss': VSS, 'vcc': VCC, 'phi1': PHI1, 'phi2': PHI2}
    for i, n in enumerate(inputs):
        names[f'in{i}'] = n
    return (pullup, tdefs, names), dict(inputs=inputs, latches=latch_src, gates=gates, n_nodes=nxt[0])


def simulate_nmos(definitions, meta, half_steps, seed=0, input_persistence=0.9):
    """Two-phase non-overlapping clock; inputs change once per cycle with persistence."""
    rng = np.random.default_rng(seed)
    chip = Chip(definitions=definitions)
    chip.state = [False] * chip.n_nodes; chip.state[chip.pwr] = True
    chip.recalc(range(chip.n_nodes))
    vals = rng.integers(0, 2, len(meta['inputs']))
    out = np.zeros((half_steps, len(chip.used_nodes)), np.uint8)
    for t in range(half_steps):
        phase = t % 2
        if phase == 0:
            keep = rng.random(len(vals)) < input_persistence
            vals = np.where(keep, vals, rng.integers(0, 2, len(vals)))
            for i, n in enumerate(meta['inputs']):
                chip.pullup[n] = bool(vals[i]); chip.pulldown[n] = not vals[i]
            chip.recalc(meta['inputs'])
        chip.set_pin('phi1', phase == 0)
        chip.set_pin('phi2', phase == 1)
        out[t] = [chip.state[n] for n in chip.used_nodes]
    return out, chip
