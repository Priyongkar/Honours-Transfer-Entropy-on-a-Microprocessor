"""Supplementary (post hoc, after the whole-thesis marking): do lag-1 TE scores single out the
recording's 41 oriented latch pairs among the other pass-neighbour pairs? Answers the specialist
marker's circularity point. Netlist code as in recording.py. Writes output/supp_latch_vs_passnb.json."""
import ast
import numpy as np
from common import RUN, ROOT, OUT, tk, save, boot_mult, auc_with_ci
VSS, VCC = 558, 657

X = np.asarray(np.load(RUN / 'thesis_v5/idtxl/data_int8.npy')).astype(np.uint8)
with np.load(RUN / 'te.npz') as z:
    TE1 = z['te_bits'].copy()
p2g = np.load(RUN / 'groups.npz')['physical_to_group'].astype(int)
ref = np.load(OUT / 'recording_reference.npz'); Lat, Wd = ref['latch'], ref['wired']
G = X.shape[0]
defs = []
for line in (ROOT / 'data/visual6502/transdefs.js').read_text().splitlines():
    if line.strip().startswith("['t"):
        defs.append(ast.literal_eval(line.strip().rstrip(',')))
defs.sort(key=lambda d: int(d[0][1:]))
gate, c1, c2 = (np.asarray([d[i] for d in defs]) for i in (1, 2, 3))
gated_by = {}
for t, g in enumerate(gate):
    gated_by.setdefault(int(g), []).append(t)
groups_of = lambda nd: {int(p2g[t]) for t in gated_by.get(int(nd), []) if p2g[t] >= 0}
rate = X.mean(1); clock_g = int(np.argmax(rate))
clock_gates = {int(gate[t]) for t in np.flatnonzero(p2g == clock_g)}
PN = np.zeros((G, G), bool); PNclk = np.zeros((G, G), bool)
for t in range(len(gate)):
    if c1[t] in (VSS, VCC) or c2[t] in (VSS, VCC):
        continue
    for ga in groups_of(c1[t]):
        for gb in groups_of(c2[t]):
            if ga != gb:
                PN[ga, gb] = PN[gb, ga] = True
                if int(gate[t]) in clock_gates:
                    PNclk[ga, gb] = PNclk[gb, ga] = True
PN &= ~(Wd | Wd.T); PNclk &= PN
res = dict(pass_neighbour_pairs=int(PN.sum()), clock_gated=int(PNclk.sum()), latch=int(Lat.sum()),
           latch_within_pass_neighbours=int((Lat & PN).sum()))
mult = boot_mult(G, 1000, 2026)
off = ~np.eye(G, dtype=bool)
ii, jj = np.nonzero(off)
for name, neg in [('vs_non_clock_gated_pass_neighbours', PN & ~PNclk),
                  ('vs_all_other_pass_neighbours', PN & ~Lat)]:
    pos = Lat[ii, jj]; m = pos | neg[ii, jj]
    r, _ = auc_with_ci(TE1[ii, jj][m], pos[m], np.zeros(m.sum(), int), jj[m], mult)
    r['n_neg'] = int(neg.sum()); res[name] = r
save('supp_latch_vs_passnb', res)
print(res)
