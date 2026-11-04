# Settling versus co-drive: design and predictions

Written before any statistic below was computed. The SHA-256 of this file is
printed at the start of the analysis output, so later edits are detectable.

## Question

On the recording, lag-0 CMI ranks wired pairs above unwired ones (H2,
comparative clause). Two accounts predict that:

- **Settling (direct coupling).** The source changes and, through the wire,
  settles the target within the same phase.
- **Co-drive.** The source and the target are both settled in that phase by a
  common upstream change; the wire adds nothing beyond that.

The recording cannot separate them (no timing below one phase, no
interventions). The simulated chip can, because its clamping reference E0
says which node changes which within the same half-step.

## Data

Simulated 6502, `results/.../thesis_v5/simulated_chip/states.npz` (60,000
half-steps of every node state, Fibonacci/LFSR program) and
`interventions.npz` (E0: clamp each of the 1,100 non-constant nodes at 40
moments; fraction of moments at which each node then differs within the same
half-step). Nothing is re-simulated.

## Pair classes (ordered node pairs i -> j, both intervened, both non-constant)

- **Direct:** i -> j wired (terminal-to-gate, node image) and E0[i, j] > 0.
- **Co-driven only:** not wired in either direction, E0[i, j] = 0 and
  E0[j, i] = 0, and at least one common cause: an intervened node c (not i or
  j) with E0[c, i] > 0 and E0[c, j] > 0.
- **Wired but idle:** i -> j wired and E0[i, j] = 0 (secondary contrast).

Co-drive strength of a pair = Jaccard overlap of the cause sets of i and j
(causes of u = intervened nodes c != i, j with E0[c, u] > 0).

## Representations and scores

Change coding (Δ_t = s_t XOR s_{t-1}; the recording's coding) is primary;
node states are secondary. Nodes with identical series are grouped as on the
recording and pairs within one group are excluded (no measure can order
them). Score: lag-0 CMI I(X_t; Y_t | Y^(k)) with the AIS-selected k ≤ 6 of
the target's group, plug-in, bits (toolkit `te_all_sources`, delay 0).
Secondary scores: lag-0 CMI conditioned on half-step parity, and lag-1 TE.

## Primary statistic

P1: AUC of lag-0 CMI separating **direct** from **co-driven-only** pairs,
with comparisons only within strata defined by: phase class of source and of
target (even-locked, odd-locked, mixed; locked = at least 95% of changes on
one parity), source change-rate quintile, target change-rate quintile and
co-drive-strength quintile. 95% target-cluster bootstrap, 1,000 draws,
seed 2026.

## Decision rule (fixed now)

- Lower 95% limit of P1 above 0.5: direct coupling leaves lag-0 dependence
  beyond what matched co-drive produces (supports the settling reading of
  H2's causal clause on the chip). Called **material** if the point estimate
  is at least 0.60.
- Upper limit below 0.5: co-driven pairs score higher than direct pairs
  matched on co-drive (against the settling reading).
- Otherwise: lag-0 CMI does not distinguish settling from co-drive on the
  chip.

## Secondary (descriptive, no verdict)

- P1 under node states.
- P1 with parity-conditioned lag-0 CMI.
- P2: Spearman correlation of lag-0 CMI with co-drive strength among
  co-driven-only pairs (does co-drive alone generate lag-0 dependence?).
- P3: stratified AUC of lag-0 CMI, direct versus wired-but-idle.
- Unstratified AUCs and class means, for scale.

## Expectation (stated so it can be wrong)

I expect P1 above 0.5 under states and closer to 0.5 under change coding,
because change coding merges inverter inputs and outputs into one group and
removes the strongest direct pairs.

## Limits known in advance

Different program from the recording; E0 has 40 moments, so a pair coded
"no influence" may have influence at moments not sampled; the clock nodes
are common causes of almost everything, which the co-drive-strength strata
absorb rather than exclude.
