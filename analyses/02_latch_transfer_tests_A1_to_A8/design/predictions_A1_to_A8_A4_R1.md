# New tests: design and predictions (A1, A1b, A1c, A2, A3, A4, A5, A6, A8, R1)

Written on 2 November 2026 before any statistic below was computed. Each
script prints this file's SHA-256 at the start of its log. The author
approved running before reviewing (review on 3 November); the hash shows
whether the file changed after the runs. Full rationale: `THESIS_PLAN.md`
§4–5 and `AGENT_RUNBOOK.md` §5.

## Shared definitions

- Change coding Δ_t = s_t XOR s_{t−1} is primary; node states secondary.
- Units: groups of identical series; pairs within a group are excluded.
- Phase class: even- or odd-locked if ≥ 95% of events fall on one parity of
  the sample (recording) or half-step (chip, NMOS) index; otherwise mixed.
- Matching strata: phase-pair class × source change-rate quintile × target
  change-rate quintile; on the chip, × co-drive-strength quintile (Jaccard
  overlap of E0 cause sets) where stated.
- Stratified AUC: positives compared with negatives only within the same
  stratum, pooled by number of comparisons.
- Uncertainty: target-cluster bootstrap, 1,000 draws, seed 2026 (NMOS:
  circuit-cluster bootstrap).
- Decision rule for a primary AUC: lower 95% limit > 0.5 = supported;
  supported and estimate ≥ 0.60 = material; upper limit < 0.5 =
  contradicted; otherwise not supported. For a paired difference: lower
  limit > 0 = supported.
- Latch pair: a clock-gated transistor whose two channel nodes are both
  signal nodes. Clock-gated: on the recording, gate node is a gate node of
  the clock group (the group with events in ~99% of samples); on the chip,
  gate node whose state series equals `clk0` or its complement. Driven side:
  channel node with a pull-up or a pull-down transistor to ground; storage
  side: neither. Oriented latch pair: driven → storage. On the recording a
  latch pair maps to the groups gated by its two channel nodes.

## A1 — lag-1 TE against delayed influence (simulated chip; core test)

Data: `simulated_chip/states.npz`, `interventions.npz` (E0 same half-step,
E1 after the next, 40 moments). Ordered intervened node pairs i → j.
Classes: **D** (delayed-only) E1[i,j] > 0 and E0[i,j] = 0; **I**
(immediate) E0[i,j] > 0; **N** (none) E0 and E1 zero in both directions.
Scores (AIS-selected k ≤ 6): lag-1 TE, lag-0 CMI, lag-1 TE given half-step
parity. Strata include co-drive strength.

- **P1 (primary):** stratified AUC of lag-1 TE, D vs N, change coding.
- P2: AUC(lag-1, D vs N) − AUC(lag-0, D vs N), paired. Prediction > 0.
- P3: AUC(lag-0, I vs N) − AUC(lag-1, I vs N), paired. Prediction > 0.
- P4: AUC of parity-conditioned lag-1 TE, D vs N. Prediction > 0.5.
- P5: P1 under node states. Prediction > 0.5.
- Expected: P1 > 0.5 if latch transfers dominate D. Risk: D may be small
  or dominated by clock-driven pairs; class sizes are reported first.

## A1b — structural latch pairs (chip)

Oriented latch node pairs vs N pairs; strata without co-drive. Primary:
stratified AUC of lag-1 TE. Secondary: lag-0 CMI (prediction: lower than
lag 1); reverse direction (storage → driven); share of latch pairs in D.

## A1c — IDTxl's chip selections (states; no change-coded run exists)

The 30-target IDTxl runs (lag-1 sources; lags 0 and 1). Precision of
selected group-pair edges, over the node pairs they stand for, against D,
I and wiring; comparator: pairwise lag-1 TE given the same number of
sources per target. Prediction: IDTxl's lag-1 selections more enriched for
D than the comparator. Descriptive (30 targets).

## A2 — latch pairs on the recording (exploratory)

Reference: oriented latch group pairs, recomputed from the netlist (thesis
count for comparison: 82 clock-gated pass-neighbour pairs). Negatives: all
other ordered group pairs. Scores: lag-1 TE, lag-0 CMI, IDTxl lag-1 network
membership, source-group size, minus layout distance, source event rate.
Statistics: stratified AUC (primary for lag-1 TE, with bootstrap), AP lift,
precision at IDTxl's per-target allocation; whole recording and
clock-aligned (samples ≥ 1,408). Prediction: lag-1 TE ranks latch pairs
above matched controls and above each baseline; lag-0 CMI does not.
Exploratory because the reference was defined after the multivariate errors
were seen.

## A3 — Proposition 2 pair by pair (recording)

B_ij = I(ΔX_{t−1}; P_t | Y^(k)) with the target's selected k. For strictly
phase-locked same-parity pairs, TE ≤ B must hold exactly for plug-in
estimates (a violation indicates a bug). Reported: maximum excess; TE/B by
phase-pair class; for 95%-locked pairs, any excess; among same-phase pairs,
stratified AUC of lag-1 TE wired vs unwired within B deciles (prediction
≈ 0.5). Exploratory (bound derived after the main results).

## A4 — IDTxl with the phase in the conditioning set (recording)

Sample parity added as process 519; sources at lag 1 only (parity excluded
from candidates); `add_conditionals = [(519, 1)]`; other settings as the
main run (200 analytic surrogates, for comparability). Targets: the same
100 stratified targets as R1. Comparator: the main run's lag-1 network on
those targets. Primary: group-level precision difference against the
netlist (paired target bootstrap). Secondary: precision against oriented
latch pairs; parent-set sizes; share of opposite-phase sources.
Prediction (honestly uncertain): with the phase supplied, fewer same-phase
sources are admitted; precision against latch pairs rises; precision
against the netlist may not.

## A5 — delay profile

TE at u = 1 … 8 (recording: samples; chip: half-steps, change coding) for
oriented latch pairs, pull-down wired pairs and matched N pairs;
stratified AUC per u against N. Prediction: latch AUC peaks at u = 1;
pull-down pairs near chance at every u.

## A6 — orientation of latch pairs

Fraction of oriented latch pairs with TE(driven → storage) >
TE(storage → driven), against the rule "lower event rate is the target" and
the rule "higher conditional entropy given own past is the target", paired
bootstrap over pairs; recording and chip. Prediction: TE orients latch
pairs better than both rules.

## A8 — NMOS circuits (latches known by construction)

Ten circuits (seeds 500–509, 80 gates, 6 inputs, 16 latches), 6,000 and
20,000 half-steps, change coding. Positives: latch input → storage, and
storage → reader-gate output. Negatives: pairs with no wiring either way
and no latch relation. Stratified AUC (strata include circuit), circuit
bootstrap. Prediction: AUC > 0.5 for lag-1 TE; lag-1 above lag-0 for these
pairs.

## R1 — the stratified H6 design, re-drawn and re-run

100 targets drawn with seed 2026 from the AIS-k bands {1}, {2–4}, {5–8},
{9–11}, {12}, proportional (largest remainder). IDTxl multivariate TE,
sources at lags 0 and 1, settings as the main run. Compare with the thesis
(236 edges, 41 wired, 17.37%). The original target list was not saved, so
this may be a second stratified sample rather than an exact reproduction;
the target list is recorded.
