# New tests A9 and A10: design and predictions

Written on 3 November 2026 before any statistic below was computed on the
recording, the simulated chip or the new NMOS circuits. Each script prints
the SHA-256 of this file, of itself and of `common.py` at the start of its
log. The only code checks run before the full runs use synthetic data
(`--selftest`) and a load check that prints array shapes only
(`--check-load`); neither computes or prints a statistic of these tests.

Shared definitions are those of `PREDICTIONS.md` (change coding; groups;
strict and 95 % phase locking; stratified AUC; target-cluster bootstrap,
1,000 draws, seed 2026; decision rule: lower 95 % limit > 0.5 supported,
upper < 0.5 contradicted; a paired difference is supported if its lower
limit exceeds 0, contradicted if its upper limit is below 0).

---

## A9 — an exact decomposition of lag-1 TE for strictly locked same-phase pairs

**Status.** Hash-logged and exploratory: the rationale was derived after
A3's matched-bound result (0.400) and A7's wiring result (0.300) had been
seen. Nothing in A9 is held out.

### The identity (derived before the run; checked on synthetic data)

Let target Y and source X be strictly locked to the same parity p, let
Δ denote the change series, k the target's AIS-selected history, and use
exactly the samples t = k … T−1 that the lag-1 TE estimator uses.

- If the target's history (ΔY_{t−1}, …, ΔY_{t−k}) contains a change, the
  history fixes the parity of t, and by Proposition 2 applied cell by cell
  the cell contributes nothing to TE or to the bound B. Only the quiet
  cell Q (all k history samples zero) contributes.
- In Q: if t has parity p, then ΔX_{t−1} = 0; otherwise ΔY_t = 0. The joint
  of (P_t, ΔX_{t−1}, ΔY_t) in Q is therefore fixed by three numbers:
  π = Pr(t has parity p | Q), c = Pr(ΔY_t = 1 | Q, parity p) and
  a = Pr(ΔX_{t−1} = 1 | Q, parity ≠ p).
- With q₀ = Pr(Q), h the binary entropy and m = 1 − (1 − π)a, the plug-in
  estimates satisfy exactly

      TE = q₀ · [ h(πc) − m · h(πc / m) ]          (F)
      B  = q₀ · [ h(π)  − m · h(π / m) ]           (G)

- Only a depends on the source. Within a target, TE is a strictly
  increasing function of a (for 0 < π < 1, 0 < c), so lag-1 TE ranks a
  target's strictly locked same-phase sources exactly by a.
- Writing s = t − 1 and H = {s of parity p with ΔY_{s−1} = … = ΔY_{s−k+1} = 0},
  exactly
      a = r̃ (1 − κ̃) / (1 − c̃),
  where r̃ = Pr(ΔX_s = 1 | H) is the source's activity, κ̃ = Pr(ΔY_s = 1 |
  ΔX_s = 1, H) its **followership** (how often the target changes in the
  same sample as the source), and c̃ = Pr(ΔY_s = 1 | H) depends on the
  target alone. A source that its target copies loses exactly the samples
  in which lag-1 TE could see it.
- Across targets, at matched B, TE/B = F/G. On a grid of π and c it is
  decreasing in a and is governed mainly by c (approximately equal to c
  for small a). Matched-bound comparisons across targets therefore rank
  pairs mainly by the target's c.

### Data and sets

1. **Recording** (primary): the reduced recording, 519 groups, saved
   lag-1 TE (`te.npz`), saved bound B (`A3_bound_arrays.npz`), wiring as in
   A3. Set S = ordered pairs of distinct groups that are strictly locked to
   the same parity. Positives: wired in that direction; negatives: all other
   pairs in S (as in A3).
2. **Simulated chip, program 1 and program 2** (replication): change coding,
   nodes and groups exactly as in test C2 of `chip.py` (intervened nodes U,
   AIS k ≤ 6, node image of the terminal-to-gate reference). S = node pairs
   whose groups are distinct and strictly locked to the same parity.
   Positives: wired; negatives: unwired in both directions (as in C2). TE is
   recomputed with the toolkit estimator. Bootstrap over nodes, as in C2.

### Checks (a failure means a bug, not a result)

- C1: for every pair in S, |TE − (F)| < 1e-9 bits and |B − (G)| < 1e-9 bits
  (recording: against the saved TE and B; chip: against recomputed TE and
  `bound_all_sources`).
- C2: |a − r̃(1 − κ̃)/(1 − c̃)| < 1e-9 for every pair in S.
- C3: within every target, TE is non-decreasing in a.

### Tests (each on the recording, chip program 1 and chip program 2)

- **A9.1 Followership.** Within-target AUC (strata = target) of κ̃, wired
  against unwired. Prediction: above 0.5 (wiring makes the target change
  with its source).
- **A9.2 Followership lowers the rank.** Within-target paired difference
  AUC(TE) − AUC(r̃), wired against unwired: how far lag-1 TE ranks wired
  sources below where their activity alone would place them. Prediction:
  below 0.
- **A9.3 Target composition at matched bound.** Within deciles of B over S,
  AUC of the target's c, wired against unwired. Prediction: below 0.5
  (derived from the identity and the A3 result: wired pairs at matched
  bound go into targets that rarely change after a quiet stretch).
- **A9.4 Accounting (primary).** Paired difference between the AUC of
  lag-1 TE within B-decile × c-quintile strata and its AUC within
  B-decile strata alone. Prediction: above 0. In addition, the
  B × c-stratified AUC is predicted to have an interval that includes 0.5.
  Reading: both hold → the matched-bound inversion on S is accounted for by
  target composition; the difference is above 0 but the stratified AUC's
  upper limit is below 0.5 → partly accounted for, residual reported; the
  difference's interval includes 0 → not accounted for.

### Reported without a prediction

- AUC of lag-1 TE on S: within phase strata only; within B deciles; within
  targets (equal to the within-target AUC of a, by C3).
- Within-target AUCs of r̃ and of 1 − κ̃; within-B-decile AUCs of q₀ and π.
- Scope: number of pairs and wired pairs in S, against all same-phase
  pairs (recording: 315 of 770 wired same-phase pairs are in S).
- For 95 %-locked pairs that are not strictly locked, the distribution of
  |TE − (F)|, as the extent to which the identity approximates them.

---

## A10 — latch transfer in 50 new NMOS circuits (replication of A8)

**Status.** Hash-logged and held out: the circuits (seeds 600–649, input
seeds 613–662) have never been simulated. The predictions for the reader
relation were formed from A8's result and are tested here out of sample.

Design: identical to A8 (`nmos.py`): 80 gates, 6 inputs, 16 latches,
change coding, AIS k ≤ 4, 6,000 and 20,000 half-steps (6,000 is primary),
negatives = pairs with no wiring either way and no latch relation, strata
= circuit × phase classes × rate quintiles, circuit-cluster bootstrap
(1,000 draws, seed 2026).

- **A10.1** Latch input → storage, lag-1 TE: AUC above 0.5.
- **A10.2 (primary)** Latch input → storage, paired lag-1 minus lag-0:
  above 0 (A8's prediction, retested with five times as many circuits;
  A8 gave +0.057 [−0.109, +0.124]).
- **A10.3** Storage → reader output, paired lag-0 minus lag-1: above 0
  (the reader settles within the phase; formed after A8).
- **A10.4** Storage → reader output, lag-1 TE: not above 0.5.
- Reported without prediction: the same four statistics at 20,000
  half-steps; numbers of scorable pairs; pooled A8 + A10 (60 circuits).
