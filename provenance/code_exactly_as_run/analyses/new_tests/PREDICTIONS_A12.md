# New tests A12–A19: predictions of the followership account, fixed before they are computed

Written on 3 November 2026 (night), before any statistic below was computed. The author
approved running these tests in advance, before reading this file. Each script prints the
SHA-256 of this file, of itself and of `common.py` at the start of its log. Checks before the
full runs use synthetic data, or print only timings and array shapes; none computes a statistic
of these tests on the recording or the chip.

**Status.** All are hash-logged. Their rationale comes from Proposition 2 and Remark 3, which were
derived after the main results; several predictions are partly informed by descriptive results
already reported, and each says so. Shared definitions, decision rule and bootstrap are those of
`PREDICTIONS.md` (target-cluster bootstrap, 1,000 draws, seed 2026; an AUC is supported if its
lower limit exceeds 0.5 and contradicted if its upper limit is below 0.5; a paired difference is
supported if its lower limit exceeds 0).

---

## A12 — the locking dose–response on the recording

Partly anticipated: Table 6.2 already reports lag-1 TE at chance among same-parity pairs and 0.630
among pairs involving a mixed group (clock group excluded). A12 tests the graded, target-level
version.

Targets are binned by their locking index L (share of events in the dominant parity):
**mixed** L < 0.95, **locked** 0.95 ≤ L < 1, **strict** L = 1. The clock group is excluded as source
and target (as in Table 6.2's clock-excluded row). Within-target AUC of lag-1 TE, wired against
pairs unwired in either direction, in each bin.

- **A12.1 (primary).** AUC(mixed) − AUC(strict) > 0 (paired, same bootstrap draws).
- A12.2. AUC(strict) has its upper limit below 0.5.
- A12.3. AUC(mixed) has its lower limit above 0.5.
- Reported: AUC(locked); the same three bins for lag-0 CMI; numbers of targets and wired pairs.

## A13 — the lag 2–7 account of the latch deficit

Partly anticipated: the account was formed to fit the lag-8 values already reported. Its specific
test, splitting by history length, has not been computed.

For latch pairs at source lag u = 2 … 7, the target's history **covers** the copy if k ≥ u − 1.
Lag-u TE, latch pairs against pairs with no influence (chip: no measured influence; recording: no
wiring and no latch relation either way), stratified as in A5 (phase-pair class × rate quintiles)
and additionally by u; covered and uncovered targets scored separately, negatives into targets of
the same class.

- **A13.1 (primary; chip programs 1 and 2 separately).** AUC(covered) has its upper limit below 0.5.
- **A13.2 (primary).** AUC(uncovered) − AUC(covered) > 0 (paired).
- A13.3. AUC(uncovered) is not below chance: its upper limit is at or above 0.5 (the deficit
  disappears when the history does not cover the copy).
- The recording is reported the same way; if fewer than 10 uncovered latch pairs exist on a
  dataset, its uncovered statistics are descriptive.

## A14 — the between-target remainder at matched bound

Partly anticipated: the descriptive AUCs of q₀ (0.549) and π (0.446) within bound deciles are
already reported (A9).

On the strictly locked same-phase pairs (recording; chip programs 1 and 2), lag-1 TE wired against
unwired, within B deciles × q₀ quintiles × π quintiles of the target.

- **A14.1 (primary, recording).** AUC(B × q₀ × π) − AUC(B) > 0 (paired).
- A14.2. The B × q₀ × π AUC's interval includes 0.5.
- Reading: both → the remainder is accounted for by the target's quiet-cell terms; A14.1 only →
  partly; neither → not accounted for.

## A15 — non-idealised generality and coding

New synthetic networks as in A11 with three departures: 10 hidden units (never scored), each
driving 4 observed units of its phase within the sample with probability 0.5, which creates
co-driven siblings; spurious events in either phase with probability 0.01 per sample per unit,
which breaks strict locking; and coupling ρ = 0.6 throughout. λ ∈ {0.5, 0.7, 0.85, 1.0}, 20 networks
each (seeds 20,000 + 100·condition + replicate). Each network scored as events (change coding) and
as states (the cumulative parity of the events).

- **A15.1 (primary).** Events, λ = 1: within-target lag-1 TE AUC for fast drivers against same-phase
  non-drivers has its upper limit below 0.5.
- **A15.2 (primary).** Events: AUC(λ = 1) − AUC(λ = 0.5) < 0 (upper limit below 0).
- A15.3. Coding: at λ = 1, AUC(states) − AUC(events) > 0 (paired; state coding reduces the penalty).
- A15.4. Remedy: lag-0 CMI given parity, events, lower limit above 0.5 in all four conditions.

## A16 — timed clamps: latch transfer by intervention

A1 clamped at random moments and its delayed class held no latch pair. Here each program is
re-simulated (warm-up 200, as originally; checked against the saved states) and clamped at 40
moments when clk0 is high and 40 when it is low (seeds 3 and 4, moments in 500 … 59,998). A latch
pair is scored in the set whose moments find its gate low (the latch closed), so a clamp on its
driven side cannot cross before the next half-step. Classes per set as A1: immediate (E₀ > 0),
delayed-only (E₁ > 0, E₀ = 0), none (both zero in both sets, both directions).

- **A16.1 (primary, each program).** More than half of the latch pairs are delayed-only in their
  closed-latch set.
- A16.2. In the closed-latch sets, pairs that are delayed-only in at least one set and immediate in
  neither: lag-1 TE AUC against none above 0.5 (stratified as A1).
- A16.3. On the same class, lag-1 minus lag-0 CMI > 0 (the A1 secondary prediction that failed).

## A17 — the H6 design on the clock-aligned samples

IDTxl multivariate TE, source lags 0 and 1, every other setting as R1, on R1's 100 targets, on
samples 1,408 onwards. Comparators: the main run's lag-1 network on the same targets (whole
recording), and R1.

- **A17.1 (primary).** Precision against the netlist minus the main run's lag-1 network > 0
  (paired over targets).
- A17.2. Precision minus R1's: interval includes 0 (the gain is not an artefact of the misaligned
  samples).

## A18 — a source block of two samples

Consequence 3 of Proposition 2: a source block of length 2 escapes the bound. The account predicts
the escape still does not recover the wiring, because the wiring acts at lag 0.

TE with source block (ΔX_{t−1}, ΔX_{t−2}) and the target's AIS history, on the recording.

- **A18.1 (primary).** Group-level AUC against the netlist, all ordered group pairs: lower limit
  ≤ 0.5 (no recovery).
- A18.2. Strictly locked same-phase pairs, within-target AUC: lower limit ≤ 0.5.
- Reported: the same for lag-1 TE on the same samples, for comparison.

## A19 — phase-resolved local TE on the latch pairs

For each target, lag-1 TE is split into the share contributed by samples of the target's own phase
(where it can change) and of the other phase (where the source can only announce the phase).
Prediction: transfer arises where the target can change; phase information where it cannot.

- **A19.1 (primary).** Within the storage targets of the 41 latch pairs, the in-phase share of
  lag-1 TE ranks latch sources above the other sources with positive TE into the same targets
  (within-target AUC, lower limit above 0.5).
- A19.2. Median in-phase share: above 0.5 for latch pairs; below 0.5 for strictly locked same-phase
  pairs.
