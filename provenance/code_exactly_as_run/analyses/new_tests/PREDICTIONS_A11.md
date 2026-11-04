# New test A11: does the followership penalty generalise beyond the 6502?

Written on 3 November 2026 before any network below was simulated. The author
approved the design in advance, before reading this file. The script prints the
SHA-256 of this file, of itself, of `common.py` and of `a9_decompose.py` at the
start of its log. The only check before the full run (`--selftest`) simulates one
network with a seed outside the test set, checks the Remark 3 identity on it, and
prints timing only; no AUC of the test is computed before the full run.

**Status.** Hash-logged. Synthetic: the system is new and no prediction was
formed from it, but the mechanism tested (Remark 3, A9) was derived from the
recording, so A11 is a test of generality, not an independent discovery.

## Question

On the 6502, wiring acts within the sample, units lock to the clock phase, and
lag-1 transfer entropy penalises the sources its targets follow (A9). Is that a
property of the chip, or of any system in which (i) coupling completes within one
sample and (ii) units lock to a shared rhythm sampled at half its period, as
neurons lock to an oscillation? And does the strength of locking, which can be
measured without knowing the wiring, predict when a lagged measure misleads?

## Model

Networks of N = 100 binary event units over T = 50,000 samples. Sample parity
plays the rhythm's phase. Each unit is assigned a preferred phase (50 per phase)
and a base rate r ~ U(0.02, 0.10) per sample.

- **Locking strength λ.** A unit's spontaneous events fall in its preferred phase
  with probability 2λr and in the other with probability 2(1 − λ)r per sample,
  so a fraction λ of them fall in the preferred phase. λ = 0.5 is no locking;
  λ = 1 is strict locking.
- **Layers.** 40 source units (L0), 30 units of L1 and 30 of L2.
- **Fast edges (within the sample).** Each L1 unit has 2 drivers in L0, and each
  L2 unit 2 drivers in L0 ∪ L1, all of its own preferred phase. When a driver has
  an event at t, the target has one at t with probability ρ (the coupling, which
  sets the followership).
- **Delayed edges (one sample).** Each L1 and L2 unit, with probability 0.5, has
  one delayed driver in L0 of the opposite phase: a driver event at t − 1 gives a
  target event at t with probability ρ. These are the latch analogue.
- Events combine by OR. Units are generated in layer order, so every event is
  final before it is propagated.

**Conditions.** λ ∈ {0.5, 0.7, 0.85, 0.95, 1.0} × ρ ∈ {0.3, 0.6, 0.9};
20 networks per condition (300 in all); network seeds 10,000 + 100·(condition
index) + replicate.

## Scores and comparisons

Per network, as on the chip: AIS-selected target history k ≤ 6; lag-1 TE; lag-0
CMI; lag-0 CMI conditioned on parity (the remedy). Pairs are pooled over a
condition's 20 networks. Units for the target-cluster bootstrap are
(network, target), 2,000 per condition; 1,000 draws, seed 2026.

- **Fast edges:** positives are fast driver → target pairs; negatives are pairs of
  the same preferred phase with no edge in either direction (co-driven siblings
  and two-step pairs included, as in a real network).
- **Delayed edges:** positives are delayed driver → target pairs; negatives are
  pairs of opposite preferred phase with no edge in either direction.
- All AUCs are within target (strata = network × target).

## Predictions

- **A11.1 (primary).** With strict locking (λ = 1) and coupling ρ = 0.6 and 0.9,
  lag-1 TE ranks fast drivers *below* same-phase non-drivers: within-target AUC
  with upper 95 % limit below 0.5, in both conditions.
- **A11.2 (primary).** Without locking (λ = 0.5) lag-1 TE is blind to fast edges:
  the AUC interval includes 0.5 for every ρ. And locking makes it worse: for
  ρ = 0.6 and 0.9, AUC(λ = 1) − AUC(λ = 0.5) has upper limit below 0
  (independent bootstraps of the two conditions).
- **A11.3 (remedy).** Lag-0 CMI conditioned on parity ranks fast drivers above
  same-phase non-drivers in all 15 conditions (lower limit above 0.5).
- **A11.4 (positive control).** Lag-1 TE ranks delayed drivers above
  opposite-phase non-drivers in all 15 conditions (lower limit above 0.5).

Reading: if A11.1 and A11.2 hold, locking turns a lagged measure that is merely
blind to within-sample coupling into one that is misleading, outside the 6502.

## Reported without prediction

- Fast-edge lag-1 AUC for ρ = 0.3 and at λ = 0.7, 0.85, 0.95 (the dose–response).
- Unconditioned lag-0 CMI for fast edges.
- **Diagnostic.** Per condition, the mean phase-locking index (each unit's
  fraction of events in its dominant parity, measured from the events alone,
  averaged over units) against the fast-edge lag-1 AUC; Spearman correlation over
  the 15 conditions.
- **Check.** On every λ = 1 network, the Remark 3 identity on strictly locked
  same-phase pairs (|TE − formula| < 1e-9 bits; a failure is a bug), and the
  within-target AUC of followership κ̃ for fast edges.
