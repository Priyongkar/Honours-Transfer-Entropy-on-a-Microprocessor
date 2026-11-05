# R7: the locking gradient within one wiring class, and the matched-bound inversion matched on history

Written on 5 November 2026, before either statistic was computed. `run_locking_within_class.py` prints
the SHA-256 of this file and of itself at the start of its log.

## Why

Two limits were added to the thesis after its independent read. (1) The locking gradient (A12: targets
active in both phases score above strictly locked ones under lag-1 TE) compares targets that also differ
in what reaches them, since targets active in both phases receive more cross-phase input, so the gradient
could be a wiring-class effect. (2) The plug-in bias of lag-1 TE varies across targets with the selected
history k, and the across-target inversion at matched bound (A3: 0.400) compares pairs from targets with
different k, so bias could contribute to it.

## What has been seen

A12 (mixed 0.604, locked 0.449, strict 0.449 within-target AUC on Donkey Kong), the held-out gradients,
A3 (0.400 at matched bound) and A9/A14. Neither statistic below has been computed.

## Design

Donkey Kong, reduced recording, 519 groups, lag-1 TE at the selected histories (the main run's matrix).
(a) As A12 (clock group excluded; negatives the unwired pairs into the target; within-target AUC with
    the target as stratum; 95% target-cluster intervals, 1,000 draws, seed 2026), but with positives
    restricted to one class of wired pair: **pull-down** (at least one pull-down edge, from
    `recording_reference.npz`) or **pass-only** (wired, no pull-down edge). Locking bins as A12.
(b) As A3 (strictly locked same-phase pairs; lag-1 TE against wired within deciles of the bound), with
    the strata refined to decile of the bound x the target's selected k.

## Predictions (the account's)

- **R7.1** Pull-down pairs: within-target AUC for mixed targets minus strictly locked targets, lower
  limit above 0.
- **R7.2** Pass-only pairs: the same, lower limit above 0.
- **R7.3** Matched on bound and k, the AUC's upper limit stays below 0.5 (the inversion is not an
  artefact of bias differing with k).

Descriptive: each class's strict-bin AUC and interval; the k-matched minus bound-matched difference.
Outcomes: R7.1 or R7.2 failing means the gradient depends on the wiring class and F5 must be narrowed;
R7.3 failing means plug-in bias across histories contributes to the across-target inversion.
