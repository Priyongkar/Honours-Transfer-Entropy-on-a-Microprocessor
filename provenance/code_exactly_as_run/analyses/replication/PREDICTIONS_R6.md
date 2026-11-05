# R6: zero-lag partial correlation as a comparator on the recording, fixed before it is computed

Written on 5 November 2026, before partial correlation was computed on any of the recordings.
`run_partial_correlation.py` prints the SHA-256 of this file and of itself at the start of its log.

## Why

The benchmarks the thesis cites (simulated fMRI; simulated calcium imaging) found zero-lag partial
correlation or inverse covariance the best detector of connections. The thesis's contemporaneous scores
are information-theoretic (zero-lag MI, lag-0 CMI, MI given the phase). Partial correlation conditions
each pair on every other unit's simultaneous sample, so it is the natural contemporaneous comparator
and the one an examiner expects.

## What has been seen

Everything printed in the thesis, including the clock-aligned Donkey Kong scores from sample 1,408
(zero-lag MI 0.516, MI given the phase 0.641, phase relation alone 0.647, layout distance 0.928,
source-group size 0.864, lag-0 CMI 0.561). Partial correlation has not been computed on any recording.

## Design

Donkey Kong, reduced recording, clock-aligned samples (from 1,408), 519 groups. Scores, all symmetric
and applied to every ordered pair:
- `partial`: |partial correlation| from the pseudo-inverse of the 519 x 519 covariance of the binary
  change series.
- `partial|parity`: the same after removing each unit's parity-specific mean (conditioning on the phase).
- `correlation`: |Pearson correlation| (descriptive).
Comparators recomputed on the same samples: zero-lag MI; zero-lag MI given the phase; layout distance.
Group-level AUC against the wiring, 95% target-cluster intervals (1,000 draws, seed 2026), paired
differences from the same draws; physical and within-target AUCs and AP lift reported.

## Predictions

- **R6.1** `partial` minus zero-lag MI: lower limit above 0 (conditioning on the other units removes
  indirect co-activity that MI keeps).
- **R6.2** `partial|parity` minus `partial`: lower limit above 0 (the phase is the shared driver).
- **R6.3** Layout distance minus `partial|parity`: lower limit above 0 (F4 holds against the best
  contemporaneous comparator).
- **R6.4 (descriptive)** `partial|parity` against MI given the phase: no prediction; reported.

A prediction is met if the interval condition holds. If R6.3 fails, F4's claim that activity-free
baselines beat every information measure would need to be narrowed to information-theoretic measures.
