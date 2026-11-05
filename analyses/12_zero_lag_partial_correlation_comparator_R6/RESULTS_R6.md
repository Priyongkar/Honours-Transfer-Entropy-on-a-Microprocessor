# R6 results: zero-lag partial correlation on the clock-aligned Donkey Kong recording

Design `PREDICTIONS_R6.md` (sha256 23767735…). Run of record `output/R6_partial_correlation_run_of_record.{log,json}`
(12:19–12:24, 5 Nov 2026); a rerun at 12:26 added descriptive quantities (physical AP lift, precision at the
multivariate network's allocation, the paired difference from plain correlation) and left every designed
statistic unchanged (`output/R6_partial_correlation_results.json`).

Group AUC [95% target-cluster], AP lift, physical AUC, within-target AUC, precision at allocation (%):
- partial correlation: 0.790 [0.772, 0.808], 10.83, 0.659, 0.854, 23.2 [21.3, 24.9]
- partial correlation given the phase: 0.721 [0.704, 0.739], 3.93, 0.660, 0.723, 9.9
- correlation: 0.669 [0.642, 0.696], 6.82, 0.634, 0.772, 11.3
- zero-lag MI: 0.516; MI given the phase: 0.641; layout distance: 0.928 (all as printed)

Predictions: R6.1 met (+0.273 [+0.247, +0.300] over zero-lag MI); R6.2 **not met** (given the phase minus plain
−0.069 [−0.080, −0.058]); R6.3 met (distance minus partial|phase +0.207 [+0.188, +0.226]).
Descriptive: partial minus plain correlation +0.121 [+0.096, +0.145]; partial|phase minus MI|phase +0.079 [+0.060, +0.098].

Reading: partial correlation is the best activity-based score in the thesis at the group level, consistent with the
fMRI and calcium benchmarks, and still below both activity-free baselines; physically it is weak (AP lift 1.51,
precision at allocation 0.10 %). Conditioning on the phase lowers it: a reading made afterwards is that conditioning
on every other unit, the clock group among them, already removes the phase, so removing it again discards
dependence that carries wiring.
