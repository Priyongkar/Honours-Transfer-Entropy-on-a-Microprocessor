# R7 results: the locking gradient within one wiring class; the matched-bound inversion matched on history

Design `PREDICTIONS_R7.md` (sha256 be66a622f7390222…). Run of record `output/R7_locking_within_class_run_of_record.log`
(12:21, 5 Nov 2026); a rerun with one descriptive block added (`output/R7_locking_within_class.log`, 12:22) left
every designed statistic unchanged. Results `output/R7_locking_within_class_results.json`.

(a) Locking gradient within a class (within-target AUC of lag-1 TE; mixed / locked / strict targets; targets with a
wired source of the class in brackets):
- Pull-down pairs (783 wired): 0.551 [0.395, 0.697] (50) / 0.395 (41) / 0.441 [0.389, 0.480] (165);
  mixed minus strict +0.110 [−0.056, +0.263]. **R7.1 not met**: same direction, interval includes 0.
- Pass-only pairs (389 wired): 0.628 [0.598, 0.667] (39) / 0.593 (41) / 0.485 [0.391, 0.590] (102);
  mixed minus strict +0.142 [+0.027, +0.240]. **R7.2 met.**

(b) Matched-bound inversion matched on the target's selected history k:
- Designed scope, strictly locked same-phase pairs (77,006 pairs, 315 wired): bound deciles 0.452 [0.419, 0.487];
  deciles x k 0.436 [0.407, 0.469]; difference −0.015 [−0.052, +0.021]. **R7.3 met** (upper limit below 0.5).
- A3's own scope, all same-phase pairs of locked units (770 wired), added after the first run because the design
  file described A3's scope as strict by mistake (descriptive): bound deciles 0.400 [0.375, 0.429] (A3: 0.400);
  deciles x k 0.465 [0.438, 0.496]; difference +0.065 [+0.033, +0.098].

Reading: the locking gradient is not a wiring-class artefact for pass-only pairs; for pull-down pairs it points the
same way but few targets active in both phases have pull-down wired sources, so it is not established there.
Among strictly locked pairs the inversion owes nothing to history-dependent bias; among all same-phase locked
pairs, part of it does (the AUC rises by 0.065 when histories are matched) and the rest remains below chance.
