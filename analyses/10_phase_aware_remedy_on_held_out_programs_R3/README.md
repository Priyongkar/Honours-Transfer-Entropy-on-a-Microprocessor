# 10 The phase-aware remedy on the held-out programs (R3)

Zero-lag mutual information conditioned on the phase was formed on Donkey Kong after the main results
(thesis Section 5.7.2). R3 tests it on Pitfall and Space Invaders, whose other results had been seen (R2)
but on which this score had not been computed. Predictions were written and hash-logged first.

- `design/predictions_R3.md`: the five predictions and the decision rule (sha256 282eeba2…, printed at
  the start of both logs).
- `remedy_on_held_out_programs_R3.py`: the R2 pipeline (reduction, clock alignment, groups, reference,
  bootstrap draws) with the remedy, unconditioned zero-lag MI, the phase relation alone and source-group
  size added.
- `RESULTS_R3.md`: the results. R3.1, R3.2, R3.3 and R3.5 met on both programs; R3.4 (remedy above the
  phase relation alone) met on neither.
- Post hoc, after R3.4 failed, descriptive only:
  - `posthoc_remedy_against_phase_class_and_ap_lift.py`: AP lifts on all three programs and the Donkey
    Kong clock-aligned difference from the phase relation.
  - `posthoc_clock_aligned_start_1407_against_1408.py`: the Donkey Kong clock-aligned scores from sample
    1,407 (used in analysis 09) against 1,408 (the alignment rule). Only the two scores without a target
    history move; the thesis prints the values from 1,408.

Outputs: `results/10_remedy_held_out_R3/`.
