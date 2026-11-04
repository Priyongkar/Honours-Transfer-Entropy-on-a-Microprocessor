# 09 Recomputing the supplementary results

The original outputs of the supplementary analyses (baselines, parity conditioning, the RQ2 decomposition,
enrichment) were lost. These two scripts recompute them from the main run's saved score matrices and
networks. They are not hash-logged: the target values were already printed in the thesis.

- `recompute_supplementary_results_pass1.py`: blocks B (baselines, within-covariate deciles, covariate model), P (parity), D (RQ2 networks and increments, target and two-way bootstrap), E (enrichment). `--selftest` reproduces the main run's AUCs.
- `recompute_supplementary_results_pass2_definitions_matched.py`: items where pass 1's definitions differed from the thesis (phase class, physical-fitted covariate model, rates model, self-information and entropy deciles, joint Holm).
- `COMPARISON_WITH_THESIS.md`: every value against the printed one.

Outputs: `results/09_recomputed_supplementary/`.
