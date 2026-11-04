# 09 Recomputing the supplementary results

The original outputs of the supplementary analyses (baselines, parity conditioning, the RQ2 decomposition,
enrichment) were lost. These two scripts recompute them from the main run's saved score matrices and
networks. They are not hash-logged: the target values were already printed in the thesis.

- `recompute_supplementary_results_pass1.py`: blocks B (baselines, within-covariate deciles, covariate model), P (parity), D (RQ2 networks and increments, target and two-way bootstrap), E (enrichment). `--selftest` reproduces the main run's AUCs.
- `recompute_supplementary_results_pass2_definitions_matched.py`: items where pass 1's definitions differed from the thesis (phase class, physical-fitted covariate model, rates model, self-information and entropy deciles, joint Holm).
- `COMPARISON_WITH_THESIS.md`: every value against the printed one.

Outputs: `results/09_recomputed_supplementary/`.

Added 4 November (evening):
- `recompute_supplementary_results_pass3_verdict_statistics_H1_to_H4.py`: H1-H4 point values, target-cluster, two-way and vertex intervals and Holm p-values; H3 without the clock; H4 within parity strata; main H2 robustness rows. `COMPARISON_PASS3.md` compares each with the thesis (point values and target-cluster intervals reproduce exactly; re-drawn two-way and vertex limits within 0.015).
- `h1_robustness_checks_top1_bias_exact_p.py`: each target's strongest source wired (lag-1 TE 4.8 %, lag-0 CMI 21.0 %), bias-subtracted and exact-p group AUCs.
- `h6_run_of_record_R1_networks_and_intervals.py`, `h6_run_of_record_R1_lag_split_baselines_classes.py`: the H6 table on R1, the run of record (the original H6 target list was not saved), from R1's saved per-target IDTxl outputs.
- `results/figures/scripts/fig_followership_and_locking_three_programs.py`: the thesis figure of the followership penalty on five data sets, its gradient with locking on three programs, and the synthetic networks.
