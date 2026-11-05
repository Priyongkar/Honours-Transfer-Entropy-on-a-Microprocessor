# 12 Zero-lag partial correlation as a comparator (R6)

The benchmarks the thesis cites found zero-lag partial correlation the best detector of connections. R6 computes it,
plain and given the clock phase, on the clock-aligned Donkey Kong recording, with predictions hash-logged first.
- `design/predictions_R6.md`: three predictions (sha256 23767735…). R6.1 and R6.3 met; R6.2 (given the phase above
  plain) not met: conditioning on the phase lowers it.
- `zero_lag_partial_correlation_comparator_R6.py`: the test. Run of record `R6_partial_correlation_run_of_record.*`;
  a rerun added descriptive quantities and changed no designed statistic.
- `RESULTS_R6.md`: the numbers. Partial correlation is the best activity-based score (group AUC 0.790; precision 23.2 %
  at the multivariate allocation), still below layout distance and source size; weak at the physical level.
Outputs: `results/12_partial_correlation_R6/`.
