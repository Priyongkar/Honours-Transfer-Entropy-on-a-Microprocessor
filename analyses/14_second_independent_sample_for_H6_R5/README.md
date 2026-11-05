# 14 A second, independent sample for the H6 design (R5)

H6 (admitting lag-0 sources raises the multivariate search's precision about fourfold) rested on one stratified sample
of 100 targets. R5 draws a second sample of 100 with a new seed, excluding every target of the first, with predictions
hash-logged before the search ran, and scores it as the run of record.
- `design/predictions_R5.md` (sha256 8b79f8c8…): four predictions; R5.1–R5.3 met, R5.4 descriptive (interval includes 0).
- `h6_design_second_independent_sample_R5.py`: the search (IDTxl, lags 0 and 1, every other setting as R1); it reads
  the saved reference from the data folder and writes per-target outputs there; copies of those outputs are in
  `results/14_h6_second_sample_R5/R5_targets/`.
- `score_second_sample_against_run_of_record_R5.py`: precision, recall, allocation-matched comparators (lag-0 CMI and
  source-group size), target and two-way intervals, for R5, R1 and both pooled.
- `RESULTS_R5.md`: 17.5 % against 3.3 % on the same targets (+14.21 [+7.68, +21.43]); pooled 17.4 % against 3.7 %.
Outputs: `results/14_h6_second_sample_R5/`.
