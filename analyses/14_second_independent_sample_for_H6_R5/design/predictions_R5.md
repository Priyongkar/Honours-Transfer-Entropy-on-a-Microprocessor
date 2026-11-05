# R5: a second, independent draw of the H6 design, fixed before the run

Written on 5 November 2026, before any target of this draw was analysed. `run_h6_second_draw.py` prints
the SHA-256 of this file, of itself and of the shared modules at the start of its log.

## What has been seen

The H6 result rests on one stratified sample of 100 targets (R1, the run of record: group precision
17.2 per cent against 4.1 for the main run's lag-1 network on the same targets, +13.13 [+7.16, +19.23];
equal to lag-0 CMI at the same allocation; below the source-size baseline). The thesis names the single
sample as an external-validity threat. Nothing about any other target's lag-0 search has been seen.

## Design (identical to R1 unless stated)

- Targets: a fresh stratified draw of 100 by `idtxl_runs.draw_targets` (same bands of selected history,
  same quotas), seed 2027, **excluding the 100 targets of R1**, so the two samples do not overlap.
- Search: IDTxl multivariate TE, source lags 0 and 1, target history up to 12, 200 surrogates, alpha 0.05,
  the analytic-bits fix, every other setting as R1; the same data matrix (reduced recording, 519 groups).
- Scoring: as the run of record (`h6_run_of_record_R1.py`): omnibus network FDR; group precision of the
  selected edges on the R5 targets; the main run's lag-1 network on the same targets; lag-0 CMI allocated
  as the search allocated; source-group size allocated the same way (added); target-cluster and two-way
  bootstraps, 1,000 draws, seed 2026/2027.

## Predictions

- **R5.1 (primary, the H6 statistic)** Group precision of the lags-{0,1} search minus the main run's lag-1
  network on the same targets: lower limit of the target-cluster interval above 0.
- **R5.2** More than half of the selected sources are at lag 0.
- **R5.3** Source-group size at the same allocation minus the search's precision: lower limit above 0
  (the baseline stays ahead).
- **R5.4 (descriptive)** The search minus lag-0 CMI at the same allocation: interval includes 0. A lower
  limit above 0 would count against the reading that conditioning adds nothing detectable once lag 0 is
  admitted; this is not a test of equivalence.

A prediction is met if the stated interval condition holds. Also reported: R5's precision beside R1's,
and R1 and R5 pooled (200 targets).
