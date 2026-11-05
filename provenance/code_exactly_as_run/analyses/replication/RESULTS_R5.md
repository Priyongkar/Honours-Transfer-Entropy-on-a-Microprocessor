# R5 results: a second, independent draw of the H6 design

Design `PREDICTIONS_R5.md` (sha256 8b79f8c8114e3581…, saved 12:17 on 5 Nov 2026); search `run_h6_second_draw.py`
12:19–14:04 (`output/R5_search.log`; targets in `output/R5_targets.json`, seed 2027, none in R1); scored afterwards by
`supplementary_recompute/score_h6_second_draw.py` (`output/h6_second_draw_R5.json`).

R5 (100 targets, 98 kept by the omnibus test): 200 sources, 35 wired; group precision 17.50 % [11.33, 24.21] against
3.29 % for the main run's lag-1 network on the same targets; +14.21 [+7.68, +21.43] (two-way [+4.32, +26.20]).
132 of 204 selected sources at lag 0 (33 wired); 72 at lag 1 (2 wired). Lag-0 CMI at the same allocation 15.50 %
[11.17, 20.86]; search minus lag-0 CMI +2.00 [−3.11, +7.41]. Source-group size at the same allocation 30.00 %
[24.39, 35.39]; size minus search +12.50 [+3.30, +21.05]. Physical: search minus lag-1 +0.37 [−0.04, +1.94].
Pooled with R1 (200 targets): 17.35 % [13.10, 21.66] against 3.73; +13.62 [+9.46, +18.16] (two-way [+7.13, +20.38]);
minus lag-0 CMI +0.91 [−2.84, +4.48]; size minus search +10.05 [+3.87, +16.01]; physical +0.41 [+0.08, +1.18].

Predictions: R5.1 met; R5.2 met (64.7 % at lag 0); R5.3 met; R5.4 (descriptive) interval includes 0.
Reading: the H6 result replicates on targets it was not drawn from; the physical agreement check includes zero on the
second sample alone and is met pooled, so the verdict stays fragile at the physical level.
