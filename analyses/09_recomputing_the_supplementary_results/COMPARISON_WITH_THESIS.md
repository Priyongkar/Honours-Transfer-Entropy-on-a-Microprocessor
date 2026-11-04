# Recomputed supplementary results against the values printed in the thesis (v10)

The original output files for these results (supp_r3_baselines, supp_parity_*, supp_rq2_decomposition_*,
supp_h6_*) are not on any machine or backup. They were recomputed from the main run's saved score matrices
and networks by `recompute_supplementary.py` (pass 1) and `recompute_supplementary_extras.py` (pass 2).
Both scripts use the target-cluster bootstrap (1,000 draws, seed 2026). The two-way bootstrap uses a
second, independent multinomial for sources (seed 2027).

## Reproduced (equal at the precision printed, or within 0.001)

- Baselines: source size, layout distance, size product and the clock-excluded rows (AUCs, AP lifts and allocation-matched precisions)
- Panel B within-covariate deciles, including source self-information (0.377/0.460/0.586/0.557 vs printed 0.378/0.460/0.587/0.557) and source entropy (within 0.001)
- Parity-conditioned scores, the reversed lag-1 row and the strata; clock-aligned lag-0 CMI|parity 0.620 / within-target 0.585 (printed 0.619 / 0.585)
- Phase class only 0.647, AP lift 1.42, within-target 0.644; phase class then lag-0 0.700
- RQ2 networks at the allocation: 56/20/29/29 wired; precisions 3.88/1.39/2.01/2.01/13.93/13.72/19.47/25.1/20.0 %
- Increments (target bootstrap): total +2.49 [1.62, 3.36] (printed [1.62, 3.28]); conditioning +1.87 [1.06, 2.66] (printed [1.05, 2.63])
- Enrichment: two-step 2/8 and 0/6; shared parent 63/47 and 17/1, ratio 1.38, Holm 0.001; Holm across all eight classes 0.188 and 0.083 (printed 0.19 and 0.083)

## Differ

| Item | Printed | Recomputed | Reading |
|---|---|---|---|
| Covariate model, group-fitted AUC; lag-0 increment | 0.974; +0.0043 [0.0028, 0.0059] | 0.969; +0.0041 [0.0025, 0.0060] | The feature specification differs slightly. The conclusion is unchanged: a small positive increment. |
| Covariate model, physical-fitted AUC; increment | 0.894 [0.877, 0.911]; -0.0007 [-0.0049, +0.0033] | 0.890 [0.872, 0.908]; +0.005 [-0.002, +0.012] | The conclusion is unchanged: the interval includes zero. |
| Rates model (group AUC, AP lift, within-target) | 0.623, 10.39, 0.715 | 0.625, 8.75, 0.721 | The specification differs. Not used in any verdict. |
| Two-way interval, conditioning increment | [-0.06, +4.48], 97.3 % positive | [+0.09, +4.41], 98.2 % positive | Same scheme with different draws. The lower bound sits on zero, and its sign depends on the draws. |
| Two-way interval, total | [+0.57, +5.20], 99.8 % | [+0.50, +5.05] | Same conclusion. |

## Consequence for the text

rq2.tex currently says the conditioning increment is "not established under resampling of sources". That
sentence rests on the sign of a bound that changes with the draws. The suggested wording: *its two-way
interval reaches zero (lower bound between -0.06 and +0.09 across two independent sets of 1,000 draws;
97-98 % of draws positive), so the conditioning share is marginal under resampling of sources.* The H5
verdict is unaffected.

Not recomputable: the original H6 run (236 edges, 17.37 %), because its target list was not saved. It
rests on R1 (17.2 %) and A17 (17.3 %), with a disclosure.
