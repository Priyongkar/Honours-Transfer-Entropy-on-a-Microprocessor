# R3 results: the phase-aware remedy on the held-out programs

Design: `PREDICTIONS_R3.md` (sha256 282eeba2395a0706a65bac07cf1ba4e449b53e4140fa40ef5ab710bbbd8226f3, printed
at the start of both logs). Run on 5 November 2026; outputs `output/R3_Pitfall_results.json`,
`output/R3_SpaceInvaders_results.json`.

Clock-aligned, group level, 95% target-cluster intervals (1,000 draws).

| | Pitfall | Space Invaders | Met |
|---|---|---|---|
| R3.1 MI\|parity group AUC above 0.5 | 0.648 [0.631, 0.666] | 0.661 [0.644, 0.679] | yes ×2 |
| R3.2 MI\|parity minus lag-1 TE | +0.252 [+0.227, +0.276] | +0.262 [+0.239, +0.288] | yes ×2 |
| R3.3 MI\|parity minus lag-0 CMI | +0.076 [+0.063, +0.090] | +0.075 [+0.061, +0.090] | yes ×2 |
| R3.4 MI\|parity minus phase class | +0.000 [−0.019, +0.018] | +0.014 [−0.004, +0.033] | no ×2 |
| R3.5 Layout distance minus MI\|parity | +0.281 [+0.259, +0.304] | +0.268 [+0.246, +0.288] | yes ×2 |

Reported, not predicted: whole-recording MI|parity 0.735 and 0.745; physical AUC 0.582 and 0.587;
within-target AUC 0.682 and 0.692.

## Post hoc, after R3.4 failed (descriptive)

- `../supplementary_recompute/posthoc_r3_phase_class.py`: AP lift of MI|parity against phase class only,
  clock-aligned: Pitfall 4.90 vs 1.42, Space Invaders 5.74 vs 1.42, Donkey Kong 4.31 vs 1.42. Donkey Kong
  MI|parity minus phase class, from sample 1,408: −0.005 [−0.023, +0.013] (200 draws).
- `../supplementary_recompute/posthoc_align_offset.py`: the Donkey Kong clock-aligned scores from sample
  1,407 (used by `recompute_supplementary.py`) against 1,408 (the alignment rule used here). Scores with a
  target history agree to 0.002; zero-lag MI falls from 0.524 to 0.516 and MI|parity from 0.661 to 0.641.
  The thesis now prints the values from 1,408.

## Reading

The remedy transfers to programs it was not formed on: above chance, above lag-1 TE and above
unconditioned lag-0 CMI, below layout distance. In AUC it does no better than the phase relation alone;
its advantage is at the top of each ranking (AP lift three to four times the phase relation's). The
whole-recording values are raised by the samples before the reduction locks to the clock.
