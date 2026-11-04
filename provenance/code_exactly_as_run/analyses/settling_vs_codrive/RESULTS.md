# Settling versus co-drive: results

Design and decision rule: `PREDICTIONS.md` (SHA-256 e37c0c32…f8a0e, printed
at the start of both runs). Script: `settling_codrive.py`. Outputs:
`output/results.json`, `output/auc_table.csv`, logs in `output/`. Run on
2 November 2026 from the saved chip outputs of the final run; nothing was
re-simulated. Full run about 5 minutes.

## Verdict under the pre-specified rule

**Direct coupling leaves lag-0 dependence beyond matched co-drive, and the
effect is material.** Under change coding, lag-0 CMI separates directly
influencing wired pairs from co-driven, non-influencing, unwired pairs
matched on phase class, change rates and co-drive strength with a stratified
AUC of **0.620 [0.579, 0.658]** (1,328 direct pairs; 285,854 co-driven pairs
in 232 shared strata). The lower limit is above 0.5 and the point estimate
is at least 0.60, the threshold set in advance. The margin over 0.60 is
small.

## Table

Stratified AUC, direct versus comparison class, 95% target-cluster bootstrap
(1,000 draws). P1: against co-driven-only pairs (pre-specified; primary is
change coding, lag-0 CMI). P3: against wired-but-idle pairs (secondary).
Same-phase rows are post hoc: only pairs whose members are locked to the
same phase.

| Representation | Score | P1, all pairs | P1, same phase (post hoc) | P3, all pairs |
|---|---|---|---|---|
| Change | lag-0 CMI | **0.620 [0.579, 0.658]** | 0.759 [0.717, 0.797] | 0.642 [0.562, 0.715] |
| Change | lag-0 CMI given parity | 0.639 [0.594, 0.680] | 0.788 [0.744, 0.825] | 0.646 [0.565, 0.720] |
| Change | lag-1 TE | 0.310 [0.278, 0.344] | 0.290 [0.242, 0.342] | 0.580 [0.487, 0.668] |
| State | lag-0 CMI | 0.826 [0.799, 0.848] | 0.903 [0.886, 0.918] | 0.744 [0.677, 0.801] |
| State | lag-0 CMI given parity | 0.772 [0.742, 0.800] | 0.904 [0.887, 0.919] | 0.739 [0.679, 0.791] |
| State | lag-1 TE | 0.322 [0.297, 0.350] | 0.273 [0.237, 0.311] | 0.542 [0.431, 0.645] |

Class means of lag-0 CMI (millibits), change coding: direct 53.3, wired idle
15.8, co-driven only 22.3, other pairs 12.9. P2, Spearman correlation of
lag-0 CMI with co-drive strength among co-driven-only pairs: 0.14 (change),
0.12 (state).

## Reading

1. **Co-drive is real but smaller.** Co-driven pairs score above other
   pairs (22.3 against 12.9 millibits), and more strongly co-driven pairs
   score higher (P2 positive). So co-drive does produce contemporaneous
   dependence, as Section 6.4.2 of the thesis feared.
2. **Direct coupling adds to it.** Matched on co-drive strength, phase and
   rates, direct pairs still rank above co-driven ones. The effect is
   largest exactly where the recording is ambiguous, among same-phase pairs
   (0.759 under change coding, 0.903 under states).
3. **Change coding costs about half the signal.** The stratified AUC falls
   from 0.826 (states) to 0.620 (change), consistent with the thesis's
   Section 5.4.2.
4. **The wire must carry influence.** Wired pairs that never transmit an
   effect at the sampled moments score like co-driven pairs (15.8 against
   22.3 millibits), and direct pairs beat them (P3 0.642). Lag-0 CMI tracks
   influence through the wire, not the wire itself.
5. **Lag 1 is anti-informative here too.** At lag 1, direct pairs rank
   below co-driven pairs (0.31), the aliasing pattern of Section 6.3.2 on a
   second program.

The expectation written in `PREDICTIONS.md` (above 0.5 under states, closer
to 0.5 under change coding) held for the direction; under change coding the
AUC stayed clearly above 0.5.

## Limits

- Simulated chip and a different, short program; the recording's pairs are
  not tested. The result shows that lag-0 CMI *can* separate settling from
  co-drive under change coding on this netlist, not how much of the
  recording's H2 effect is settling.
- E0 samples 40 moments. A pair coded "no influence" may transmit at
  unsampled moments, which would put true influence among the co-driven
  pairs and bias P1 towards 0.5, so P1 is conservative in that respect.
- Almost every unwired pair has a common cause (the clock nodes reach
  nearly everything): 902,950 of 1.2 million pairs are co-driven only. The
  co-drive-strength strata, not exclusion, handle that.
- 30% of direct pairs under change coding (582 of 1,910) fall within one
  activity group and cannot be scored, as on the recording.
- The same-phase restriction was chosen after the primary result was seen.

## Suggested thesis changes

- **§6.4.4 and Table 4.2 (H2).** The causal clause stays undecidable *on the
  recording*, but replace "not decidable" with: "On the simulated chip,
  where an interventional reference separates the two, lag-0 CMI ranks
  directly influencing wired pairs above co-driven, non-influencing pairs
  matched on phase, rates and co-drive strength (stratified AUC 0.620
  [0.579, 0.658] under change coding; 0.759 among same-phase pairs, post
  hoc). Within-phase coupling therefore contributes beyond co-drive under
  the recording's coding, on this netlist and a different program." Keep
  the verdict at "partially supported": this is supplementary evidence and
  cannot change a verdict under Section 4.1.2.
- **§6.4.2.** Soften "no analysis that separates them was run" and point to
  the new section.
- **Table 8.1, mechanism 1, Standing.** "Wiring-related dependence is
  contemporaneous; on the simulated chip direct coupling contributes beyond
  matched co-drive; not separated on the recording."
- **§8.4 future work item 2.** Remove the part now done; keep sub-phase
  timing and running the recorded program.
- **New subsection 5.4.5** (Methods + result, about a page) with the table
  above, the pre-registration hash and the limits. Mark it supplementary but
  pre-specified in its own file, dated before it ran.
