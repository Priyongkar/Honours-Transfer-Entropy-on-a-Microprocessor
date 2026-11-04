# R2: replication on two further programs, fixed before any statistic is computed

Written on 4 November 2026 (afternoon), after the recordings were generated and before any statistic
of them was computed. The scripts that run these tests print the SHA-256 of this file, of themselves
and of the shared module at the start of their logs.

## What has and has not been seen

- **Generated, not analysed.** `data/regen/pitfall.regen.ts.spikes.npy` and
  `data/regen/spaceinvaders.regen.ts.spikes.npy`, 639,999 rows × 3,510 transistors each, produced by
  the authors' simulator (Sim2600, ported to Python 3; `analyses/pitfall/record_game.py`). The same code
  regenerates the published Donkey Kong recording identically in every analysed sample
  (`analyses/pitfall/CONVERSION_NOTE.md`). Of the two new recordings only the file sizes and the
  generation logs (row counts, run times) have been seen.
- **Seen, from a different source.** A Pitfall recording from Wang et al. (2023), produced by a different
  simulator and only 340 phases long, showed 51 of 286 groups phase-locked (Donkey Kong: 445 of 519).
  That makes the phase-structure prediction (R2.3) genuinely uncertain; it is stated here because it was
  seen.
- **Status.** Hash-logged and **held out**: no prediction below was formed from either recording. The
  predictions are those of the Donkey Kong results and of Propositions 2 and 3, which were derived from
  the Donkey Kong results.

## Pipeline (identical to the main run unless stated)

Per game, separately:
1. **Reduction.** The main run's rule: crop to the last nonzero row; keep every third row at the offset
   with the most events (`row_stride 3`, `phase_rule max_total_events_first_tie`).
2. **Clock alignment** (the rule applied to Donkey Kong, fixed here before it is applied): find the
   transistors gated by `cp1` and `cclk`; the analysis starts at the first reduced sample after the last
   raw clock event that falls off the kept offset. All statistics are on these clock-aligned samples;
   the whole-recording values are reported beside them.
3. **Groups.** Identical non-constant reduced series (`group_identical`).
4. **Histories.** Per target, the k from 1 to 12 maximising AIS minus its analytic null mean.
5. **Scores.** Lag-1 TE and lag-0 CMI (the thesis's estimator, `te_all_sources`), zero-lag MI.
6. **References.** The terminal-to-gate netlist at the group and physical levels (same chip, same
   netlist); latch group pairs and phase classes as in `recording.py`; layout distance and source-group
   size as in the main run.
7. **Inference.** Target-cluster bootstrap, 1,000 draws, seed 2026. Stratified and within-target AUCs as
   in A2, A9 and A12. Decision rule as for A1–A19: an AUC is supported if its lower limit exceeds 0.5 and
   contradicted if its upper limit is below 0.5; a paired difference is supported if its lower limit
   exceeds 0.

No multivariate (IDTxl) analysis is run.

## Predictions (each game separately; primary in bold)

| ID | Prediction | Statistic | Met if |
|---|---|---|---|
| **R2.1** | Lag-1 TE does not rank the wiring above chance, and ranks it below chance | Group-level AUC against the netlist, clock-aligned | Upper limit < 0.5 |
| R2.1b | The same at the physical level (H1's deciding level) | Physical AUC | Lower limit ≤ 0.5 |
| **R2.2** | Lag-0 CMI ranks the wiring better than lag-1 TE (H2's comparative clause) | Paired group-AUC difference, lag-0 minus lag-1 | Lower limit > 0 |
| **R2.3** | Most activity groups are phase-locked | Share of groups with ≥ 95 % of events in one parity | > 0.5 (point estimate; descriptive count reported) |
| **R2.4** | Within strictly locked same-phase targets, targets follow wired sources more than unwired ones, and lag-1 TE ranks wired sources below where their activity places them (A9.1, A9.2) | Within-target followership AUC; lag-1 TE minus activity | AUC lower limit > 0.5; difference upper limit < 0 |
| R2.5 | The penalty is graded by the target's locking (A12.1) | Within-target lag-1 TE AUC, mixed minus strict targets, clock group excluded | Lower limit > 0 |
| **R2.6** | Lag-1 TE identifies latch pairs above matched controls (A2) | Stratified AUC, oriented latch pairs vs pairs alike in phase and rate | Lower limit > 0.5 |
| R2.6b | Lag-0 CMI does not (A2) | The same for lag-0 CMI | Upper limit < lag-1 TE's lower limit |
| R2.7 | Layout distance ranks the wiring above every information measure | Group AUC of distance vs lag-0 CMI (paired) | Lower limit > 0 |

**Across programs (descriptive, three points; no test).** Donkey Kong, Pitfall and Space Invaders
ordered by the share of strictly locked same-phase wired group pairs, beside their lag-1 group AUCs. The
account predicts that the more locked program has the lower AUC.

## What would count against the account

R2.1 failing with a lag-1 AUC above chance on a locked program; R2.4 failing (wired sources not followed
more, or not penalised); R2.6 failing (lag-1 TE not finding latch transfer). R2.3 failing (most groups
mixed) would not count against the account but would change what it predicts: R2.5's mixed bin would
dominate and R2.1 would be expected to move towards or above chance.

## Reporting

Every prediction, met or not, goes into the thesis's scoreboard and the failed-predictions table, with
both games in one table. Whole-recording values are reported beside the clock-aligned ones. The
replication tests dependence on the program, not on the system: it is the same chip, netlist and
simulator.
