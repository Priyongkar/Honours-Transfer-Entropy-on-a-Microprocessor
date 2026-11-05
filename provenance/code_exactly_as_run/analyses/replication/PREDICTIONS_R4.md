# R4: physical delay against phase relation on the simulated chip, fixed before any statistic is computed

Written on 5 November 2026, before any of the statistics below was computed. `run_delay_vs_phase.py`
prints the SHA-256 of this file, of itself and of the shared modules at the start of its log.

## The question

Lag-1 TE finds the latch pairs (one-step physical delay, units locked to opposite phases) and misses the
wiring (immediate physical influence, mostly between units locked to the same phase). Two accounts fit:

- **Phase account.** For units locked to the clock, the phase relation fixes the delay *in samples*: a
  target locked to the other phase can show a change only one sample later, whatever the physical delay;
  a target locked to the same phase only zero or two samples later. Lag-1 TE then follows the phase
  relation, not the physical delay.
- **Physical-delay account.** Lag-1 TE follows the physical delay: it finds one-step influence and misses
  immediate influence, whatever the phase relation.

The two coincide for the latch pairs and for most wiring, which is why the thesis cannot yet separate
them. They come apart for immediate influence between units locked to opposite phases, and for one-step
influence between units locked to the same phase. Clamps measure the physical delay; locking is read from
the activity; so the two can be crossed.

## What has been seen

- A16 (`A13_A16_chip_results.json`): class sizes of immediate, delayed-only and no influence, and that
  lag-1 TE ranks delayed-only influence above none at 0.589 and 0.537 (programs 1 and 2), not above lag 0.
  Those AUCs pooled all phase relations, stratifying by phase class. Nothing below has been computed.
- On the recording (thesis Table A.6), the 30 wired group pairs of opposite parity rank at 0.752 under
  lag 1. That result informs R4.1 and R4.2, which are therefore partly anticipated.

## Data and classes (identical to A16 unless stated)

Both simulated programs; change coding; groups of identical change series; node pairs in different
groups. Physical delay from the saved timed clamps (`output/A16_program{1,2}_timed_clamps.npz`, 40 moments
per clock level): **immediate** = the clamp changes the target in the same half-step at some moment;
**one-step** = it changes the target only after one half-step (delayed-only); **none** = no change in
either direction at either delay. Phase relation from the change series (locked: at least 95 per cent of
changes in one phase): **same** (both locked, same phase), **opposite** (both locked, different phases),
**mixed** (either unit not locked).

Each cell's positives (immediate or one-step pairs with that phase relation) are compared with the none
pairs of the same phase relation, stratified by the source's and target's event-rate quintiles; AUC with
95% target-cluster intervals (1,000 draws, seed 2026). Scores: lag-1 TE and lag-0 CMI with each group's
selected history (k up to 6), as in A16; lag-2 TE reported.

## Predictions (phase account)

- **R4.1** Opposite phase, immediate influence: lag-1 TE ranks it above none (lower limit > 0.5).
- **R4.2** Opposite phase, immediate influence: lag-1 TE minus lag-0 CMI > 0.
- **R4.3** Lag-1 TE AUC for opposite-phase immediate influence minus that for same-phase one-step
  influence > 0 (the phase relation, not the physical delay, decides what lag 1 finds).
- **R4.4** Mixed pairs: lag-1 TE AUC for one-step influence minus that for immediate influence > 0 (where
  the phase does not fix the sample delay, the physical delay decides).

A prediction is met on a program if the interval excludes the null in the predicted direction; met ×2 if
on both programs. A cell with fewer than 10 positive pairs is reported descriptively and decides nothing.

## Reported, not predicted

Cell sizes; lag-0 CMI and lag-2 TE in every cell; lag-1 TE for same-phase immediate and opposite-phase
one-step influence.

## What the outcomes would mean

R4.1–R4.3 met: the thesis's "delay matches lag" holds for delay *in samples*, which under locking the
phase relation sets; the latch control and the wiring failure differ in phase relation and that is the
operative difference. R4.1–R4.3 failed with the physical-delay pattern instead: lag 1 follows the
physical delay, and the phase relation is not the cause. R4.4 tests whether physical delay matters where
locking does not fix the sample delay.
