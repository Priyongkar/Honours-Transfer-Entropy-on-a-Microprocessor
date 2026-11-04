# A7 — second simulated program: design and predictions

Written on 3 November 2026 (about 01:00) before the second program was
simulated. The first program's and the recording's results for A1, A1b, A2,
A5 and A6 were known when this was written (see `PREDICTIONS.md` for those
tests); this file fixes what must replicate on a program not used to form
them. The script prints this file's SHA-256 first.

## Program

A 6502 program written for this test, exercising what the Fibonacci/LFSR
program does not: SBC, the stack (TXS, PHA/PLA, JSR/RTS), (zp),Y indirect
addressing, indexed compare-and-swap loops. Each iteration fills a 16-byte
table with a linear congruential sequence, bubble-sorts it, calls a checksum
subroutine (subtraction through an indirect pointer, with stack use), and
writes the checksum and the table maximum to memory. Accepted only if every
byte written matches an independent 6502 emulator (py65) running the same
binary, and a deliberately altered byte fails the check.

Simulation: same switch-level simulator and netlist; reset, warm-up 200
half-steps, 60,000 recorded half-steps; clamp interventions E0/E1 at 40
moments (seed 2) on every non-constant node, as for the first program.

## Analyses and predictions (change coding primary; same definitions,
## strata, bootstrap and decision rule as `PREDICTIONS.md`)

| ID | Test | Prediction |
|---|---|---|
| A7-C2 | Lag-1 TE, wired vs unwired node pairs (stratified AUC) | Interval includes or lies below 0.5 |
| A7-A1b | Lag-1 TE, oriented latch pairs vs no-influence pairs | > 0.5, material (≥ 0.60) |
| A7-A1b-rev | Same, reverse direction (storage → driven) | < 0.5 |
| A7-A1b-lag0 | Lag-0 CMI, latch pairs (change coding) | Below the lag-1 value |
| A7-A5 | Delay profile of latch pairs | Highest AUC at delay 1 |
| A7-A6 | Orientation of latch pairs by TE asymmetry | Above the event-rate rule |
| A7-A1 | Lag-1 TE, delayed-only influence vs none (P1) | Reported; the first program showed this class contains no latch pairs, so no strong prediction (expected weakly above 0.5) |

A7 replicates the latch result if A7-A1b, A7-A1b-rev, A7-A5 and A7-A6 all
come out as predicted.
