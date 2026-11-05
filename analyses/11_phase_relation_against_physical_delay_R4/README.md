# 11 Phase relation against physical delay (R4)

Lag-1 TE finds the latch pairs and misses the wiring, but latch pairs differ from the wiring both in
physical delay and in phase relation. For locked units the phase relation fixes the delay in samples, so
the two can only be separated by measuring the physical delay by intervention. R4 crosses the physical
delay from the saved timed clamps (analysis 05, A16) with the phase relation read from the activity, on
both simulated programs, with predictions written and hash-logged first.

- `design/predictions_R4.md`: the four predictions and the decision rule (sha256 39a3c0cc…).
- `phase_relation_against_physical_delay_R4.py`: the test. It reads the A16 timed clamps and the
  program-2 run from the data folder (`HONOURS_DATA_ROOT/analyses/new_tests/output/`), which are too large
  for the repository.
- `RESULTS_R4.md`: all four predictions met on both programs. Lag-1 TE follows the delay in samples, which
  for locked units the phase relation sets, not the physical delay.

Outputs: `results/11_delay_vs_phase_R4/`. The log of a first start that stopped before any statistic (an
assembler not installed) is kept as `R4_delay_vs_phase_failed_start.log`.
