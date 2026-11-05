# Transfer Entropy on a Microprocessor

Code, designs and results for the Honours thesis *Transfer Entropy on a Microprocessor: Pairwise and
Multivariate Inference of Signal Flow in the MOS 6502* (Priyongkar Kundu, School of Computer Science,
University of Sydney, 2026; supervisor Dr Joseph Lizier).

The thesis scores pairwise and multivariate transfer entropy against the netlist of the MOS 6502, using
a transistor-level recording of the chip running Donkey Kong, a switch-level simulator with
interventions, synthetic circuits and networks, two further recorded programs, and clamps applied to the
recorded run itself.

## How the repository is organised

The analyses are numbered in the order they were done; each folder holds its design files (the
predictions fixed before the run, in `design/`), its scripts and a results summary. Outputs are in
`results/`, mirrored by number.

```
src/te6502/                      the shared library
  transfer_entropy_estimators.py   plug-in MI, AIS, TE and CMI, null tests, metrics, IDTxl wrapper
  switch_level_simulator.py        Python port of the Visual6502 switch-level simulation, with clamps
  test_helpers.py                  bootstrap, stratified AUC, decision rule, paths for A1-A19
  paths.py                         where the code, the large data and the outputs live

analyses/
  00_main_run_notebooks/                                  the pre-specified main run (H1-H8) and supplementary notebook
  01_settling_versus_codrive_on_the_simulated_chip/       an early chip analysis whose estimator the tests reuse
  02_latch_transfer_tests_A1_to_A8/                       what lag-1 TE recovers: latch transfer (chip, recording, NMOS)
  03_multivariate_idtxl_runs_A4_R1_A17/                   IDTxl with the phase supplied, the H6 re-draw, the aligned rerun
  04_followership_mechanism_A9_to_A11/                    the exact within-target account and its synthetic generalisation
  05_consequences_of_the_account_A12_to_A19/              dose-response, history coverage, remainder, timed clamps, block of two
  06_regenerating_recordings_from_the_original_simulator/ Sim2600 reproduces the published recording; Pitfall, Space Invaders
  07_replication_on_pitfall_and_space_invaders_R2/        out-of-sample replication on two further programs
  08_timed_clamps_on_the_recorded_program_A20/            the recording's latch pairs confirmed by intervention
  09_recomputing_the_supplementary_results/              lost supplementary outputs recomputed and compared with the thesis
  10_phase_aware_remedy_on_held_out_programs_R3/         the phase-aware remedy tested on the two held-out programs
  11_phase_relation_against_physical_delay_R4/          lag 1 follows the phase relation, not the physical delay
  12_zero_lag_partial_correlation_comparator_R6/        the benchmarks' best contemporaneous method, on the recording
  13_locking_gradient_within_wiring_class_R7/           the locking gradient within one wiring class; bias matched on history

results/                         outputs (JSON and logs) of each analysis, numbered as above; figures/
simulator/                       the patch that ports Sim2600 to Python 3 and adds forced-wire clamps
provenance/code_exactly_as_run/  the scripts byte-for-byte as they were run (see below)
docs/                            file renames, data layout, how to run
```

## Status of each result

Every result in the thesis carries one status: **pre-specified** (H1-H8, fixed before the main run),
**exploratory**, **supplementary**, **hash-logged** (design written in `design/` before the run; the run's
log records the file's SHA-256), **held out** (run on data no prediction was formed from), or **partly
anticipated**. The design files were not lodged with a third-party registry, so hash-logged tests are
weaker than pre-registration.

## Provenance

The scripts in `analyses/` were reorganised and renamed for this repository; their imports and paths
changed, so their SHA-256 differs from the hashes recorded in the run logs. The files exactly as run, with
their original names, are in `provenance/code_exactly_as_run/`, and they reproduce the logged hashes.
`docs/file_renames.md` maps every original name to its new one.

## Data

Large inputs are not in the repository: the Donkey Kong recording (Jonas & Kording 2017), the regenerated
recordings, the Visual6502 netlist, the main run's score matrices and the patched Sim2600. Set
`HONOURS_DATA_ROOT` to the folder holding `data/`, `results/v5_.../` and `vendor/Sim2600/`; see
`docs/data_layout.md`. The recordings can be regenerated from Sim2600 with
`analyses/06_.../regenerate_game_recording_with_sim2600.py`.

## Environment

Python 3.10 (`environment/requirements.txt`); JIDT and IDTxl 1.6.1 for the multivariate runs. Every
analysis script with a `--selftest` option runs a check that computes no statistic of the test.
