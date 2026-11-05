# 13 The locking gradient within one wiring class; the matched-bound inversion matched on history (R7)

Two limits named in the thesis after its independent read, tested with predictions hash-logged first.
- `design/predictions_R7.md` (sha256 be66a622…). R7.1 (gradient within pull-down pairs) not met: same direction,
  interval includes 0, few targets active in both phases have a pull-down wired source. R7.2 (within pass-only
  pairs) met. R7.3 (matched on bound and history, strictly locked pairs stay below 0.5) met.
- `locking_gradient_within_wiring_class_R7.py`: the test; it reads the saved reference and bound arrays from the
  data folder. Run of record `R7_locking_within_class_run_of_record.log`; a rerun added one descriptive block on
  A3's own scope (all same-phase locked pairs), where matching on history raises the AUC from 0.400 to 0.465.
- `RESULTS_R7.md`: the numbers and the reading.
Outputs: `results/13_locking_within_class_R7/`.
