#!/bin/zsh
# Runs A12-A19 (design: PREDICTIONS_A12.md). IDTxl (A17) in its own chain; the rest sequentially beside it.
cd "${0:A:h}"
PY=~/opt/anaconda3/envs/honours6502/bin/python
export JAVA_HOME=~/opt/anaconda3/envs/honours6502/lib/jvm
echo "start $(date)" > output/run_a12_a19.log
shasum -a 256 PREDICTIONS_A12.md a12_recording.py a13_a16_chip.py a15_synthetic.py a17_idtxl_aligned.py a9_decompose.py common.py idtxl_runs.py >> output/run_a12_a19.log
caffeinate -i $PY -u a17_idtxl_aligned.py --workers 3 > output/A17.log 2>&1 &
( caffeinate -i $PY -u a12_recording.py > output/A12.log 2>&1
  echo "A12 recording finished $(date)" >> output/run_a12_a19.log
  caffeinate -i $PY -u a15_synthetic.py > output/A15.log 2>&1
  echo "A15 finished $(date)" >> output/run_a12_a19.log
  caffeinate -i $PY -u a13_a16_chip.py > output/A13_A16.log 2>&1
  echo "A13-A16 chip finished $(date)" >> output/run_a12_a19.log ) &
wait
echo "all finished $(date)" >> output/run_a12_a19.log
for f in A12 A15 A13_A16 A17; do
  echo "== $f: $(grep -c -E 'Traceback|Error' output/$f.log) error lines; last: $(grep -v -i warning output/$f.log | tail -1)" >> output/run_a12_a19.log
done
