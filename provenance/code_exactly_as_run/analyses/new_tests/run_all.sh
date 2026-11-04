#!/bin/zsh
# Runs every new test in order (design: PREDICTIONS.md). Logs in output/.
cd "${0:A:h}"
PY=~/opt/anaconda3/envs/honours6502/bin/python
export JAVA_HOME=~/opt/anaconda3/envs/honours6502/lib/jvm
echo "start $(date)" > output/run_all.log
# Fast analyses in parallel (about 30 minutes)
$PY -u chip.py      > output/chip.log 2>&1 &
$PY -u recording.py > output/recording.log 2>&1 &
$PY -u nmos.py      > output/nmos.log 2>&1 &
wait
echo "fast analyses finished $(date)" >> output/run_all.log
# IDTxl jobs, one at a time with 3 workers (about 3-4 hours)
$PY -u idtxl_runs.py --workers 3 > output/idtxl.log 2>&1
echo "finished $(date)" >> output/run_all.log
for f in chip recording nmos idtxl; do
  echo "== $f: $(grep -c -E 'Traceback|Error' output/$f.log) error lines; last: $(grep -v -i warning output/$f.log | tail -1)" >> output/run_all.log
done
