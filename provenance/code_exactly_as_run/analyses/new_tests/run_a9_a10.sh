#!/bin/zsh
# Runs A9 and A10 in parallel (design: PREDICTIONS_A9.md). Logs in output/.
cd "${0:A:h}"
PY=~/opt/anaconda3/envs/honours6502/bin/python
echo "start $(date)" > output/run_a9_a10.log
shasum -a 256 PREDICTIONS_A9.md a9_decompose.py a10_nmos.py common.py nmos.py >> output/run_a9_a10.log
caffeinate -i $PY -u a9_decompose.py > output/A9.log 2>&1 &
caffeinate -i $PY -u a10_nmos.py     > output/A10.log 2>&1 &
wait
echo "finished $(date)" >> output/run_a9_a10.log
for f in A9 A10; do
  echo "== $f: $(grep -c -E 'Traceback|Error' output/$f.log) error lines; last: $(tail -1 output/$f.log)" >> output/run_a9_a10.log
done
