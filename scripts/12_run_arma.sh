#!/usr/bin/env bash
# Launch Arm A: one Exomiser batch worker per batch file, under tmux.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate rdacmg
E=$HOME/rdacmg/resources/exomiser/exomiser-cli-15.1.0
R=$HOME/rdacmg/repo
B=$HOME/rdacmg/batches
L=$HOME/rdacmg/logs; mkdir -p $L

n=0
for f in $B/*.batch; do
  bn=$(basename "$f" .batch)
  prop=${bn%_w*}
  tmux kill-session -t "ex_$bn" 2>/dev/null
  tmux new-session -d -s "ex_$bn" \
    "java -Xmx6g -Dspring.config.location=file:$R/config/${prop}.properties \
     -jar $E/exomiser-cli-15.1.0.jar batch $f > $L/${bn}.log 2>&1; echo DONE_$bn >> $L/done.txt"
  n=$((n+1))
done
echo "launched $n Exomiser batch workers"
