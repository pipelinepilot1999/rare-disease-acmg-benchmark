#!/usr/bin/env bash
set -uo pipefail
bash "$(dirname "${BASH_SOURCE[0]}")/render_config.sh" >/dev/null
E=$HOME/rdacmg/resources/exomiser/exomiser-cli-15.1.0
R=$HOME/rdacmg/repo
L=$HOME/rdacmg/logs; mkdir -p $L; rm -f $L/done_pen.txt
n=0
for f in $HOME/rdacmg/batches_pen/*.batch; do
  bn=$(basename "$f" .batch)
  tmux kill-session -t "pen_$bn" 2>/dev/null
  tmux new-session -d -s "pen_$bn" \
    "java -Xmx6g -Dspring.config.location=file:$R/config/informed.properties \
     -jar $E/exomiser-cli-15.1.0.jar batch $f > $L/pen_${bn}.log 2>&1; echo DONE >> $L/done_pen.txt"
  n=$((n+1))
done
echo "launched $n penetrant-contrast workers"
