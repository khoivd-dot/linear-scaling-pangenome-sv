#!/bin/bash
# After run_queue.sh: E5 (recommended human PGGB settings). Rerun after a restart.
set -uo pipefail
cd /home/user/work/t3
until grep -q "queue done" queue.log 2>/dev/null; do sleep 60; done
bash /home/user/linear-scaling-pangenome-sv/pipeline/benchmark/tier3/run_e5.sh >> e5.log 2>&1
