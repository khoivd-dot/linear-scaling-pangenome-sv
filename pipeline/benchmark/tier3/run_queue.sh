#!/bin/bash
# Tier 3 queue, as run on 2026-11-24/25: E1/E2 rep1, E4, E3 scaling rep2, then E1/E2 rep2 — one at a time
# so timings do not overlap. Every step is resumable; rerun this script after a container restart.
set -uo pipefail
T=/home/user/linear-scaling-pangenome-sv/pipeline/benchmark/tier3
cd /home/user/work/t3 && bash $T/run_e1e2.sh >> e1e2.log 2>&1
cd /home/user/work/t3 && bash $T/run_e4.sh >> e4.log 2>&1
cd /home/user/work/scale && bash $T/scaling_rep2.sh project >> scale_build.log 2>&1 && bash $T/scaling_rep2.sh build >> scale_build.log 2>&1
cd /home/user/work/t3 && REPS=2 bash $T/run_e1e2.sh >> e1e2.log 2>&1
echo "queue done $(date -u +%FT%TZ)"
