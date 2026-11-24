#!/usr/bin/env bash
# run_block.sh — build one PGGB graph for one haploblock (AWS Batch array task).
#
# Runs INSIDE the image built from aws/Dockerfile (pinned PGGB v0.6.0 + awscli).
# Each array index picks one line of BLOCK_LIST (an https:// or s3:// FASTA).
#
# Environment:
#   BLOCK_LIST        s3:// or https:// URL of a text file, one block FASTA URL per line (required)
#   RESULTS_S3        s3://bucket/prefix for outputs, or a local directory for a test (required)
#   N_HAPS            keep a seeded random subset of N haplotypes (default 0 = all)
#   DEDUPE            1 = collapse identical haplotype sequences first (default 1)
#   SEED              subset seed (default 42)
#   PIPELINE_CONFIG   PGGB parameter file (default config/haploblocks.yaml)
#   THREADS           default: nproc
#   AWS_BATCH_JOB_ARRAY_INDEX  set by Batch; default 0 for a single local run
set -euo pipefail

: "${BLOCK_LIST:?set BLOCK_LIST}" "${RESULTS_S3:?set RESULTS_S3}"
N_HAPS="${N_HAPS:-0}"; DEDUPE="${DEDUPE:-1}"; SEED="${SEED:-42}"
THREADS="${THREADS:-$(nproc)}"; IDX="${AWS_BATCH_JOB_ARRAY_INDEX:-0}"
export PIPELINE_CONFIG="${PIPELINE_CONFIG:-config/haploblocks.yaml}"
cd /opt/pangenome

fetch() { case "$1" in s3://*) aws s3 cp --only-show-errors "$1" "$2" ;; /*) cp "$1" "$2" ;; *) curl -fsSL --retry 5 -o "$2" "$1" ;; esac; }
# RESULTS_S3 may also be a local directory, for a test run without AWS.
put() { case "$2" in s3://*) aws s3 cp --only-show-errors "$1" "$2" ;; *) mkdir -p "$(dirname "$2")" && cp "$1" "$2" ;; esac; }

WORK=$(mktemp -d /tmp/block.XXXX); trap 'rm -rf "$WORK"' EXIT
fetch "$BLOCK_LIST" "$WORK/list.txt"
SRC=$(grep -v '^\s*$' "$WORK/list.txt" | sed -n "$((IDX + 1))p")
[ -n "$SRC" ] || { echo "FATAL: no line $((IDX + 1)) in $BLOCK_LIST"; exit 1; }
BLOCK=$(basename "$SRC" | sed -E 's/_(all_seqs|rep_seq)\.fasta(\.gz)?$//')
RUN="n${N_HAPS}_dd${DEDUPE}"
echo "=== block $BLOCK (index $IDX) run $RUN threads $THREADS ==="

T0=$(date +%s)
fetch "$SRC" "$WORK/raw.fasta"
case "$SRC" in *.gz) gunzip -c "$WORK/raw.fasta" > "$WORK/raw.tmp" && mv "$WORK/raw.tmp" "$WORK/raw.fasta" ;; esac
T_FETCH=$(( $(date +%s) - T0 ))

ARGS=(--n-haps "$N_HAPS" --seed "$SEED"); [ "$DEDUPE" = 1 ] && ARGS+=(--dedupe)
python3 pipeline/haploblocks/prepare_block.py "$WORK/raw.fasta" "$WORK/in.fa" "${ARGS[@]}"
rm -f "$WORK/raw.fasta"
samtools faidx "$WORK/in.fa"
N_SEQS=$(wc -l < "$WORK/in.fa.fai")

# Parameters come from the same config loader as every other runner.
read -r S P K KK J E < <(python3 -c '
import sys; sys.path.insert(0, "scripts")
from gen_pggb_config import load_config
c = load_config()
print(c["segment_length"], c["minimum_identity"], c["mash_kmer"], c["match_length"], c["path_jump_max"], c["edge_jump_max"])')

T1=$(date +%s)
/usr/bin/time -v -o "$WORK/time.txt" \
  pggb -i "$WORK/in.fa" -o "$WORK/out" -t "$THREADS" -n "$N_SEQS" \
       -p "$P" -s "$S" -K "$K" -k "$KK" -j "$J" -e "$E" > "$WORK/pggb.log" 2>&1 \
  || { tail -40 "$WORK/pggb.log"; put "$WORK/pggb.log" "$RESULTS_S3/$RUN/$BLOCK/pggb.FAILED.log"; exit 1; }
T_PGGB=$(( $(date +%s) - T1 ))

GFA=$(ls "$WORK"/out/*final.gfa | head -1)
read -r NODES EDGES PATHS GBP < <(awk -F'\t' '$1=="S"{n++;bp+=length($3)} $1=="L"{e++} $1=="P"||$1=="W"{p++} END{print n+0,e+0,p+0,bp+0}' "$GFA")
RSS_KB=$(awk -F': ' '/Maximum resident/{print $2}' "$WORK/time.txt")
IMDS=http://169.254.169.254/latest
TOKEN=$(curl -fs -m 2 -X PUT "$IMDS/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 60" || true)
ITYPE=$(curl -fs -m 2 -H "X-aws-ec2-metadata-token: $TOKEN" "$IMDS/meta-data/instance-type" || echo unknown)

cat > "$WORK/metrics.json" <<EOF
{"block": "$BLOCK", "source": "$SRC", "n_haps_requested": $N_HAPS, "dedupe": $DEDUPE, "seed": $SEED,
 "n_sequences": $N_SEQS, "threads": $THREADS, "instance_type": "$ITYPE",
 "segment_length": $S, "fetch_s": $T_FETCH, "pggb_wall_s": $T_PGGB, "pggb_max_rss_kb": ${RSS_KB:-0},
 "nodes": $NODES, "edges": $EDGES, "paths": $PATHS, "graph_bp": $GBP}
EOF
cat "$WORK/metrics.json"

gzip -c "$GFA" > "$WORK/final.gfa.gz"
DEST="$RESULTS_S3/$RUN/$BLOCK"
for f in metrics.json final.gfa.gz pggb.log time.txt; do put "$WORK/$f" "$DEST/$f"; done
if [ -f "$WORK/in.fa.dedupe.tsv" ]; then put "$WORK/in.fa.dedupe.tsv" "$DEST/dedupe.tsv"; fi
echo "=== done: $DEST ==="
