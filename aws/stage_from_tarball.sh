#!/usr/bin/env bash
# stage_from_tarball.sh — pull chosen blocks out of a haploblocks.org per-chromosome tarball.
#
# The chr21/chr22 *_all_seq.tar.gz files (about 40 GB each) are not ordered by position,
# so the whole stream is read once and only the listed blocks are written, then copied
# to DEST (s3:// or a local directory) together with a block list for `submit`.
#
#   bash aws/stage_from_tarball.sh TARBALL_URL BLOCK_NAMES_FILE DEST
#   e.g. https://data.haploblocks.org/haploblock_sequences/1000G/chr21_all_seq.tar.gz \
#        aws/blocks/chr21_pilot_20-21Mb.txt s3://BUCKET/blocks/chr21_pilot
set -euo pipefail
URL="${1:?TARBALL_URL}"; NAMES="${2:?BLOCK_NAMES_FILE}"; DEST="${3:?DEST}"
WORK=$(mktemp -d /tmp/stage.XXXX); trap 'rm -rf "$WORK"' EXIT
grep -v '^\s*$' "$NAMES" | sed 's#.*#*/&_all_seqs.fasta#' > "$WORK/patterns.txt"
echo "Streaming $URL for $(wc -l < "$WORK/patterns.txt") blocks"
# tar exits non-zero when a listed block is absent; the count check below reports that.
curl -fsSL --retry 5 "$URL" | tar -xz -C "$WORK" --wildcards --no-anchored -T "$WORK/patterns.txt" 2> "$WORK/tar.err" || true
: > "$WORK/list.txt"
while read -r f; do
  name=$(basename "$f")
  case "$DEST" in s3://*) aws s3 cp --only-show-errors "$f" "$DEST/$name" ;; *) mkdir -p "$DEST" && cp "$f" "$DEST/$name" ;; esac
  echo "$DEST/$name" >> "$WORK/list.txt"
done < <(find "$WORK" -name '*_all_seqs.fasta' | sort)
FOUND=$(wc -l < "$WORK/list.txt")
case "$DEST" in s3://*) aws s3 cp --only-show-errors "$WORK/list.txt" "$DEST/list.txt" ;; *) cp "$WORK/list.txt" "$DEST/list.txt" ;; esac
echo "Staged $FOUND of $(wc -l < "$WORK/patterns.txt") blocks; block list at $DEST/list.txt"
grep -v "Exiting with failure" "$WORK/tar.err" || true
[ "$FOUND" -gt 0 ] || { echo "FATAL: nothing staged"; exit 1; }
