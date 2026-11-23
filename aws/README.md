# Testing on AWS with haploblocks.org recombination-split blocks

[data.haploblocks.org](https://data.haploblocks.org/) (CC0) splits the 1000 Genomes
haplotypes (2,548 samples, 5,096 haplotypes) at recombination hotspots and publishes one
multi-FASTA per block. Each block is a natural, independent chunk for this project's
parallel build: one AWS Batch array task builds one PGGB graph per block.

## What is in the data

| Path | What it is | Size |
|---|---|---|
| `chr22_small_dataset/` | 12 random chr22 blocks, each `_all_seqs.fasta` (all 5,096 haplotypes), `_rep_seq.fasta`, `_cluster.tsv` | 23 to 279 MB per block |
| `haploblock_sequences/1000G/chr21_all_seq.tar.gz` | every chr21 block, all haplotypes | 41 GB |
| `haploblock_sequences/1000G/chr22_all_seq.tar.gz` | every chr22 block, all haplotypes | 38 GB |
| `haploblock_sequences/1000G/chr*_rep_seq.tar.gz` | cluster representatives only | 3 to 28 GB |
| `bidirectional_blast_samples/block_stats.tsv` | per-block coordinates, cluster counts, entropy | 3 MB |

Facts checked on 2026-11-23 that matter for the pipeline:

- Headers look like `HG00096_chr22_region_44303006-44307622_hap0`. Consecutive blocks share
  their boundary coordinate, so blocks touch but do not overlap.
- Some files carry header-only records (2 in the 4.6 kb chr22 block, 68 in the 11 kb one).
  `pipeline/haploblocks/prepare_block.py` drops them.
- Sequences are consensus sequences built from the phased small-variant panel: lengths vary
  only by small indels (4,617 to 4,620 bp in one block), so most structural variants are
  not represented.
- Many haplotypes are identical within a block (5,096 haplotypes collapse to 318 distinct
  sequences in the 4.6 kb block). `--dedupe` builds the graph once per distinct sequence and
  writes a table mapping every haplotype to its representative.
- Block lengths: chr21 597 blocks (median 41 kb), chr22 669 blocks (median 30 kb), with
  25% and 40% under 25 kb. The repo's default `-s 5000` makes wfmash fail outright on blocks
  shorter than 5 kb ("reference sketch is empty"), so these runs use
  `config/haploblocks.yaml` (`-s 1000`). On the 4.6 kb block, `-s 1000`, `500` and `300`
  gave the same graph (119 nodes, 157 edges, 32 paths).
- The repo's 1 Mb pilot region, chr21:20,000,000-21,000,000, is covered by 20 blocks
  (`blocks/chr21_pilot_20-21Mb.txt`).

## Run it

Needs the AWS CLI v2 with credentials for your account, and Docker.

```bash
export AWS_REGION=us-west-2          # same region as the HPRC bucket
bash aws/haploblocks_aws.sh deploy   # Spot Batch compute env, queue, job def, ECR repo, S3 bucket
bash aws/haploblocks_aws.sh push     # build aws/Dockerfile (pinned PGGB v0.6.0 + awscli) and push

# 1. Smoke test: 12 chr22 blocks, 64 haplotypes each, identical sequences collapsed
bash aws/haploblocks_aws.sh submit aws/blocks/chr22_small_dataset.txt 64 1
bash aws/haploblocks_aws.sh status

# 2. Scaling: same blocks at 256 and 1024 haplotypes, then all 5,096
bash aws/haploblocks_aws.sh submit aws/blocks/chr22_small_dataset.txt 256 1
bash aws/haploblocks_aws.sh submit aws/blocks/chr22_small_dataset.txt 0 1

# 3. The 1 Mb pilot region, from the chr21 tarball (one job streams 41 GB once)
bash aws/haploblocks_aws.sh stage chr21 aws/blocks/chr21_pilot_20-21Mb.txt
# when that job succeeds, fetch the list it prints and submit it like step 1

bash aws/haploblocks_aws.sh collect   # every metrics.json -> results/aws_metrics.tsv
bash aws/haploblocks_aws.sh teardown  # deletes compute, queue, ECR; keeps the results bucket
```

Each task writes `metrics.json` (wall time, peak memory, instance type, nodes, edges,
paths, graph bp), `final.gfa.gz`, the PGGB log and the dedupe table to
`s3://<bucket>/results/n<N>_dd<0|1>/<block>/`.

To try the job container without AWS, set `RESULTS_S3` to a local directory:

```bash
docker build -f aws/Dockerfile -t pangenome-haploblocks .
echo https://data.haploblocks.org/chr22_small_dataset/chr22_27020234-27031341_all_seqs.fasta > /tmp/list.txt
docker run --rm -v /tmp:/tmp -e BLOCK_LIST=/tmp/list.txt -e RESULTS_S3=/tmp/hb-results -e N_HAPS=256 pangenome-haploblocks
```

## Cost controls

- Compute is EC2 Spot through AWS Batch, scaled from 0 vCPUs, capped by `MAX_VCPUS`
  (default 128). AWS Batch itself has no extra charge; you pay for the instances, S3 and ECR.
- Spot interruptions are retried up to 3 times (`Host EC2*` status reason). AWS
  recommends Spot for jobs of roughly 30 minutes or less, which per-block jobs are.
- Every task has a 6-hour timeout. `teardown` removes everything that costs money by the
  hour; the results bucket is kept on purpose.
