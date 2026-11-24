#!/usr/bin/env bash
# haploblocks_aws.sh — deploy, run and tear down the haploblock PGGB test on AWS Batch.
#
#   bash aws/haploblocks_aws.sh deploy            # CloudFormation stack (Spot Batch, ECR, S3)
#   bash aws/haploblocks_aws.sh push              # build aws/Dockerfile and push to ECR
#   bash aws/haploblocks_aws.sh submit LIST [N_HAPS] [DEDUPE]
#                                                 # one array task per line of LIST
#   bash aws/haploblocks_aws.sh stage CHROM NAMES_FILE
#                                                 # pull named blocks from the CHROM all_seq tarball
#                                                 # into S3 (one job; reads ~40 GB once)
#   bash aws/haploblocks_aws.sh status [JOB_ID]
#   bash aws/haploblocks_aws.sh collect           # all metrics.json -> results/aws_metrics.tsv
#   bash aws/haploblocks_aws.sh teardown          # delete the stack (results bucket is kept)
#
# Needs: aws CLI v2 with credentials, docker. Region from AWS_REGION (default us-west-2,
# where the HPRC bucket lives). STACK (default pangenome-haploblocks), MAX_VCPUS (default 128).
set -euo pipefail
cd "$(dirname "$0")/.."

STACK="${STACK:-pangenome-haploblocks}"
export AWS_REGION="${AWS_REGION:-us-west-2}"
TAG="v0.6.0-hb1"

out() { aws cloudformation describe-stacks --stack-name "$STACK" \
          --query "Stacks[0].Outputs[?OutputKey=='$1'].OutputValue" --output text; }

case "${1:-}" in
deploy)
  VPC=$(aws ec2 describe-vpcs --filters Name=isDefault,Values=true --query 'Vpcs[0].VpcId' --output text)
  [ "$VPC" != "None" ] || { echo "No default VPC: pass SubnetIds/SecurityGroupId to the template yourself."; exit 1; }
  SUBNETS=$(aws ec2 describe-subnets --filters Name=vpc-id,Values="$VPC" Name=map-public-ip-on-launch,Values=true \
              --query 'Subnets[].SubnetId' --output text | tr '\t' ',')
  SG=$(aws ec2 describe-security-groups --filters Name=vpc-id,Values="$VPC" Name=group-name,Values=default \
         --query 'SecurityGroups[0].GroupId' --output text)
  aws cloudformation deploy --stack-name "$STACK" --template-file aws/batch-stack.yaml \
    --capabilities CAPABILITY_IAM \
    --parameter-overrides "SubnetIds=$SUBNETS" "SecurityGroupId=$SG" "MaxvCpus=${MAX_VCPUS:-128}" "ImageTag=$TAG"
  aws cloudformation describe-stacks --stack-name "$STACK" --query 'Stacks[0].Outputs' --output table ;;
push)
  REPO=$(out RepositoryUri)
  aws ecr get-login-password | docker login --username AWS --password-stdin "${REPO%%/*}"
  docker build --platform linux/amd64 -f aws/Dockerfile -t "$REPO:$TAG" .
  docker push "$REPO:$TAG" ;;
submit)
  LIST="${2:?usage: submit LIST [N_HAPS] [DEDUPE]}"; N_HAPS="${3:-0}"; DEDUPE="${4:-1}"
  BUCKET=$(out ResultsBucket)
  N=$(grep -cv '^\s*$' "$LIST")
  KEY="lists/$(basename "$LIST")"
  aws s3 cp "$LIST" "s3://$BUCKET/$KEY"
  ENV="[{\"name\":\"BLOCK_LIST\",\"value\":\"s3://$BUCKET/$KEY\"},{\"name\":\"N_HAPS\",\"value\":\"$N_HAPS\"},{\"name\":\"DEDUPE\",\"value\":\"$DEDUPE\"}]"
  ARRAY=(); [ "$N" -gt 1 ] && ARRAY=(--array-properties "size=$N")
  aws batch submit-job --job-name "hb-$(basename "$LIST" .txt)-n$N_HAPS-dd$DEDUPE" \
    --job-queue "$(out JobQueue)" --job-definition "$(out JobDefinition)" "${ARRAY[@]}" \
    --container-overrides "{\"environment\":$ENV}" --query jobId --output text ;;
stage)
  CHROM="${2:?usage: stage CHROM NAMES_FILE}"; NAMES="${3:?usage: stage CHROM NAMES_FILE}"
  BUCKET=$(out ResultsBucket); SET=$(basename "$NAMES" .txt)
  aws s3 cp "$NAMES" "s3://$BUCKET/lists/names-$SET.txt"
  CMD="aws s3 cp s3://$BUCKET/lists/names-$SET.txt /tmp/names.txt && bash aws/stage_from_tarball.sh https://data.haploblocks.org/haploblock_sequences/1000G/${CHROM}_all_seq.tar.gz /tmp/names.txt s3://$BUCKET/blocks/$SET"
  aws batch submit-job --job-name "hb-stage-$SET" --job-queue "$(out JobQueue)" --job-definition "$(out JobDefinition)" \
    --container-overrides "{\"command\":[\"bash\",\"-c\",\"$CMD\"]}" --query jobId --output text
  echo "When it succeeds: aws s3 cp s3://$BUCKET/blocks/$SET/list.txt aws/blocks/$SET.s3.txt, then submit that list." ;;
status)
  if [ -n "${2:-}" ]; then
    aws batch describe-jobs --jobs "$2" --query 'jobs[0].{status:status,array:arrayProperties.statusSummary,reason:statusReason}'
  else
    for s in SUBMITTED PENDING RUNNABLE STARTING RUNNING SUCCEEDED FAILED; do
      printf '%-10s %s\n' "$s" "$(aws batch list-jobs --job-queue "$(out JobQueue)" --job-status "$s" --query 'length(jobSummaryList)')"
    done
  fi ;;
collect)
  BUCKET=$(out ResultsBucket); mkdir -p results/aws
  aws s3 cp --recursive --only-show-errors "s3://$BUCKET/results/" results/aws/ --exclude '*' --include '*/metrics.json'
  python3 - <<'PY'
import glob, json, csv
rows = [json.load(open(p)) for p in sorted(glob.glob("results/aws/*/*/metrics.json"))]
if rows:
    with open("results/aws_metrics.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
print(f"{len(rows)} runs -> results/aws_metrics.tsv")
PY
  ;;
teardown)
  aws cloudformation delete-stack --stack-name "$STACK"
  aws cloudformation wait stack-delete-complete --stack-name "$STACK"
  echo "Stack deleted. The results bucket was kept; delete it yourself when done." ;;
*) sed -n '2,13p' "$0"; exit 1 ;;
esac
