# Pinned PGGB v0.6.0 (same digest as config/pipeline.yaml) plus awscli and this repo's scripts.
# The base image already has python3 (with pyyaml), samtools, /usr/bin/time and curl.
# Build from the repository root:  docker build -f aws/Dockerfile -t pangenome-haploblocks .
FROM ghcr.io/pangenome/pggb@sha256:44f3563c5c3e2032c293fe309eed178b754021b07963ca4857e4ae5efda86474
RUN python3 -m pip install --no-cache-dir awscli==1.46.3
WORKDIR /opt/pangenome
COPY config/ config/
COPY scripts/gen_pggb_config.py scripts/
COPY pipeline/haploblocks/ pipeline/haploblocks/
COPY aws/run_block.sh aws/stage_from_tarball.sh aws/
CMD ["bash", "aws/run_block.sh"]
