#!/usr/bin/env python3
"""HG002 truth per haplotype from the GIAB HG002 Q100 v5.0q GRCh38 benchmark (small variants + SVs).

GIAB built v5.0q with dipcall from the HG002 Q100 v1.1 assembly; its dipcall hap1 is the paternal
assembly (chr*_PATERNAL), so GT a|b gives HG002#1 allele a and HG002#2 allele b. Only PASS records.

Writes, for HAP in 1 2:
  OUT_HG002_HAP.vcf   haploid calls (smvar + stvar), for truvari on all sizes
  OUT_HG002_HAP.bed   smvar benchmark regions inside the region (all-sizes scoring)
  OUTSV_HG002_HAP.vcf the same calls; OUTSV_HG002_HAP.bed the stvar benchmark regions (SV scoring)

    giab_truth.py SMVAR.vcf.gz STVAR.vcf.gz SMVAR.bed STVAR.bed CHROM START END OUT OUTSV
"""
import sys

import pysam


def records(vcf, chrom, start, end):
    for r in pysam.VariantFile(vcf).fetch(chrom, start, end):
        if list(r.filter.keys()) not in ([], ["PASS"]):
            continue
        gt = r.samples[0]["GT"]
        if len(gt) != 2:
            continue
        for h, a in enumerate(gt, 1):
            if a in (None, 0):
                continue
            alt = r.alleles[a]
            if alt == "*" or set(alt.upper()) - set("ACGTN"):
                continue
            yield h, (r.pos, r.ref.upper(), alt.upper())


def bed(path, chrom, start, end, out):
    with open(out, "w") as fh:
        for line in open(path):
            c, s, e = line.split()[:3]
            s, e = max(int(s), start), min(int(e), end)
            if c == chrom and e > s:
                fh.write(f"{c}\t{s}\t{e}\n")


def main():
    smvar, stvar, smbed, stbed, chrom, start, end, out, outsv = sys.argv[1:]
    start, end = int(start), int(end)
    calls = {1: set(), 2: set()}
    for vcf in (smvar, stvar):
        for h, rec in records(vcf, chrom, start, end):
            if start < rec[0] <= end:
                calls[h].add(rec)
    for h in (1, 2):
        for prefix in (out, outsv):
            with open(f"{prefix}_HG002_{h}.vcf", "w") as fh:
                fh.write("##fileformat=VCFv4.2\n")
                fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n")
                for p, r, a in sorted(calls[h]):
                    fh.write(f"{chrom}\t{p}\t.\t{r}\t{a}\t60\tPASS\t.\tGT\t1\n")
        bed(smbed, chrom, start, end, f"{out}_HG002_{h}.bed")
        bed(stbed, chrom, start, end, f"{outsv}_HG002_{h}.bed")
        print(f"HG002#{h}: {len(calls[h])} calls, {sum(len(r) != len(a) and max(len(r), len(a)) >= 50 for _, r, a in calls[h])} >= 50 bp")


if __name__ == "__main__":
    main()
