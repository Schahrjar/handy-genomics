#!/usr/bin/env python3



# Shahryar Alavi

# Department of Neurodegenerative Diseases,
# UCL Queen Square Institute of Neurology,
# University College London,
# London, UK

# 2025-10-27
# s.alavi@ucl.ac.uk



import argparse
import gzip

def parse_args():
    parser = argparse.ArgumentParser(description="Identify samples carrying a specific phased haplotype by defining a set of SNPs")
    parser.add_argument("--vcf", required=True, help="Phased VCF file (.vcf or .vcf.gz); the VCF output of SHAPEIT")
    parser.add_argument("--snps", required=True, help="Comma-separated dbSNP IDs (rsIDs)")
    parser.add_argument("--pattern", required=True, help="Comma-separated pattern of 0 or 1 indicating REF or ALT corresponding to the defined SNPs")
    return parser.parse_args()

def open_vcf(vcf_path):
    return gzip.open(vcf_path, "rt") if vcf_path.endswith(".gz") else open(vcf_path)

def main():
    args = parse_args()
    snps = args.snps.split(",")
    pattern = [int(x) for x in args.pattern.split(",")]

    # Store genotypes for the selected SNPs
    genotypes = {}
    samples = []
    with open_vcf(args.vcf) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                for s in snps:
                    genotypes[s] = {}
                continue

            fields = line.strip().split("\t")
            snp_id = fields[2]
            if snp_id not in snps:
                continue

            gt_values = [x.split(":")[0] for x in fields[9:]]
            genotypes[snp_id] = dict(zip(samples, gt_values))

    print("sample\tstatus")

    for s in samples:
        hap1, hap2 = [], []
        for snp in snps:
            gt = genotypes[snp].get(s, "./.")
            a1, a2 = gt.split("|") if "|" in gt else (".", ".")
            hap1.append(0 if a1 == "0" else 1)
            hap2.append(0 if a2 == "0" else 1)

        match1 = hap1 == pattern
        match2 = hap2 == pattern

        if match1 and match2:
            status = "homozygous"
        elif match1 or match2:
            status = "heterozygous"
        else:
            status = "none"

        print(f"{s}\t{status}")

if __name__ == "__main__":
    main()
