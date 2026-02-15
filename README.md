[![GitHub release](https://img.shields.io/github/v/release/Schahrjar/handy-genomics)](https://github.com/Schahrjar/handy-genomics/releases/latest)
[![last commit](https://img.shields.io/github/last-commit/Schahrjar/handy-genomics)](https://github.com/Schahrjar/handy-genomics/commits/main)
[![Downloads](https://img.shields.io/github/downloads/Schahrjar/handy-genomics/total?style=flat-square)](https://github.com/Schahrjar/handy-genomics/releases)
[![License](https://img.shields.io/badge/license-MIT-blue)](https://opensource.org/license/mit)

# Handy Genomics
Do a more robust genomics data analysis, sometimes by a surprising simple script!

High-throughput sequencing devices are overtaking analysis capabilities, and sequence data may not get analysed as they should be. Here are some simple scripts for geneticists having some human sequence data for analysis.

## RNU4-2 scanner
The RNU4-2 is a novel disease causing gene. This non-coding gene is not covered by exome kits. However, there is a small chance of off-target reads that could have captured the gene. Causative variants of RNU4-2 are condensed in a 18-bp region, called stem loop. RNU4-2 scanner uses Samtools to search for any possible RNU4-2 stem loop variants.

> [!IMPORTANT]
> **update 2025-06-11:**
>
> It is developed to [RNUscanner](https://github.com/Schahrjar/RNUscanner), a more comprehensive tool with easy to use instructions, and scalability to include other RNU genes (RNU2-2, RNU5B-1, RNU5A-1, etc.).

## RNA-Seq plotting
RNA-Seq is widely used for transcriptomics, from research to clinics. Performing appropriate analyses and creating high quality plots unleashes the great value of RNA-Seq data.\
**update 2025-10-28:** I introduced enhancements in plotting and pathway analysis.

## Haplotype scanner
Phasing variants of WGS data using SHAPEIT gives a BCF file with samples' genotypes. Then next step would be infering which samples carry a specific haplotype, defined by a set of SNPs. The HapScanner gets a set of SNPs and scans the genotypes from SHAPEIT output and quickly returns all samples' genotypes for the defined haplotype block.
