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
RNA-Seq is widely used for transcriptomics, from research to clinics. Performing appropriate analyses and creating high quality plots unleashes the great value of RNA-Seq data.
