# handy-genomics
Do a more robust genomics data analysis, sometimes by a surprising simple script!

High-throughput sequencing devices are overtaking analysis capabilities, and sequence data may not get analysed as they should be.

## RNU4-2 scanner
The RNU4-2 is a novel disease causing gene. This non-coding gene is not covered by exome kits. However, there is a small chance of off-target reads that could have captured the gene. Causative variants of RNU4-2 are condensed in a 18-bp region, called stem loop. RNU4-2 scanner uses Samtools to search for any possible RNU4-2 stem loop variants.
