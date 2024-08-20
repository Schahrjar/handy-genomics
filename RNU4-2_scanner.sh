#!/bin/bash



# The RNU4-2 is a novel disease causing gene in NDD:
# https://www.nature.com/articles/s41588-024-01882-9
# This non-coding gene is not covered by exome kits.
# However, there is a small chance of off-target reads that could have captured the gene.
# Causative variants of RNU4-2 are condensed in a 18-bp region, called stem loop.
# Here, Samtools is applied to screen RNU4-2 stem loop for any possible variants.
#
# By Shahryar Alavi,
# August 2024



# ========== ========== ========== ========== ==========
# Before running the script, define these variables:

# list paths to BAM files in a text file (one BAM file directory in each line). this could be done using, for example:
# ls -1 /path/to/BAMs > list_BAM-dirs.txt
listBAMdirs=list_BAM-dirs.txt

# Path to your hg38 reference genome:
RefGenomeF=/path/to/Homo_sapiens_assembly38.fasta
# If BAM files are mapped to annother reference assembly (e.g. hg19) then the coordinate of stem loop region (below) must be adjusted.

# Set the coordinate of stem loop region:
RegionCoordinate=chr12:120291825-120291842
# This is chr12:120729628-120729645 if BAMs are mapped to hg19.

# Outputs will be saved in RNU4-2_StemLoop_scanner.tsv file.
# ========== ========== ========== ========== ==========



printf "BAM\tchr12position\tVariant\tDepth\n" > RNU4-2_StemLoop_scanner.tsv

for BAMdir in $(cat $listBAMdirs); do
    SampleScreeningResult=$(samtools mpileup -f $RefGenomeF -r $RegionCoordinate $BAMdir | awk -F '\t' '{print $2 "\t" $5 "\t" $4}' | grep "+\|*\|#\|A\|C\|G\|T\|a\|c\|g\|t")
    printf "$BAMdir\t$SampleScreeningResult\n" >> RNU4-2_StemLoop_scanner.tsv
done

exit
