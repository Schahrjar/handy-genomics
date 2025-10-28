#!/usr/bin/env Rscript



# Shahryar Alavi

# Department of Neurodegenerative Diseases
# UCL Queen Square Institute of Neurology
# University College London
# London, UK

# Updated 2025-07-15
# s.alavi@ucl.ac.uk



# Rare genetic diseases are mostly caused by causative variants. Among them are splicing variants or variants of genes that play major roles in cellular pathways. We have used RNA-Seq method to study the consequences of such genomic variants at transcriptomics level.
# To analyse RNA-Seq data, we have two approaches: first, we QC the data to make sure if it can reliably reflect our traits of study. Second, we not only seek for differentially expressed genes, but also look for cellular pathways affected by the deficient gene. The later could potentially suggest treatments and patient care.




# Define variables:
work.dir <- "Path/to the folder/where each gene count file/ is stored in a folder / with folder name as sample name/"

# Load (and instaall) required libraries:
library(data.table)
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
library(DESeq2)
library(dplyr) # For data manipulation like arrange, head, etc.
library(ggplot2) # For plotting
library(ggrepel) # For non-overlapping labels
library(patchwork) # For combining plots
library(org.Hs.eg.db) # For gene annotation
library(AnnotationDbi) # For gene annotation
library(clusterProfiler) # For gene enrichment
library(enrichplot) # For GSEA plots
library(ReactomePA) # For Reactome pathways


# Load genes counts data (analysed using STAR):
# Initialize a list to store data.tables
gene_data_list <- list()

for (GeneRead.path in list.files(path = work.dir, pattern = "ReadsPerGene.out.tab", full.names = TRUE, recursive = TRUE)) {
  sample.name <- basename(dirname(GeneRead.path)) # This extracts the directory name as the sample ID
  
  # Load gene reads count data (V1=Gene, V4=stranded-reverse counts)
  gene.reads <- fread(GeneRead.path, sep = "\t", header = FALSE, skip = 4, select = c(1, 4))
  
  
  # Rename columns to 'Gene' and the actual sample.name.
  # data.table::setnames preserves the exact sample.name even if it starts with a number.
  setnames(gene.reads, c("V1", "V4"), c("Gene", sample.name))
  
  # Store the data.table in the list
  gene_data_list[[sample.name]] <- gene.reads
}

# Combine all samples' read counts more efficiently and keep original names
# Use Reduce with merge to combine all data.tables in the list
# This ensures that column names are kept as they were set
if (length(gene_data_list) > 0) {
  combined_dt <- Reduce(function(dt1, dt2) merge(dt1, dt2, by = "Gene", all = TRUE), gene_data_list)
} else {
  stop("No gene count files found to process.")
}

# Fill NA values with 0 (for genes not detected in all samples)
combined_dt[is.na(combined_dt)] <- 0

# Convert the combined data.table to a matrix.
# The `rownames = "Gene"` argument sets the 'Gene' column as row names.
df.expression.samples <- as.matrix(combined_dt, rownames = "Gene")

# Ensure counts are integers
df.expression.samples <- round(df.expression.samples)

# At this point, colnames(df.expression.samples) should contain the original sample names

# Load required libraries
library(rtracklayer)     # For importing GTF
library(GenomicRanges)   # For working with GRanges

# Load samples info
df.samples <- fread("list-samples.tsv", sep = "\t", header = TRUE) # User provided; should include a column as exact sample names and the main category (e.g. patient/control) and covariate (if any).
df.samples <- as.data.frame(df.samples)

# === Critical Sample Order Check ===
# Ensure sample order in df.samples matches column names of df.expression.samples
if (!all(colnames(df.expression.samples) %in% df.samples$Sample)) {
  stop("Not all sample IDs from expression data found in sample information table. Check 'Sample' column in list-samples.tsv and ensure it matches file names.")
}
# Reorder df.samples based on the column names of df.expression.samples
df.samples <- df.samples[match(colnames(df.expression.samples), df.samples$Sample), ]
stopifnot(all(colnames(df.expression.samples) == df.samples$Sample))
# ===   ===   ===

# Build DESeq dataset
dataset.DESeq <- DESeqDataSetFromMatrix(countData = df.expression.samples, colData = df.samples, design = ~ Phenotype + Gender)

# Run DESeq analysis
DiffExp.DESeq <- DESeq(dataset.DESeq)


# Apply LFC shrinkage for more stable and accurate log2 fold changes
# You might need to adjust 'coef' based on your DESeq2 contrast.
# Example: If your categories are 'CASE' and 'CONTROL', and 'CONTROL' is the reference level.
# resultsNames(DiffExp.DESeq) will show something like "Category_CASE_vs_CONTROL"
resultsNames(DiffExp.DESeq)
result.DiffExp <- lfcShrink(DiffExp.DESeq, coef="Phenotype_patient_vs_control", type="ashr")


# Run a quick PCA to check the case-control distances:
# Variance-stabilizing transformation
vsd <- vst(dataset.DESeq, blind = FALSE)

# Assign sample names for plotting: This should directly use the clean column names from vsd
colData(vsd)$Sample <- colnames(vsd)

# PCA data
# Plot PCA also by Carrier (or other relevant metadata) as shape
pcaData <- plotPCA(vsd, intgroup = c("Phenotype"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$Sample <- colData(vsd)[rownames(pcaData), "Sample"]

# PCA plot
plot.pca <- ggplot(pcaData, aes(PC1, PC2, color = Phenotype, shape = Gender, label = Sample)) + # Label uses clean Sample
  geom_point(size = 3) +
  geom_text_repel(
    size = 3.5,
    box.padding = 0.5,
    max.overlaps = Inf,
    segment.color = "grey60",
    color = "black"  # sample label text color
  ) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  ggtitle("Case/Control PCA") +
  coord_fixed() +
  scale_color_brewer(palette = "Dark2") +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = c(0.8, 0.8), # legend position relative to the plot area
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_blank()
  )

print(plot.pca)
# Save high-res image as PDF (vector graphic)
ggsave("plots/PCA_plot_highres.png", plot = plot.pca, device = "png",
       width = 6, height = 5, units = "in")




# If cases and controls don’t separate on a PCA, that is not automatically bad:
# 1- it can mean the biological effect is small compared with technical/biological noise
# 2- sample size is too small
# 3- we simply don’t have the right signal in the data.





# Create a log2FC table
# 1. Create a data frame from DESeq2 results
log2fc <- data.frame(
  GeneID = rownames(result.DiffExp),
  Log2FC = result.DiffExp$log2FoldChange,
  Pvalue = result.DiffExp$pvalue,
  Padj = result.DiffExp$padj,
  stringsAsFactors = FALSE
)

# 2. Remove rows with NA values
log2fc.clean <- na.omit(log2fc)

# 3. Filter based on adjusted p-value threshold (e.g., FDR < 0.1)
log2fc.filtered <- log2fc.clean[log2fc.clean$Padj < 0.1, ]

# 4. Map Ensembl IDs to HGNC symbols using org.Hs.eg.db
GeneSymbol <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(log2fc.filtered$GeneID),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

# --- Handle one-to-many mappings by keeping only one symbol per Ensembl ID ---
# If an ENSEMBL ID maps to multiple SYMBOLS, keep the first SYMBOL encountered.
# This ensures that each ENSEMBL ID has only one SYMBOL associated with it for the merge.
GeneSymbol <- GeneSymbol[!duplicated(GeneSymbol$ENSEMBL), ]

# Check for mapping loss (this check will now be accurate as GeneSymbol should have unique ENS IDs)
cat(sprintf("Warning: %d of %d significant genes (%.2f%%) were lost during Ensembl to Symbol mapping for table export.\n",
            length(unique(log2fc.filtered$GeneID)) - nrow(GeneSymbol), length(unique(log2fc.filtered$GeneID)),
            (length(unique(log2fc.filtered$GeneID)) - nrow(GeneSymbol)) / length(unique(log2fc.filtered$GeneID)) * 100))

# 5. Merge gene symbols with the filtered DE table
merged.table <- merge(
  log2fc.filtered,
  GeneSymbol,
  by.x = "GeneID",
  by.y = "ENSEMBL",
  all.x = TRUE
)

# 6. Reorder columns: HGNC symbol (or blank), Ensembl ID, Log2FC, P-value, Adjusted P-value
merged.table <- merged.table[, c("SYMBOL", "GeneID", "Log2FC", "Pvalue", "Padj")]

# 7. Sort by adjusted p-value
merged.table <- merged.table[order(merged.table$Padj), ]

# 8. Write table
write.table(
  merged.table,
  file = "Log2FC_with_HGNC.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


# Set thresholds
fc_threshold <- 1.5
pval_threshold <- 0.01

# Copy data
volcano_data <- data.frame(
  GeneID = rownames(result.DiffExp),
  Log2FC = result.DiffExp$log2FoldChange,
  Pvalue = result.DiffExp$pvalue,
  Padj = result.DiffExp$padj,
  stringsAsFactors = FALSE
)
volcano_data <- na.omit(volcano_data)

# Map Ensembl IDs to HGNC symbols for volcano plot labeling
volcano_data_symbols <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(volcano_data$GeneID),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

# --- RESOLUTION: Handle one-to-many mappings by keeping only one symbol per Ensembl ID ---
# If an ENSEMBL ID maps to multiple SYMBOLS, keep the first SYMBOL encountered.
volcano_data_symbols <- volcano_data_symbols[!duplicated(volcano_data_symbols$ENSEMBL), ]

volcano_data <- merge(volcano_data, volcano_data_symbols, by.x = "GeneID", by.y = "ENSEMBL", all.x = TRUE)


# Assign significance labels
volcano_data$Significance <- "Not significant"
volcano_data$Significance[volcano_data$Log2FC > log2(fc_threshold) & volcano_data$Padj < pval_threshold] <- "Upregulated"
volcano_data$Significance[volcano_data$Log2FC < -log2(fc_threshold) & volcano_data$Padj < pval_threshold] <- "Downregulated"


# Plot
plot.volcano <- ggplot(volcano_data, aes(x = Log2FC, y = -log10(Padj), color = Significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "darkgreen", "Downregulated" = "red", "Not significant" = "grey")) +
  geom_vline(xintercept = c(-log2(fc_threshold), log2(fc_threshold)), linetype = "dashed", color = "black") + # Use log2(fc_threshold)
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_line(color = "grey90"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = c(0.05, 0.9),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  labs(
    title = "Differential Gene Expression",
    x = "Log2 Fold Change",
    y = "-log10(Adjusted p-value)"
  ) +
  guides(color = guide_legend(title = NULL)) +
  # Automate ggrepel on top-N significant genes
  geom_text_repel(
    data = subset(volcano_data, Padj < pval_threshold & abs(Log2FC) > log2(fc_threshold)) %>%
      arrange(Padj) %>%
      head(15), # Select top 15 significant genes by Padj
    aes(label = SYMBOL),
    size = 4,
    color = "black",
    box.padding = 0.5,
    max.overlaps = Inf
  )

print(plot.volcano)

# Save as high-res image as PNG
ggsave("Volcano_plot_highres.png", plot.volcano, device = "png", width = 7, height = 5, units = "in")


# Combined PCA and volcano plots:
# Optional: ensure both plots use consistent theme tweaks
plot.pca <- plot.pca +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = c(0.8, 0.8),  # inside plot box
    legend.background = element_rect(fill = "white", color = "black")
  )

plot.volcano <- plot.volcano +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = c(0.5, 0.9),  # inside plot box
    legend.background = element_rect(fill = "white", color = "black")
  )

# Combine plots side by side with guides = "collect"
combined_plot <- plot.pca + plot_spacer() + plot.volcano +
  plot_layout(ncol = 3, widths = c(0.7, 0.07, 1.1), guides = "collect") + 
  plot_annotation(title = "SHTN1 LoF Blood Transcriptome",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)))

print(combined_plot)

# Save high-resolution output as PNG
ggsave("Combined_PCA_Volcano.png", combined_plot, width = 12, height = 5.5, device = "png")


# Trying to setup clusterProfiler for gene enrichment analysis.
# Filter significant genes
sig.genes <- result.DiffExp[which(result.DiffExp$padj < 0.05 & !is.na(result.DiffExp$padj)), ]

# Extract Ensembl IDs (already cleaned during data ingestion)
ensembl.ids <- rownames(sig.genes)

# Confirm some look right
cat("First few cleaned Ensembl IDs:\n")
print(head(ensembl.ids))

# Check which are valid keys in org.Hs.eg.db
valid.ensembl.ids <- intersect(ensembl.ids, keys(org.Hs.eg.db, keytype = "ENSEMBL"))

# Map to Entrez IDs
gene.df <- bitr(valid.ensembl.ids,
                fromType = "ENSEMBL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

# Check for mapping loss
cat(sprintf("Warning: %d of %d significant genes (%.2f%%) were lost during Ensembl to Entrez mapping for ORA.\n",
            length(valid.ensembl.ids) - nrow(gene.df), length(valid.ensembl.ids),
            (length(valid.ensembl.ids) - nrow(gene.df)) / length(valid.ensembl.ids) * 100))

# Extract Entrez IDs
entrez.sig.genes <- unique(gene.df$ENTREZID)
# Summary
cat(length(entrez.sig.genes), "Entrez IDs found for ORA.\n")


# Define universe for ORA (all genes that passed DESeq2 filtering and could be mapped to Entrez)
all_genes_in_analysis_ensembl <- rownames(counts(dataset.DESeq)) # Genes after low count filtering
universe_map <- bitr(all_genes_in_analysis_ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe_entrez <- unique(universe_map$ENTREZID)
cat(sprintf("Universe for ORA contains %d Entrez IDs.\n", length(universe_entrez)))


#Run Functional Enrichment
#Use either GO (Gene Ontology) or KEGG/Reactome pathway analysis.
# (A) GO Biological Process (BP) Enrichment

ego <- enrichGO(gene = entrez.sig.genes,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",       # "BP", "MF", or "CC"
                pAdjustMethod = "BH",
                pvalueCutoff = 0.1,
                qvalueCutoff = 0.1, # Can be removed if pvalueCutoff applies to adjusted P
                universe = universe_entrez) # Explicitly setting universe

# KEGG Pathway Enrichment
ekegg <- enrichKEGG(gene = entrez.sig.genes,
                    organism = "hsa",     # Human KEGG code
                    pvalueCutoff = 0.1,
                    universe = universe_entrez) # Explicitly setting universe

# Reactome Pathways
ereactome <- enrichPathway(gene = entrez.sig.genes,
                           organism = "human",
                           pvalueCutoff = 0.1,
                           readable = TRUE,
                           universe = universe_entrez) # Explicitly setting universe


# Plot the Results
# Dotplot
plot_ego_bp <- dotplot(ego, showCategory = 15, font.size = 12, title = "GO: Biological Processes")
ggsave("GO_BP_Dotplot.png", plot = plot_ego_bp, device = "png", width = 8, height = 7)

# Barplot
plot_ekegg <- barplot(ekegg, showCategory = 10, title = "KEGG Pathways")
ggsave("KEGG_Barplot.png", plot = plot_ekegg, device = "png", width = 8, height = 6)

# Reactome dotplot
plot_ereactome <- dotplot(ereactome, showCategory = 10, title = "Reactome Pathways")
ggsave("Reactome_Dotplot.png", plot = plot_ereactome, device = "png", width = 8, height = 6)


# Run ridgeplot using clusterProfiler:
# Start from your DESeq2 results
# 1- Prepare a Ranked Gene List
# Start from your result object
gene_list <- result.DiffExp$log2FoldChange
names(gene_list) <- rownames(result.DiffExp) # Already cleaned Ensembl IDs

# Remove NAs
gene_list <- gene_list[!is.na(gene_list)]

# Check if keys are valid
valid_keys <- intersect(names(gene_list), keys(org.Hs.eg.db, keytype = "ENSEMBL"))

# Keep only valid Ensembl IDs
gene_list <- gene_list[names(gene_list) %in% valid_keys]

# Now do the mapping
gene_map <- bitr(
  names(gene_list),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
# Check for mapping loss
cat(sprintf("Warning: %d of %d genes (%.2f%%) were lost during Ensembl to Entrez mapping for GSEA.\n",
            length(names(gene_list)) - nrow(gene_map), length(names(gene_list)),
            (length(names(gene_list)) - nrow(gene_map)) / length(names(gene_list)) * 100))


# Merge fold changes with Entrez mapping
gene_df <- merge(
  data.frame(ENSEMBL = names(gene_list), log2FC = gene_list),
  gene_map,
  by = "ENSEMBL"
)

# Remove duplicate Entrez IDs (keep max abs(log2FC))
gene_df <- gene_df[order(abs(gene_df$log2FC), decreasing = TRUE), ]
gene_df <- gene_df[!duplicated(gene_df$ENTREZID), ]

# Rebuild named vector for GSEA
gene_list_unique <- gene_df$log2FC
names(gene_list_unique) <- gene_df$ENTREZID
gene_list_unique <- sort(gene_list_unique, decreasing = TRUE)

# 2. Run GSEA (Gene Set Enrichment Analysis)
gsea.go <- gseGO(
  geneList = gene_list_unique,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE,
  by = "fgsea" # Added for faster permutations
)

# 3. Plot with ridgeplot()
plot_gsea_go <- ridgeplot(gsea.go, showCategory = 15, fill = "pvalue") +
  ggtitle("GSEA Ridgeplot: GO Biological Processes")
ggsave("GSEA_GO_BP_Ridgeplot.png", plot = plot_gsea_go, device = "png", width = 10, height = 7)

# also perform KEGG GSEA
gsea.kegg <- gseKEGG(
  geneList      = gene_list_unique,
  organism      = "hsa",     # KEGG code for human
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  verbose       = FALSE,
  by = "fgsea" # Added for faster permutations
)
# plot
plot_gsea_kegg <- ridgeplot(gsea.kegg, showCategory = 15, fill = "pvalue") +
  ggtitle("GSEA Ridgeplot - KEGG Pathways")
ggsave("GSEA_KEGG_Ridgeplot.png", plot = plot_gsea_kegg, device = "png", width = 10, height = 7)

# use gseaplot2() to zoom in on a pathway of interest:
plot_gsea_kegg_zoom <- gseaplot2(gsea.kegg, geneSetID = gsea.kegg$ID[1], pvalue_table = TRUE)
ggsave("GSEA_KEGG_Pathway_Zoom.png", plot = plot_gsea_kegg_zoom, device = "png", width = 10, height = 8)

# Reactome GSEA
gsea.reactome <- gsePathway(
  geneList      = gene_list_unique,
  organism      = "human",
  pvalueCutoff  = 0.05,
  minGSSize     = 10,
  maxGSSize     = 500,
  verbose       = FALSE
)

# Plot
plot_gsea_reactome <- ridgeplot(gsea.reactome, showCategory = 15, fill = "pvalue") +
  ggtitle("GSEA Ridgeplot - Reactome Pathways")
ggsave("GSEA_Reactome_Ridgeplot.png", plot = plot_gsea_reactome, device = "png", width = 10, height = 7)

# Merge the 3 plots:
# 1. Create the individual ridgeplots
p.go <- ridgeplot(gsea.go, showCategory = 15, fill = "pvalue") +
  ggtitle("GSEA Ridgeplot: GO Biological Processes") +
  theme(plot.title = element_text(face = "bold", size = 14))

p.kegg <- ridgeplot(gsea.kegg, showCategory = 15, fill = "pvalue") +
  ggtitle("GSEA Ridgeplot: KEGG Pathways") +
  theme(plot.title = element_text(face = "bold", size = 14))

p.reactome <- ridgeplot(gsea.reactome, showCategory = 15, fill = "pvalue") +
  ggtitle("GSEA Ridgeplot: Reactome Pathways") +
  theme(plot.title = element_text(face = "bold", size = 14))

# 2. Combine them using patchwork with guides = "collect"
combined_gsea_plot <- p.go / p.kegg / p.reactome +
  plot_annotation(title = "GSEA Ridgeplots for GO, KEGG, and Reactome",
                  theme = theme(plot.title = element_text(size = 16, face = "bold"))) +
  plot_layout(guides = "collect") # Added guides = "collect"

# 3. Save to high-resolution output as PNG
ggsave("GSEA_ridgeplots_combined.png", combined_gsea_plot, device = "png", width = 10, height = 18)





# End of the script
