#!/usr/bin/env Rscript



# Shahryar Alavi
# UCL Institute of Neurology
# 2025-05-08
# s.alavi@ucl.ac.uk



# Below is the part of the analysis used for our ARHGAP19 project to plot patients-healthy control RNA-Seq data. This is assumed that you already have the read counts data from your RNA-Seq.



# Gene enrichment analysis unsing pathfindR:
library(pathfindR)
# Remove NA values from the expression profile:
log2FC.clean <- na.omit(log2fc)
# Remove low confidence genes:
log2FC.clean <- log2FC.clean[log2FC.clean$Pvalue < 0.05, ]

# Run pathFindR
df.path.enrich <- run_pathfindR(log2FC.clean, p_val_threshold = 0.005, output_dir = "pathfindR")
df.processed.genes <- input_processing(log2FC.clean, p_val_threshold = 0.005)



# ==== Gene enrichment plot ====
# Term-Gene Graph
# change hsa IDs to path names:
df.path.enrich2 <- df.path.enrich
df.path.enrich2$ID <- df.path.enrich$Term_Description
# Draw the graph:
term_gene_graph(result_df = df.path.enrich2, num_terms = 3)

# Hierarchical Clustering
clustered.path.enrich <- cluster_enriched_terms(df.path.enrich, plot_dend = FALSE, plot_clusters_graph = FALSE)

# plotting only selected clusters for better visualization
selected_clusters <- subset(clustered.path.enrich, Cluster %in% 1:2)
# Plotting the enrichment bubble chart
enrichment_chart(selected_clusters, plot_by_cluster = TRUE)



# ==== Pathways heatmap ====

# select a subset of clusters:
selected_clusters2 <- subset(clustered.path.enrich, Fold_Enrichment > 2 & Status == "Representative")

# Choose the case sample:
cases <- c("sample names of cases/patients")

score_matrix <- score_terms(
  enrichment_table = selected_clusters2,
  exp_mat = expresion.data,
  cases = cases,
  use_description = TRUE,
  label_samples = FALSE,
  case_title = "ARHGAP19 deficient",
  control_title = "Healthy control",
  low = "#f7797d",
  mid = "#fffde4",
  high = "#1f4037"
)



# ==== Volcano plot ====

# Load the library for volcano plot
library(EnhancedVolcano)

# Define expresion thresholds
keyvals <- ifelse(
  log2FC.clean$Log2FC < -1.5 & log2FC.clean$Pvalue < 0.005, 'red',
  ifelse(log2FC.clean$Log2FC  > 1.5 & log2FC.clean$Pvalue < 0.005, 'darkgreen',
         'gray'))

keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == 'darkgreen'] <- 'Upregulated'
names(keyvals)[keyvals == 'gray'] <- 'Not significant'
names(keyvals)[keyvals == 'red'] <- 'Downregulated'

genes.of.interest <- c("STEAP4", "IL1RL1", "SPRY1", "KCNMB4", "CCL2", "ADH1B", "C11orf87", "ARFGEF3", "KCNA1", "VGF", "KIF20A", "WDR62", "FGFR3", "BIRC5", "FOXM1", "TK1", "IL7R", "ACTC1", "FCRLA", "MCM4", "MCM5", "ORC6", "ANAPC15", "CDK1", "CHEK1", "DYNC1I1", "ITGA6", "ITGA2", "LMNB1")

EnhancedVolcano(log2FC.clean, lab = rownames(log2FC.clean), x = 'Log2FC', y = 'Pvalue', title = 'Differential Expression', pCutoff = 0.005, FCcutoff = 1.5, pointSize = 1, labSize = 5.0, selectLab = genes.of.interest, subtitle = "", legendLabels=c('Not significant','Log2FC','p-value', 'Log2FC and p-value'), legendPosition = 'right', colCustom = keyvals) +
  theme(plot.title = element_text(angle = 0, face = "bold", size = 23, hjust = 0.5))


q()
