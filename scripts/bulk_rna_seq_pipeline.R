setwd("E:/Bulk RNA project/projects/bulk_rnaseq/")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install core packages
BiocManager::install(c("tximport", "DESeq2", "readr", "biomaRt"))

# Load them
library(tximport)
library(DESeq2)
library(readr)
library(biomaRt)
a

samples <- read.csv("samples.csv")
samples



BiocManager::install("GenomicFeatures")

BiocManager::install("txdbmaker")


library(GenomicFeatures)


# Create a TxDb from your GTF file
txdb <- makeTxDbFromGFF("ref_genome/Danio_rerio.GRCz11.111.gtf")

# Extract transcript-to-gene mapping
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")

# Preview
head(tx2gene)


# Check your samples again
print(samples)

# Add path to each quant.sf file based on sample names
files <- file.path("results/salmon_quant", samples$sample, "quant.sf")
names(files) <- samples$sample

# Quick check
head(files)


# Run tximport to load quant files (This loads transcript-level quantifications and summarizes them to gene level.)
library(tximport)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

head(txi$counts)

# Reformat condition to group replicates like this (Prepare Sample Metadata)
samples$group <- gsub("_[12]", "", samples$condition)
samples$group <- factor(samples$group, levels = c("0dpa", "1dpa", "4dpa"))  # set baseline


# Create DESeq2 Dataset Object
library(DESeq2)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ group)

# Filter Lowly Expressed Genes (This improves accuracy and reduces noise)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]


# Run the Differential Expression Pipeline (This runs normalization, dispersion estimation, and model fitting)
dds <- DESeq(dds)


# Extract and Save Results
# This will extract results for pairwise comparisons:

# 1dpa vs 0dpa
res_1_vs_0 <- results(dds, contrast = c("group", "1dpa", "0dpa"))

# 4dpa vs 1dpa
res_4_vs_1 <- results(dds, contrast = c("group", "4dpa", "1dpa"))

# 4dpa vs 0dpa
res_4_vs_0 <- results(dds, contrast = c("group", "4dpa", "0dpa"))

# View summary
summary(res_1_vs_0)
summary(res_4_vs_0)
summary(res_4_vs_1)

# To save CSV file
write.csv(as.data.frame(res_1_vs_0), "results/deseq2_1dpa_vs_0dpa.csv")
write.csv(as.data.frame(res_4_vs_1), "results/deseq2_4dpa_vs_1dpa.csv")
write.csv(as.data.frame(res_4_vs_0), "results/deseq2_4dpa_vs_0dpa.csv")


# MA Plot and QC
plotMA(res_1_vs_0, main = "1dpa vs 0dpa")
plotMA(res_4_vs_1, main = "4dpa vs 1dpa")
plotMA(res_4_vs_0, main = "4dpa vs 0dpa")



## Annotate DEGs with gene symbols and descriptions
# Load biomaRt
library(biomaRt)

# Connect to Ensembl for zebrafish
ensembl <- useMart("ensembl", dataset = "drerio_gene_ensembl")

# Extract unique gene IDs from one result (they’re the same across comparisons)
gene_ids <- rownames(res_1_vs_0)

# Query gene symbol and name
annotation <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = ensembl
)

# Merge with DE results
res_1_vs_0_annot <- merge(as.data.frame(res_1_vs_0), annotation, by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE)
res_4_vs_1_annot <- merge(as.data.frame(res_4_vs_1), annotation, by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE)
res_4_vs_0_annot <- merge(as.data.frame(res_4_vs_0), annotation, by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE)

# Save annotated versions
write.csv(res_1_vs_0_annot, "results/deseq2_1dpa_vs_0dpa_annotated.csv", row.names = FALSE)
write.csv(res_4_vs_1_annot, "results/deseq2_4dpa_vs_1dpa_annotated.csv", row.names = FALSE)
write.csv(res_4_vs_0_annot, "results/deseq2_4dpa_vs_0dpa_annotated.csv", row.names = FALSE)



## Functional Enrichment (GO, KEGG)

# Install if not already
BiocManager::install(c("clusterProfiler", "org.Dr.eg.db"))

# Load
library(clusterProfiler)
library(org.Dr.eg.db)

# Extract significant genes (padj < 0.05)
run_enrichment <- function(res_df, comparison_label) {
  sig_genes <- rownames(res_df[which(res_df$padj < 0.05), ])
  entrez_ids <- mapIds(org.Dr.eg.db, keys = sig_genes,
                       column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  
  # Remove NAs
  entrez_ids <- na.omit(entrez_ids)
  
  # GO
  go <- enrichGO(gene = entrez_ids, OrgDb = org.Dr.eg.db,
                 ont = "BP", pAdjustMethod = "BH", readable = TRUE)
  png(paste0("results/plots/GO_dotplot_", comparison_label, ".png"), width = 900, height = 600)
  print(dotplot(go, showCategory = 20, title = paste("GO BP -", comparison_label)))
  dev.off()
  
  # KEGG
  kegg <- enrichKEGG(gene = entrez_ids, organism = "dre")
  png(paste0("results/plots/KEGG_dotplot_", comparison_label, ".png"), width = 900, height = 600)
  print(dotplot(kegg, showCategory = 20, title = paste("KEGG -", comparison_label)))
  dev.off()
  
  return(list(GO = go, KEGG = kegg))
}

# Save the graphs showing enrichment
dir.create("results/plots", showWarnings = FALSE)

enrich_1_vs_0 <- run_enrichment(res_1_vs_0, "1dpa_vs_0dpa")
enrich_4_vs_1 <- run_enrichment(res_4_vs_1, "4dpa_vs_1dpa")
enrich_4_vs_0 <- run_enrichment(res_4_vs_0, "4dpa_vs_0dpa")



# Save enrichment tables
write.csv(enrich_1_vs_0$GO@result, "results/GO_enrichment_1dpa_vs_0dpa.csv")
write.csv(enrich_1_vs_0$KEGG@result, "results/KEGG_enrichment_1dpa_vs_0dpa.csv")

write.csv(enrich_4_vs_1$GO@result, "results/GO_enrichment_4dpa_vs_1dpa.csv")
write.csv(enrich_4_vs_1$KEGG@result, "results/KEGG_enrichment_4dpa_vs_1dpa.csv")

write.csv(enrich_4_vs_0$GO@result, "results/GO_enrichment_4dpa_vs_0dpa.csv")
write.csv(enrich_4_vs_0$KEGG@result, "results/KEGG_enrichment_4dpa_vs_0dpa.csv")



## Convert ENSEMBL gene IDs to gene symbols using org.Dr.eg.db

library(org.Dr.eg.db)

# Convert ENSEMBL gene IDs (rownames of mat) to SYMBOL
gene_symbols <- mapIds(org.Dr.eg.db,
                       keys = rownames(mat),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Replace ENSEMBL IDs with gene symbols
rownames(mat) <- ifelse(is.na(gene_symbols), rownames(mat), gene_symbols)




## 📌 Volcano Plots for All Comparisons 

# Required libraries
if (!requireNamespace("EnhancedVolcano", quietly = TRUE))
  BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)
library(org.Dr.eg.db)

# Create output directory
dir.create("results/plots/volcano", recursive = TRUE, showWarnings = FALSE)

# Volcano plot function with gene symbol labels
plot_volcano <- function(res_df, comparison_label) {
  
  # Drop rows with NA in log2FoldChange or padj
  res_df <- res_df[complete.cases(res_df[, c("log2FoldChange", "padj")]), ]
  
  # Convert ENSEMBL IDs to gene symbols
  gene_symbols <- mapIds(org.Dr.eg.db,
                         keys = rownames(res_df),
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  
  # Replace NA symbols with ENSEMBL fallback
  gene_labels <- ifelse(is.na(gene_symbols), rownames(res_df), gene_symbols)
  
  # Plot
  png(paste0("results/plots/volcano/volcano_", comparison_label, ".png"), width = 1000, height = 800)
  print(
    EnhancedVolcano(res_df,
                    lab = gene_labels,
                    x = "log2FoldChange",
                    y = "padj",
                    title = paste("Volcano Plot:", comparison_label),
                    pCutoff = 0.05,
                    FCcutoff = 1.5,
                    pointSize = 2,
                    labSize = 3,
                    colAlpha = 0.6)
  )
  dev.off()
  
  # Optional: save symbol-mapped results
  output_df <- res_df
  output_df$GeneSymbol <- gene_symbols
  write.csv(output_df, file = paste0("results/deseq2_", comparison_label, "_with_symbols.csv"))
}

# 🧨 Run volcano plots for all comparisons
plot_volcano(res_1_vs_0, "1dpa_vs_0dpa")
plot_volcano(res_4_vs_1, "4dpa_vs_1dpa")
plot_volcano(res_4_vs_0, "4dpa_vs_0dpa")



## 📈 Pathway Network Visualization: enrichMap (emapplot)
library(enrichplot)

# Create output folder
dir.create("results/plots/enrichMap", showWarnings = FALSE, recursive = TRUE)

# Only run if enough significant GO terms exist
if (nrow(enrich_4_vs_0$GO@result) > 1) {
  go_sim <- pairwise_termsim(enrich_4_vs_0$GO)
  
  png("results/plots/enrichMap/enrichMap_GO_4dpa_vs_0dpa.png", width = 1200, height = 1000)
  emapplot(go_sim,
           showCategory = 30,
           layout = "kk",         # Kamada-Kawai layout (for clean separation)
           color = "p.adjust")    # Color nodes by adjusted p-value
  dev.off()
}



# Reset
while (!is.null(dev.list())) dev.off()

# Ensure folder exists
dir.create("results/plots/enrichMap", showWarnings = FALSE, recursive = TRUE)

go_sim <- pairwise_termsim(enrich_4_vs_0$GO)
png("results/plots/enrichMap/enrichMap_GO_4dpa_vs_0dpa.png", width = 1200, height = 1000)
emapplot(go_sim, showCategory = 30, layout = "kk", color = "p.adjust")
dev.off()



## 📈 Automate enrichMap (emapplot) for All Comparisons

# Load library
library(enrichplot)

# Create output folder
dir.create("results/plots/enrichMap", showWarnings = FALSE, recursive = TRUE)

# Function to generate enrichMap plot
plot_enrichMap <- function(go_obj, label) {
  while (!is.null(dev.list())) dev.off()  # Ensure no open graphics devices
  
  # Only plot if sufficient GO terms exist
  if (nrow(go_obj@result) > 1) {
    go_sim <- pairwise_termsim(go_obj)
    
    tryCatch({
      png(paste0("results/plots/enrichMap/enrichMap_GO_", label, ".png"), width = 1200, height = 1000)
      print(
      emapplot(go_sim,
               showCategory = 30,
               layout = "kk",
               color = "p.adjust")
      )
      dev.off()
    }, error = function(e) {
      message(paste("⚠️ Skipped:", label, "- Reason:", e$message))
    })
  } else {
    message(paste("❗ Not enough GO terms for:", label))
  }
}

# Run for each comparison
plot_enrichMap(enrich_1_vs_0$GO, "1dpa_vs_0dpa")
plot_enrichMap(enrich_4_vs_1$GO, "4dpa_vs_1dpa")
plot_enrichMap(enrich_4_vs_0$GO, "4dpa_vs_0dpa")



## 📌 Venn Diagrams of DEGs (padj < 0.05)
# Libraries
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
if (!requireNamespace("ggvenn", quietly = TRUE)) install.packages("ggvenn")

library(VennDiagram)
library(ggvenn)
library(dplyr)

# Extract DEG lists
deg_1_0 <- rownames(subset(res_1_vs_0, padj < 0.05))
deg_4_1 <- rownames(subset(res_4_vs_1, padj < 0.05))
deg_4_0 <- rownames(subset(res_4_vs_0, padj < 0.05))

# Save DEG sets to CSV
deg_list_df <- list(
  `1dpa_vs_0dpa` = deg_1_0,
  `4dpa_vs_1dpa` = deg_4_1,
  `4dpa_vs_0dpa` = deg_4_0
)
# Save intersection sets to separate CSV
dir.create("results/venn", recursive = TRUE, showWarnings = FALSE)

# Convert to data frame with list column for saving
write.csv(data.frame(
  comparison = rep(names(deg_list_df), sapply(deg_list_df, length)),
  gene_id = unlist(deg_list_df)
), file = "results/venn/DEG_lists_all_comparisons.csv", row.names = FALSE)


venn_genes <- list(
  "1dpa_vs_0dpa" = deg_1_0,
  "4dpa_vs_1dpa" = deg_4_1,
  "4dpa_vs_0dpa" = deg_4_0
)

# Save intersecting sets
common_1_4 <- intersect(deg_1_0, deg_4_0)
common_all <- Reduce(intersect, venn_genes)

write.csv(data.frame(GeneID = common_1_4), "results/venn/DEGs_shared_1dpa_4dpa.csv", row.names = FALSE)
write.csv(data.frame(GeneID = common_all), "results/venn/DEGs_shared_all.csv", row.names = FALSE)

# Plot using base VennDiagram
png("results/plots/venn_deseq2.png", width = 1000, height = 800)
venn.plot <- venn.diagram(
  x = venn_genes,
  filename = NULL,
  fill = c("cornflowerblue", "seagreen3", "orchid"),
  alpha = 0.5,
  cat.cex = 1.2,
  main = "Overlap of DEGs (padj < 0.05)"
)
grid.draw(venn.plot)
dev.off()

# Alternative: ggvenn
png("results/plots/ggvenn_deseq2.png", width = 1000, height = 800)
ggvenn(venn_genes, fill_color = c("steelblue", "darkgreen", "orchid"))
dev.off()

## To save DEG table with gene symbol
# Load required library
library(org.Dr.eg.db)

# Create directory if not already
dir.create("results/venn", recursive = TRUE, showWarnings = FALSE)

# Extract DEG lists (padj < 0.05)
deg_1_0 <- rownames(subset(res_1_vs_0, padj < 0.05))
deg_4_1 <- rownames(subset(res_4_vs_1, padj < 0.05))
deg_4_0 <- rownames(subset(res_4_vs_0, padj < 0.05))

# Combine into a named list
deg_list <- list(
  `1dpa_vs_0dpa` = deg_1_0,
  `4dpa_vs_1dpa` = deg_4_1,
  `4dpa_vs_0dpa` = deg_4_0
)

# Create flat data frame
deg_df <- data.frame(
  comparison = rep(names(deg_list), sapply(deg_list, length)),
  ensembl_id = unlist(deg_list),
  stringsAsFactors = FALSE
)

# Map ENSEMBL to SYMBOL using org.Dr.eg.db
deg_df$gene_symbol <- mapIds(org.Dr.eg.db,
                             keys = deg_df$ensembl_id,
                             column = "SYMBOL",
                             keytype = "ENSEMBL",
                             multiVals = "first")

# Save to CSV
write.csv(deg_df, "results/venn/DEG_lists_all_comparisons_with_symbols.csv", row.names = FALSE)




## 📊 Heatmap of Top 50 DE Genes (4dpa vs 0dpa)

# Step 1: Select top 50 DE genes
res_4_vs_0_ordered <- res_4_vs_0[order(res_4_vs_0$padj), ]
res_4_vs_0_ordered <- res_4_vs_0_ordered[!is.na(res_4_vs_0_ordered$padj), ]
top50_gene_ids <- rownames(res_4_vs_0_ordered)[1:50]

# Step 2: Extract matrix and center
mat <- assay(vsd)[top50_gene_ids, ]
mat <- mat - rowMeans(mat)

# Step 3: Convert ENSEMBL IDs to gene symbols AFTER this step
library(org.Dr.eg.db)
gene_symbols <- mapIds(org.Dr.eg.db,
                       keys = rownames(mat),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
rownames(mat) <- ifelse(is.na(gene_symbols), rownames(mat), gene_symbols)

# Step 4: Create cleaner column labels (e.g., 4dpa_1)
groups <- samples$group
group_counts <- ave(seq_along(groups), groups, FUN = seq_along)
new_colnames <- paste0(groups, "_", group_counts)
colnames(mat) <- new_colnames

# Step 5: Create annotation
annotation <- data.frame(group = groups)
rownames(annotation) <- new_colnames

# Step 6: Plot
while (!is.null(dev.list())) dev.off()
png("results/plots/heatmap_top50_4dpa_vs_0dpa_pretty_symbols.png", width = 1200, height = 1000)
pheatmap(mat,
         annotation_col = annotation,
         show_rownames = TRUE, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Top 50 DE genes (gene symbols): 4dpa vs 0dpa")
dev.off()

# Save matrix
write.csv(mat, "results/heatmap_data_top50_4dpa_vs_0dpa_symbols.csv")




## 🔥 Heatmap for 1dpa vs 0dpa

# 1. Get top 50 by padj
res_1_vs_0_ordered <- res_1_vs_0[order(res_1_vs_0$padj), ]
res_1_vs_0_ordered <- res_1_vs_0_ordered[!is.na(res_1_vs_0_ordered$padj), ]
top50_gene_ids <- rownames(res_1_vs_0_ordered)[1:50]

# 2. Extract matrix
mat <- assay(vsd)[top50_gene_ids, ]
mat <- mat - rowMeans(mat)

# 3. Convert ENSEMBL to gene symbols
library(org.Dr.eg.db)
gene_symbols <- mapIds(org.Dr.eg.db,
                       keys = rownames(mat),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
rownames(mat) <- ifelse(is.na(gene_symbols), rownames(mat), gene_symbols)

# 4. Clean column names
groups <- samples$group
group_counts <- ave(seq_along(groups), groups, FUN = seq_along)
new_colnames <- paste0(groups, "_", group_counts)
colnames(mat) <- new_colnames

# 5. Annotation
annotation <- data.frame(group = groups)
rownames(annotation) <- new_colnames

# 6. Plot and save
while (!is.null(dev.list())) dev.off()
png("results/plots/heatmap_top50_1dpa_vs_0dpa.png", width = 1200, height = 1000)
pheatmap(mat,
         annotation_col = annotation,
         show_rownames = TRUE, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Top 50 DE genes: 1dpa vs 0dpa")
dev.off()

# 7. Save data
write.csv(mat, "results/heatmap_data_top50_1dpa_vs_0dpa.csv")



## 🔥 Heatmap for 4dpa vs 1dpa

# 1. Get top 50 by padj
res_4_vs_1_ordered <- res_4_vs_1[order(res_4_vs_1$padj), ]
res_4_vs_1_ordered <- res_4_vs_1_ordered[!is.na(res_4_vs_1_ordered$padj), ]
top50_gene_ids <- rownames(res_4_vs_1_ordered)[1:50]

# 2. Extract matrix
mat <- assay(vsd)[top50_gene_ids, ]
mat <- mat - rowMeans(mat)

# 3. Convert ENSEMBL to gene symbols
gene_symbols <- mapIds(org.Dr.eg.db,
                       keys = rownames(mat),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
rownames(mat) <- ifelse(is.na(gene_symbols), rownames(mat), gene_symbols)

# 4. Clean column names
groups <- samples$group
group_counts <- ave(seq_along(groups), groups, FUN = seq_along)
new_colnames <- paste0(groups, "_", group_counts)
colnames(mat) <- new_colnames

# 5. Annotation
annotation <- data.frame(group = groups)
rownames(annotation) <- new_colnames

# 6. Plot and save
while (!is.null(dev.list())) dev.off()
png("results/plots/heatmap_top50_4dpa_vs_1dpa.png", width = 1200, height = 1000)
pheatmap(mat,
         annotation_col = annotation,
         show_rownames = TRUE, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Top 50 DE genes: 4dpa vs 1dpa")
dev.off()

# 7. Save data
write.csv(mat, "results/heatmap_data_top50_4dpa_vs_1dpa.csv")



data.frame(
  SRR_ID = samples$sample,
  Group = samples$group,
  New_Label = new_colnames
)


## Extract expression values of GOI
# Get normalized counts from DESeq2 (VST)
polI_genes <- c("rrn3", "taf1b", "polr1a", "polr1b", "polr1c", "polr1d", "polr1e", "polr1f", "polr1g")

# If rownames are Ensembl IDs, convert first or match symbols
vsd_mat <- assay(vsd)
gene_symbols <- mapIds(org.Dr.eg.db, keys = rownames(vsd_mat),
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
rownames(vsd_mat) <- ifelse(is.na(gene_symbols), rownames(vsd_mat), gene_symbols)

# Subset expression
polI_expr <- vsd_mat[rownames(vsd_mat) %in% polI_genes, ]
write.csv(polI_expr, "results/polI_gene_expression.csv")


## Barplot or Heatmap of GOI
library(pheatmap)
pheatmap(polI_expr, cluster_rows = FALSE, cluster_cols = TRUE,
         annotation_col = data.frame(group = samples$group),
         main = "RNA Pol I–Related Gene Expression")


# For 4dpa vs 0dpa
res_subset <- res_4_vs_0[rownames(res_4_vs_0) %in% rownames(polI_expr), ]
write.csv(as.data.frame(res_subset), "results/polI_DE_4dpa_vs_0dpa.csv")


library(EnhancedVolcano)
library(org.Dr.eg.db)

EnhancedVolcano(res_4_vs_0,
                lab = rownames(res_4_vs_0),
                selectLab = c("rrn3", "taf1b", "polr1b", "polr1e"),
                x = 'log2FoldChange',
                y = 'padj',
                title = '4dpa vs 0dpa',
                pCutoff = 0.05,
                FCcutoff = 1.5)

