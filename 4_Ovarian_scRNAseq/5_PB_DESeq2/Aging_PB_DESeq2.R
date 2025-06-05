setwd("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/For_github/PB_DESeq2_GO")

# Load necessary libraries
library('Seurat')
library('muscat')
library('sctransform')
library('DESeq2')
library('ggplot2')
library('ComplexHeatmap')
library('pheatmap')
library('parallel')
library(RColorBrewer)

theme_set(theme_bw())

rm(list = ls())

################################################################################
# Perform PB
# Aggregate data & run DESeq2
# Aging model
################################################################################

################################################################################
# 1. Load data
################################################################################

load("/Volumes/jinho01/Benayoun_lab/Projects/Menopause_model_project/Data/Ovarian_scRNAseq/0_Annotated_Seurat_objects/Without_scTE/2024-10-15_10x_ovary_Benayoun_lab_AC_Seurat_object_with_final_annotation.RData")

DefaultAssay(ovary.AC) <- "RNA"

# Define directories
output_dirs <- list(
  results = file.path("./DESeq2_results"),
  summary = file.path("./Summary")
)

# Create directories 
for (dir in output_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

################################################################################
# 2. Aggregate counts by cell type
################################################################################

# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
sce <- as.SingleCellExperiment(ovary.AC)

rm(ovary.AC)

# Prepare data for aggregation
sce <- prepSCE(sce, 
               kid = "celltype.level2",  # Cell type variable
               gid = "Age",              # Group variable (e.g., condition)
               sid = "Library",          # Sample variable
               drop = TRUE)              # Drop other colData columns

# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(sce$cluster_id))
ns  <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
# t(table(sce$cluster_id, sce$sample_id))

# Aggregation of single-cell to pseudobulk data
pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# assayNames(pb)

# nb. of cells per cluster-sample
cell.per.samp.tab <- t(table(sce$cluster_id, sce$sample_id))

# cell types with at least 25 cells from every each sex/cohort sample
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab > 25, 2, sum) == 4]
# keep:
# "B"          "BEC"        "DNT"        "Epithelial" "Granulosa"  "LEC"        "Mesenchyme" "NK"         "Stroma"     "Theca"     

### extract pseudobulk information for samples that pass the cell number cutoff
counts.pb <- pb@assays@data[celltype.qc]

################################################################################
# 3. Perform pseudobulk DESeq2 analysis for each cell type
################################################################################

# Prepare lists to store results
vst.counts <- vector("list", length(counts.pb))
names(vst.counts) <- names(counts.pb)

deseq.res.list <- vector("list", length(counts.pb))
names(deseq.res.list) <- names(counts.pb)

significant_gene_counts <- list() 

for (cell_type in names(counts.pb)) {
  
  my.outprefix <- paste0(Sys.Date(), "_Aging_PB_DEseq2_", cell_type)
  
  # Filter counts for current cell type
  counts <- counts.pb[[cell_type]]
  
  # Filter genes expressed in at least 3 samples
  good_genes <- rowSums(counts > 0) >= 3
  counts <- counts[good_genes, ]
  
  # Create DESeq2 dataset
  dataDesign <- data.frame(
    row.names = colnames(counts),
    age = ifelse(grepl("YF", colnames(counts)), "YF", "OF")
  )
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = dataDesign,
    design = ~ age
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract DESeq2 results
  res <- results(dds, contrast = c("age", "OF", "YF"))
  res <- res[!is.na(res$padj), ]
  
  deseq.res.list[[cell_type]] <- res
  
  # Get vst data
  vst_data <- getVarianceStabilizedData(dds)
  vst.counts[[cell_type]] <- vst_data
  
  # Save significant genes (FDR < 0.05)
  significant_genes <- rownames(res)[res$padj < 0.05]
  significant_gene_counts[[cell_type]] <- length(significant_genes)  
 
}

# Save
save(deseq.res.list, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Aging_PB_DESeq2_object.RData")))
save(vst.counts, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Aging_PB_VST_counts.RData")))

######################################
# 4. Plot strip plot
######################################

set.seed(123123123)

# Set order of cell types
cell_type_order <- c(
  "Granulosa", "Theca", "Stroma", "Mesenchyme", "Pericyte", 
  "BEC", "LEC", "Epithelial", "Neutrophil", "Monocyte", 
  "Macrophage", "DC", "ILC", "NK", "NKT", 
  "CD8 NKT", "CD8 T", "CD4 T", "DNT", "DPT", "B"
)

# Filter and reorder results
available_cell_types <- intersect(cell_type_order, names(deseq.res.list))
deseq.res.list <- deseq.res.list[match(available_cell_types, names(deseq.res.list))] 

# Initialize axis labels with significant gene counts
xlab <- character(length = length(deseq.res.list))
for (i in seq_along(deseq.res.list)) {
  sig_genes <- sum(deseq.res.list[[i]]$padj < 0.05, na.rm = TRUE)  # Count significant genes
  xlab[i] <- paste(available_cell_types[i], "\n(", sig_genes, " sig.)", sep = "")
}

# Generate the strip plot
pdf(file.path(output_dirs$summary, paste0(Sys.Date(), "_Aging_stripplot_YF_vs_OF_significant_genes.pdf")), width = 6, height = 5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))

# Initialize an empty plot
plot(
  x = 1, y = 1,
  type = "n",
  xlim = c(0.5, length(deseq.res.list) + 0.5),
  ylim = c(-15, 20),
  axes = FALSE,
  xlab = "", ylab = "Log2 fold change (OF / YF)"
)

# Add horizontal reference lines
abline(h = 0)  # Line at 0
abline(h = seq(-20, 20, by = 5)[-5], lty = "dotted", col = "grey")

# Loop through each cell type and plot jittered points
for (i in seq_along(deseq.res.list)) {
  # Extract data for the current cell type
  current_result <- deseq.res.list[[i]]
  
  # Skip if the result is empty or has no rows
  if (is.null(current_result) || nrow(current_result) == 0) next
  
  # Identify significant genes (FDR < 0.05)
  sig_genes <- current_result$padj < 0.05
  
  # Assign colors based on log2 fold change direction for significant genes
  colors <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), nrow(current_result))  # Default grey
  colors[sig_genes & current_result$log2FoldChange > 0] <- "deeppink4"  # Up in OF
  colors[sig_genes & current_result$log2FoldChange < 0] <- "deeppink1"  # Up in YF
  
  # Plot the jittered points with the assigned colors
  points(
    x = jitter(rep(i, nrow(current_result)), amount = 0.2),
    y = rev(current_result$log2FoldChange),  # Plot log2 fold change as-is
    pch = 16, col = rev(colors), cex = 0.5, bg = rev(colors)
  )
}

# Add x-axis labels with the number of significant genes
axis(1, at = 1:length(deseq.res.list), tick = FALSE, las = 2, lwd = 0,
     labels = xlab, cex.axis = 0.7)

# Add y-axis labels
axis(2, las = 1, at = seq(-15, 20, 5))

box()
dev.off()

################################################################################
sink(file = paste(Sys.Date(), "_Aging_PB_DESeq2_analysis_session_Info.txt"))
sessionInfo()
sink()
