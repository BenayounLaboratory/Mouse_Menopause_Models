setwd("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/For_github/PB_DESeq2_GO")

# Load necessary libraries
library(Seurat)
library(SingleCellExperiment)
library(DESeq2)
library(muscat)
library(sva)
library(limma)
library(ggplot2)
library(dplyr)

rm(list = ls())

################################################################################
# Perform PB
# Aggregate data & run DESeq2
# Foxl2 haploinsufficiency model
################################################################################

################################################################################
# 1. Load data
################################################################################

# Load Seurat object
load("/Volumes/jinho01/Benayoun_lab/Projects/Menopause_model_project/Data/Ovarian_scRNAseq/0_Annotated_Seurat_objects/Without_scTE/2024-10-24_10x_ovary_Benayoun_lab_Foxl2_Seurat_object_with_final_annotation.RData")

DefaultAssay(ovary.Foxl2) <- "RNA"

# Extract metadata

# Subset only the relevant columns for pseudobulk analysis
metadata <- ovary.Foxl2@meta.data %>%
  dplyr::select(Library, Age, Genotype, Batch) %>%
  distinct()

# Rename columns
colnames(metadata) <- c("sample_id", "Age", "Genotype", "Batch")

# Set rownames to match sample IDs
rownames(metadata) <- metadata$sample_id

# Define output directories
output_dirs <- list(
  results = file.path("./DESeq2_results"),
  summary = file.path("./Summary")
)

# Create directories
for (dir in output_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

group_colors <- c(
  "young_wt" = "deeppink1",
  "old_wt" = "#B10D66",
  "young_het" = "#ADDFE4",
  "old_het" = "#508B8E"
)

################################################################################
# 2. Aggregate Counts by Cell Type
################################################################################

# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
sce <- as.SingleCellExperiment(ovary.Foxl2)

rm(ovary.Foxl2)

# Prepare data for aggregation
sce <- prepSCE(sce, 
               kid = "celltype.level2",  # Cell type variable
               gid = "Genotype",         # Group variable (e.g., condition)
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

assayNames(pb)
#  [1] "B"          "BEC"        "CD8 NKT"    "DC"         "DNT"        "DPT"        "Epithelial" "Granulosa"  "ILC"        "LEC"        "Macrophage" "Mesenchyme" "Monocyte"   "Neutrophil" "NK"         "NKT"       
# [17] "Pericyte"   "Stroma"     "Theca" 

# nb. of cells per cluster-sample
cell.per.samp.tab <- t(table(sce$cluster_id, sce$sample_id))

# cell types with at least 25 cells from at least 8 samples (80%)
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab > 25, 2, sum) > 8]
# keep:
# [1] "B"          "BEC"        "DNT"        "Epithelial" "Granulosa"  "NK"         "Stroma"     "Theca"     

# Drop B cells: Foxl2_wt_old_2 has 0 cells. Analyze independently
celltype.qc <- celltype.qc[celltype.qc != "B"]

### extract pseudobulk information for samples that pass the cell number cutoff
counts.pb <- pb@assays@data[celltype.qc]

################################################################################
# 3. Perform Pseudobulk DESeq2 Analysis for Each Cell Type
################################################################################

# Prepare lists to store results
sva.cleaned.counts <- vector(mode = "list", length = length(counts.pb))
names(sva.cleaned.counts) <- names(counts.pb)

vst.counts <- vector("list", length(counts.pb))
names(vst.counts) <- names(counts.pb)

deseq.res.list <- vector("list", length(counts.pb))
names(deseq.res.list) <- names(counts.pb)

significant_gene_counts <- list() 

for (cell_type in names(counts.pb)) {
  
  my.outprefix <- paste0(Sys.Date(), "_Foxl2_PB_DEseq2_", cell_type)
  
  # Extract count matrix
  counts <- counts.pb[[cell_type]]
  
  # Filter genes expressed in at least 8 samples
  good_genes <- rowSums(counts > 0) >= 8
  counts <- counts[good_genes, ]
  
  # Define sample metadata
  dataDesign <- data.frame(
    row.names = colnames(counts),
    age = ifelse(grepl("young", colnames(counts)), "young", "old"),
    genotype = ifelse(grepl("wt", colnames(counts)), "wt", "het"),
    batch = metadata[colnames(counts), "Batch"]
  )
  
  # Build design matrix (null and alternative models)
  mod1 <- model.matrix(~ age + genotype + batch, data = dataDesign)
  n.sv <- num.sv(counts, mod1, method = "be")
  
  # Apply SVAseq algorithm
  my.svseq <- svaseq(as.matrix(counts), mod1, n.sv = n.sv, constant = 0.1)
  
  if (my.svseq$n.sv > 0) {
    # Remove batch effects using surrogate variables
    my.clean <- removeBatchEffect(
      log2(counts + 0.1), 
      batch = dataDesign$batch, 
      covariates = my.svseq$sv,
      design = mod1[,1:3]
    )
  } else {
    my.clean <- removeBatchEffect(
      log2(counts + 0.1), 
      batch = dataDesign$batch, 
      design = mod1[,1:3]
    )
  }
  
  # Convert back to count format for DESeq2
  my.filtered.sva <- round(2^my.clean - 0.1)
  sva.cleaned.counts[[cell_type]] <- my.filtered.sva
  
  dds <- DESeqDataSetFromMatrix(
    countData = sva.cleaned.counts[[cell_type]],
    colData = dataDesign,
    design = ~ age + genotype
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract results
  res <- results(
    dds, 
    contrast = c("genotype", "het", "wt")
  )
  res <- res[!is.na(res$padj), ]
  
  # Get variance-stabilized data
  vst_data <- getVarianceStabilizedData(dds)
  
  # Store results
  deseq.res.list[[cell_type]] <- res
  vst.counts[[cell_type]] <- vst_data
  
  # Save significant genes (FDR < 0.05)
  significant_genes <- rownames(res)[res$padj < 0.05]
  significant_gene_counts[[cell_type]] <- length(significant_genes)  
  
}

################################################################################
# 3-1. Assess B cells
################################################################################

my.outprefix <- paste0(Sys.Date(), "_Foxl2_PB_DEseq2_B_")

# Filter counts for B cells
counts <- pb@assays@data["B"]
counts <- counts$B
counts <- counts[, colnames(counts) != "Foxl2_wt_old_2"]

# Filter genes expressed in at least 8 samples
good_genes <- rowSums(counts > 0) >= 8
counts <- counts[good_genes, ]

# Define sample metadata
dataDesign <- data.frame(
  row.names = colnames(counts),
  age = ifelse(grepl("young", colnames(counts)), "young", "old"),
  genotype = ifelse(grepl("wt", colnames(counts)), "wt", "het"),
  batch = metadata[colnames(counts), "Batch"]
)

# Build design matrix (null and alternative models)
mod1 <- model.matrix(~ age + genotype + batch, data = dataDesign)
n.sv <- num.sv(counts, mod1, method = "be")

# Apply SVAseq algorithm
my.svseq <- svaseq(as.matrix(counts), mod1, n.sv = n.sv, constant = 0.1)

my.clean <- removeBatchEffect(
    log2(counts + 0.1), 
    batch = dataDesign$batch, 
    covariates = my.svseq$sv,
    design = mod1[,1:3]
)

# Convert back to count format for DESeq2
my.filtered.sva <- round(2^my.clean - 0.1)
sva.cleaned.counts[["B"]] <- my.filtered.sva

dds <- DESeqDataSetFromMatrix(
  countData = sva.cleaned.counts[["B"]],
  colData = dataDesign,
  design = ~ age + genotype
)

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(
  dds, 
  contrast = c("genotype", "het", "wt")
)
res <- res[!is.na(res$padj), ]

# Variance-stabilized data
vst_data <- getVarianceStabilizedData(dds)

# Store results
deseq.res.list[["B"]] <- res
vst.counts[["B"]] <- vst_data

# Save significant genes (FDR < 0.05)
significant_genes <- rownames(res)[res$padj < 0.05]
significant_gene_counts[["B"]] <- length(significant_genes)  

# Save
save(vst.counts, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Foxl2_PB_VST_counts.RData")))
save(deseq.res.list, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Foxl2_PB_DESeq2_object.RData")))

################################################################################
# 4. Plot strip plot
################################################################################

set.seed(123123123)

# Set order of cell types
cell_type_order <- c(
  "Granulosa", "Theca", "Stroma", "Mesenchyme", "Pericyte", 
  "BEC", "LEC", "Epithelial", "Neutrophil", "Monocyte", 
  "Macrophage", "DC", "ILC", "NK", "NKT", 
  "CD8 NKT", "CD8 T", "CD4 T", "DNT", "DPT", "B"
)

# Filter and reorder 
available_cell_types <- intersect(cell_type_order, names(deseq.res.list))
deseq.res.list <- deseq.res.list[match(available_cell_types, names(deseq.res.list))] 

# Initialize axis labels with significant gene counts
xlab <- character(length = length(deseq.res.list))
for (i in seq_along(deseq.res.list)) {
  sig_genes <- sum(deseq.res.list[[i]]$padj < 0.05, na.rm = TRUE)
  xlab[i] <- paste(available_cell_types[i], "\n(", sig_genes, " sig.)", sep = "")
}

# Generate the strip plot
pdf(paste0(Sys.Date(), "_Foxl2_stripplot_wt_vs_het_significant_genes.pdf"), width = 6, height = 5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))

# Initialize an empty plot
plot(
  x = 1, y = 1,
  type = "n",
  xlim = c(0.5, length(deseq.res.list) + 0.5),
  ylim = c(-15, 20),
  axes = FALSE,
  xlab = "", ylab = "Log2 fold change (het / wt)"
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
  colors[sig_genes & current_result$log2FoldChange > 0] <- "springgreen"  # Up in het
  colors[sig_genes & current_result$log2FoldChange < 0] <- "deeppink1"    # Up in wt
  
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
# Save session info
sink(file = paste0(Sys.Date(), "_Foxl2_PB_DESeq2_analysis_Session_Info.txt"))
sessionInfo()
sink()
