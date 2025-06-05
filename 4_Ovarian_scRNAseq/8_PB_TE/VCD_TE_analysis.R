setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/8_PB_TE")

# Load libraries
library(Seurat)
library(sctransform)
library(muscat)
library(DESeq2)
library(clusterProfiler)
library(GSEABase)
library(ggplot2)
library(ComplexHeatmap)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(dplyr)
library(tibble)
library(tidyverse)
library(parallel)

theme_set(theme_bw())

rm(list = ls())

################################################################################
# Aggregate data & run DESeq2
# Assess TE expression via GSEA
# VCD model
################################################################################

################################################################################
# 1. Load data
################################################################################

# Load Seurat object
load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/8_PB_TE/Input/2024-03-01_10x_ovary_Benayoun_lab_VCD_Seurat_object_scTE_combined_FINAL.RData")

DefaultAssay(my.ovary.VCD.genes.TEs.combined.seurat) <- "RNA"

# Extract metadata

# Subset only the relevant columns for pseudobulk analysis
metadata <- my.ovary.VCD.genes.TEs.combined.seurat@meta.data %>%
  dplyr::select(Library, Age, Treatment, Duration, Batch) %>%
  distinct()

# Rename columns to match the required format
colnames(metadata) <- c("sample_id", "Age", "Treatment", "Duration", "Batch")

# Set rownames to match sample IDs
rownames(metadata) <- metadata$sample_id

# Define output directories
output_dirs <- list(
  results = file.path("./DESeq2_results")
)

# Create directories
for (dir in output_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

################################################################################
# 2. Aggregate Counts by Cell Type
################################################################################

# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
sce <- as.SingleCellExperiment(my.ovary.VCD.genes.TEs.combined.seurat)

rm(my.ovary.VCD.genes.TEs.combined.seurat)

# Prepare data for aggregation
sce <- prepSCE(sce, 
               kid = "celltype.level2",  # Cell type variable
               gid = "Treatment",        # Group variable (e.g., condition)
               sid = "Library"  ,        # Sample variable
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
#   [1] "B"          "BEC"        "CD4 T"      "CD8 NKT"    "DC"         "DNT"        "DPT"        "Epithelial" "Granulosa"  "ILC"        "LEC"        "Macrophage" "Mesenchyme" "Monocyte"   "Neutrophil" "NK"         "NKT"        "Pericyte"   "Stroma"    
#  [20] "Theca" 

# nb. of cells per cluster-sample
cell.per.samp.tab <- t(table(sce$cluster_id, sce$sample_id))

# cell types with at least 25 cells from at least 12 samples (75%)
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab > 25, 2, sum) >= 12]
# keep:
# [1] "BEC"        "DNT"        "Epithelial" "Granulosa"  "LEC"        "Stroma"     "Theca"     

### extract pseudobulk information for samples that pass the cell number cutoff
counts.pb <- pb@assays@data[celltype.qc]

################################################################################
# 3. Perform pseudobulk DESeq2 analysis for each cell type
################################################################################

# Prepare lists
sva.cleaned.counts <- vector(mode = "list", length = length(counts.pb))
names(sva.cleaned.counts) <- names(counts.pb)

vst.counts <- vector("list", length(counts.pb))
names(vst.counts) <- names(counts.pb)

deseq.res.list <- vector("list", length(counts.pb))
names(deseq.res.list) <- names(counts.pb)

significant_gene_counts <- list() 

for (cell_type in names(counts.pb)) {
  
  my.outprefix <- paste0(Sys.Date(), "_VCD_PB_with_scTE_DEseq2_", cell_type)
  
  # Extract count matrix
  counts <- counts.pb[[cell_type]]
  
  # Filter genes expressed in at least 14 samples
  good_genes <- rowSums(counts > 0) >= 14
  counts <- counts[good_genes, ]
  
  # Define sample metadata
  dataDesign <- data.frame(
    row.names = colnames(counts),
    age = ifelse(grepl("3m", colnames(counts)), "3m", "10m"),
    treatment = ifelse(grepl("CTL", colnames(counts)), "CTL", "VCD"),
    duration = ifelse(grepl("30d", colnames(counts)), "30d", "90d"),
    batch = metadata[colnames(counts), "Batch"]
  )
  
  # Build design matrix (null and alternative models)
  mod1 <- model.matrix(~ age + treatment + duration + batch, data = dataDesign)
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
  
  # Save SVA-corrected counts
  write.table(my.filtered.sva, file.path(output_dirs$results, paste0(my.outprefix, "_SVA_Corrected_Counts.txt")), sep = "\t", quote = FALSE)
  
  dds <- DESeqDataSetFromMatrix(
    countData = sva.cleaned.counts[[cell_type]],
    colData = dataDesign,
    design = ~ age + treatment + duration
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract results
  res <- results(
    dds, 
    contrast = c("treatment", "VCD", "CTL")
  )
  res <- res[!is.na(res$padj), ]
  
  # Get vst data
  vst_data <- getVarianceStabilizedData(dds)
  
  # Store results
  deseq.res.list[[cell_type]] <- res
  vst.counts[[cell_type]] <- vst_data
  
  # Save significant genes (FDR < 0.05)
  significant_genes <- rownames(res)[res$padj < 0.05]
  significant_gene_counts[[cell_type]] <- length(significant_genes)  
   
}

# Save
save(vst.counts, file = file.path(output_dirs$results, paste0(Sys.Date(), "_VCD_PB_with_scTE_VST_counts.RData")))
save(deseq.res.list, file = file.path(output_dirs$results, paste0(Sys.Date(), "_VCD_PB_with_scTE_DESeq2_object.RData")))

################################################################################
# 4. Assess TE expression
################################################################################

# Load TE family annotation file
my.TEs <- read.table("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/8_PB_TE/Input/rmsk.txt", header = FALSE)

# Clean TE family names (remove '?' at the end)
my.TEs <- my.TEs %>%
  mutate(V12 = str_replace(V12, "\\?$", "")) %>%  # Remove '?' if at the end
  distinct(V12, V11) %>% 
  rename(term = V12, gene = V11) 

# Prepare ranked gene lists for GSEA
gsea_geneLists <- lapply(deseq.res.list, function(res) {
  gene_names <- rownames(res)
  gene_stat <- res$stat[!is.na(res$stat)]
  gene_list <- setNames(gene_stat, gsub("^(scTE-|scTE_TE-)", "", gene_names))
  sort(gene_list, decreasing = TRUE)
})

# Run GSEA for each cell type
gsea_results_list <- lapply(gsea_geneLists, function(gene_list) {
  GSEA(geneList = gene_list, TERM2GENE = my.TEs, pvalueCutoff = 1)
})

# Extract and process GSEA results
gsea_summary_df <- do.call(rbind, lapply(names(gsea_results_list), function(cell_type) {
  gsea_results <- gsea_results_list[[cell_type]]
  if (!is.null(gsea_results) && length(gsea_results) > 0) {
    gsea_df <- as.data.frame(gsea_results)
    gsea_df %>%
      select(ID, NES, p.adjust) %>%
      rename(TE_Family = ID, padj = p.adjust) %>%
      mutate(Cell_Type = cell_type)
  } else {
    NULL
  }
}))

# Define TE family order
te_family_order <- c("LINE", "SINE", "LTR", "DNA")
gsea_summary_df$TE_Family <- factor(gsea_summary_df$TE_Family, levels = te_family_order, ordered = TRUE)

# Define cell type order
desired_order <- c("Granulosa", "Theca", "Stroma", "Mesenchyme", "Pericyte",
                   "BEC", "LEC", "Epithelial", "Neutrophil", "Monocyte",
                   "Macrophage", "DC", "ILC", "NK", "NKT", "CD8 NKT",
                   "CD8 T", "CD4 T", "DNT", "DPT", "B")

gsea_summary_df$Cell_Type <- factor(gsea_summary_df$Cell_Type, levels = desired_order, ordered = TRUE)

# Generate bubble plot
get_TE_GSEA_bubble_plot <- function(gsea_summary_df) {
  
  gsea_summary_df <- gsea_summary_df %>%
    mutate(
      minlog10fdr = -log10(padj + 1e-30),
      dot_alpha = ifelse(padj <= 0.1, 1, 0.2)
    )
  
  # Define color scale
  my.color.vector.age <- c("darkblue", "dodgerblue4", "dodgerblue3", "dodgerblue1",
                           "white", "lightcoral", "brown1", "firebrick2", "firebrick4")
  
  my.plot <- ggplot(gsea_summary_df, aes(x = Cell_Type, y = TE_Family, colour = NES, size = minlog10fdr, alpha = dot_alpha)) +
    theme_bw() +
    geom_point(shape = 16, na.rm = TRUE) +
    ggtitle("TE GSEA Analysis Across Cell Types") +
    labs(x = "Cell Type", y = "TE Family", size = "-log10(FDR)") +
    scale_x_discrete(limits = desired_order) +
    scale_y_discrete(limits = te_family_order) +
    scale_colour_gradientn(colours = my.color.vector.age, na.value = "grey50", guide = "colourbar") +
    scale_alpha(range = c(0.2, 1), guide = "none")
  
  pdf_filename <- paste(Sys.Date(), "VCD_GSEA_TE_BUBBLE_plot.pdf", sep = "_")
  ggsave(pdf_filename, plot = my.plot, width = 10, height = 6)
  
  return(my.plot)
}

# Generate and save bubble plot
gsea_bubble_plot <- get_TE_GSEA_bubble_plot(gsea_summary_df)
print(gsea_bubble_plot)

################################################################################
# Save session info
sink(file = paste0(Sys.Date(), "_VCD_TE_analysis_Session_Info.txt"))
sessionInfo()
sink()
