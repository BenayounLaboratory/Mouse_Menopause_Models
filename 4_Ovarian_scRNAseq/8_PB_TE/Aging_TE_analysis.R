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
# Aging model
################################################################################

################################################################################
# 1. Load data
################################################################################

load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/8_PB_TE/Input/2023-08-30_10x_ovary_Benayoun_lab_AC_Seurat_object_scTE_combined_FINAL.RData")

DefaultAssay(ovary.AC.genes.TEs.combined.seurat) <- "RNA"

# Define directories
output_dirs <- list(
  results = file.path("./DESeq2_results")
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
sce <- as.SingleCellExperiment(ovary.AC.genes.TEs.combined.seurat)

rm(ovary.AC.genes.TEs.combined.seurat)

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
  
  my.outprefix <- paste0(Sys.Date(), "_Aging_PB_with_scTE_DEseq2_", cell_type)
  
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
save(deseq.res.list, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Aging_PB_with_scTE_DESeq2_object.RData")))
save(vst.counts, file = file.path(output_dirs$results, paste0(Sys.Date(), "_Aging_PB_with_scTE_VST_counts.RData")))


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
  
  pdf_filename <- paste(Sys.Date(), "AC_GSEA_TE_BUBBLE_plot.pdf", sep = "_")
  ggsave(pdf_filename, plot = my.plot, width = 10, height = 6)
  
  return(my.plot)
}

# Generate and save bubble plot
gsea_bubble_plot <- get_TE_GSEA_bubble_plot(gsea_summary_df)
print(gsea_bubble_plot)

################################################################################
sink(file = paste(Sys.Date(), "_Aging_TE_analysis_session_Info.txt"))
sessionInfo()
sink()
