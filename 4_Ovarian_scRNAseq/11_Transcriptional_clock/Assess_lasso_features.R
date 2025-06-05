setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/11_Transcriptional_clock")

# Load libraries

library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(ComplexHeatmap)
library(pheatmap)
library(muscat)
library(limma)

rm(list=ls())
set.seed(123)

###############################################################
# Granulosa & Theca transcriptional clocks - assess lasso features
# Perform ORA and generate heatmap (PB data)
###############################################################

###############################################################
# 1. Load data
###############################################################

# Load Seurat objects
load("./Input/2025-03-17_10x_ovary_Granulosa_combined_Seurat_object.RData")
load("./Input/2025-03-17_10x_ovary_Theca_combined_Seurat_object.RData")

# Load lasso features
Granulosa.features <- read.table("./Input/2025-03-19_Granulosa_Lasso_Important_Features.txt", sep = "\t", header = TRUE)
Theca.features <- read.table("./Input/2025-03-19_Theca_Lasso_Important_Features.txt", sep = "\t", header = TRUE)

colnames(Granulosa.features) <- colnames(Theca.features) <- c("Feature", "Coef")

###############################################################
# 2. Perform ORA - top 500 features (abs coef)
###############################################################

Granulosa.features <- Granulosa.features[order(abs(Granulosa.features$Coef), decreasing = TRUE), ]
Theca.features <- Theca.features[order(abs(Theca.features$Coef), decreasing = TRUE), ]

top_genes.Granulosa <- head(Granulosa.features$Feature, 500)
top_genes.Theca <- head(Theca.features$Feature, 500)

universe_genes.Granulosa <- rownames(seurat.Granulosa.combined)
universe_genes.Theca <- rownames(seurat.Theca.combined)

# GO ALL (FDR < 0.05)
enriched_GO.Granulosa <- enrichGO(gene = top_genes.Granulosa, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", universe = universe_genes.Granulosa, ont = "ALL", pvalueCutoff = 0.05)       # 36 terms
enriched_GO.Theca <- enrichGO(gene = top_genes.Theca, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", universe = universe_genes.Theca, ont = "ALL", pvalueCutoff = 0.05)                   # 11 terms

# Save ORA results table
write.table(enriched_GO.Granulosa@result, file = paste0(Sys.Date(), "_Granulosa_ORA_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(enriched_GO.Theca@result, file = paste0(Sys.Date(), "_Theca_ORA_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Generate plots
pdf(paste0(Sys.Date(), "_Granulosa_ORA_dotplot.pdf"), width = 10, height = 5)
dotplot(enriched_GO.Granulosa, showCategory = 20) + ggtitle("Granulosa - ORA (GO ALL)") +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_gradientn(colors = c("darkorange2", "purple"), limits = c(0, 0.1))
dev.off()

pdf(paste0(Sys.Date(), "_Theca_ORA_dotplot.pdf"), width = 10, height = 5)
dotplot(enriched_GO.Theca, showCategory = 20) + ggtitle("Theca - ORA (GO ALL)") + 
  scale_size_continuous(range = c(3, 10)) +
  scale_color_gradientn(colors = c("darkorange2", "purple"), limits = c(0, 0.1))
dev.off()

###############################################################
# 3. Assess expression of features by age - PB data
###############################################################

# Process Pseudobulk Data for Granulosa & Theca
process_pseudobulk <- function(cell_type, seurat_obj, features) {
  
  meta <- seurat_obj@meta.data
  meta$Model <- ifelse(grepl("weeks", rownames(meta)), "AgingMA", ifelse(grepl("CTL|VCD", rownames(meta)), "VCDmodel", ifelse(grepl("Foxl2", rownames(meta)), "Foxl2model", "Agingmodel")))
  meta$Group <- ifelse(grepl("VCD|het", rownames(meta)), "Test", "Control")
  seurat_obj@meta.data <- meta
  controls <- subset(seurat_obj, subset = Group == "Control")
  
  sce <- as.SingleCellExperiment(controls)
  sce <- prepSCE(sce, kid = "Group", gid = "Model", sid = "Age_num", drop = TRUE)
  pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
  counts.pb <- pb@assays@data
  
  keep_genes <- rowSums(counts.pb$Control >= 5) >= 4
  filtered_counts <- counts.pb$Control[keep_genes, ]
  
  meta.pb <- seurat_obj@meta.data
  meta.pb$Age_num <- as.character(meta.pb$Age_num)
  meta.pb <- meta.pb[meta.pb$Group == "Control" & meta.pb$Age_num %in% colnames(filtered_counts), ]
  
  metadata.pb <- meta.pb[!duplicated(meta.pb$Age_num), c("Age_num", "Batch", "Model")]
  rownames(metadata.pb) <- metadata.pb$Age_num
  metadata.pb$Batch <- factor(metadata.pb$Batch)
  metadata.pb$Age_num <- as.numeric(as.character(metadata.pb$Age_num))
  
  metadata.pb <- metadata.pb[colnames(filtered_counts), ]
  
  dds <- DESeqDataSetFromMatrix(countData = filtered_counts, colData = metadata.pb, design = ~ Batch + Age_num)
  dds <- DESeq(dds)
  vst_data <- getVarianceStabilizedData(dds)
  
  common.genes <- intersect(rownames(filtered_counts), features$Feature)
  feature.vst <- vst_data[common.genes, ]
  metadata.pb <- metadata.pb[colnames(vst_data), ]
  
  vst_corrected <- removeBatchEffect(vst_data, batch = metadata.pb$Batch)
  feature.vst_corrected <- vst_corrected[common.genes, ]
  
  pdf(paste0(Sys.Date(), "_MeMo_", cell_type, "_clock_features_expression_PB_vst_counts_heatmap_batch_treated.pdf"), width = 5, height = 20)
  pheatmap::pheatmap(feature.vst_corrected, cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = FALSE, scale = "row", 
                     colorRampPalette(rev(c("#CC3333", "#FF9999", "#FFCCCC", "white", "#CCCCFF", "#9999FF", "#333399")))(50),
                     main = paste(cell_type, "clock features (n=", nrow(feature.vst_corrected), "genes)"),
                     cellwidth = 20, border = NA, cellheight = 0.5)
  dev.off()  
  
}

process_pseudobulk("Granulosa", seurat.Granulosa.combined, Granulosa.features)
process_pseudobulk("Theca", seurat.Theca.combined, Theca.features)

###############################################################
sink(file = paste(Sys.Date(), "_lasso_feature_analysis_session_Info.txt"))
sessionInfo()
sink()
