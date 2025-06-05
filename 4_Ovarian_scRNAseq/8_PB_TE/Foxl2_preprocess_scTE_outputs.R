setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/8_PB_TE/Foxl2_preprocess_scTE_outputs.R")

library(Seurat)
library(SeuratDisk)
library(Matrix)
library(ggplot2)

rm(list=ls())

###############################################################
# Preprocess scTE output for TE analysis
# Foxl2 dataset
###############################################################

###############################################################
# 1. Import data files and combine data
###############################################################

# Load annotated Seurat object
load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/8_PB_TE/Input/2024-10-15_10x_ovary_Benayoun_lab_AC_Seurat_object_with_final_annotation.RData")

# Extract RNA counts
ovary.Foxl2.RNA.counts <- ovary.Foxl2$RNA@counts

# Read scTE files
Foxl2_het_midage_1.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_Foxl2_het_midage_1.h5seurat")
Foxl2_het_midage_2.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_Foxl2_het_midage_2.h5seurat")
Foxl2_het_young_1.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_Foxl2_het_young_1.h5seurat")
Foxl2_het_young_2.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_Foxl2_het_young_2.h5seurat")
Foxl2_het.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_Foxl2_het.h5seurat")
Foxl2_wt_midage_1.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_Foxl2_wt_midage_1.h5seurat")
Foxl2_wt_midage_2.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_Foxl2_wt_midage_2.h5seurat")
Foxl2_wt_young_1.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_Foxl2_wt_young_1.h5seurat")
Foxl2_wt_young_2.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_Foxl2_wt_young_2.h5seurat")
Foxl2_wt.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_Foxl2_wt.h5seurat")

Foxl2_het_midage_1_TEs <- Foxl2_het_midage_1.h5@assays$RNA@data
Foxl2_het_midage_2_TEs <- Foxl2_het_midage_2.h5@assays$RNA@data
Foxl2_het_young_1_TEs <- Foxl2_het_young_1.h5@assays$RNA@data
Foxl2_het_young_2_TEs <- Foxl2_het_young_2.h5@assays$RNA@data
Foxl2_het_TEs <- Foxl2_het.h5@assays$RNA@data
Foxl2_wt_midage_1_TEs <- Foxl2_wt_midage_1.h5@assays$RNA@data
Foxl2_wt_midage_2_TEs <- Foxl2_wt_midage_2.h5@assays$RNA@data
Foxl2_wt_young_1_TEs <- Foxl2_wt_young_1.h5@assays$RNA@data
Foxl2_wt_young_2_TEs <- Foxl2_wt_young_2.h5@assays$RNA@data
Foxl2_wt_TEs <- Foxl2_wt.h5@assays$RNA@data

dim(Foxl2_het_midage_1_TEs)   # [1] 56422  1099
dim(Foxl2_het_midage_2_TEs)   # [1] 56422  3166
dim(Foxl2_het_young_1_TEs) # [1] 56422  4792
dim(Foxl2_het_young_2_TEs) # [1] 56422  4270
dim(Foxl2_het_TEs)         # [1] 56422  1798
dim(Foxl2_wt_midage_1_TEs)    # [1] 56422  5193
dim(Foxl2_wt_midage_2_TEs)    # [1] 56422  6491
dim(Foxl2_wt_young_1_TEs)  # [1] 56422  7051
dim(Foxl2_wt_young_2_TEs)  # [1] 56422  5863
dim(Foxl2_wt_TEs)          # [1] 56422  3594

# Rename colnames to match Seurat object
colnames(Foxl2_het_midage_1_TEs) <- paste0("Foxl2_het_midage_1_",colnames(Foxl2_het_midage_1_TEs),"-1")
colnames(Foxl2_het_midage_2_TEs) <- paste0("Foxl2_het_midage_2_",colnames(Foxl2_het_midage_2_TEs),"-1")
colnames(Foxl2_het_young_1_TEs) <- paste0("Foxl2_het_young_1_",colnames(Foxl2_het_young_1_TEs),"-1")
colnames(Foxl2_het_young_2_TEs) <- paste0("Foxl2_het_young_2_",colnames(Foxl2_het_young_2_TEs),"-1")
colnames(Foxl2_het_TEs) <- paste0("Foxl2_het_midage_0_",colnames(Foxl2_het_TEs),"-1")
colnames(Foxl2_wt_midage_1_TEs) <- paste0("Foxl2_wt_midage_1_",colnames(Foxl2_wt_midage_1_TEs),"-1")
colnames(Foxl2_wt_midage_2_TEs) <- paste0("Foxl2_wt_midage_2_",colnames(Foxl2_wt_midage_2_TEs),"-1")
colnames(Foxl2_wt_young_1_TEs) <- paste0("Foxl2_wt_young_1_",colnames(Foxl2_wt_young_1_TEs),"-1")
colnames(Foxl2_wt_young_2_TEs) <- paste0("Foxl2_wt_young_2_",colnames(Foxl2_wt_young_2_TEs),"-1")
colnames(Foxl2_wt_TEs) <- paste0("Foxl2_wt_midage_0_",colnames(Foxl2_wt_TEs),"-1")

# Select cell barcodes detected by both RNA and scTE
joint.bcs.Foxl2_het_midage_1  <- intersect(colnames(Foxl2_het_midage_1_TEs ), colnames(ovary.Foxl2))
joint.bcs.Foxl2_het_midage_2  <- intersect(colnames(Foxl2_het_midage_2_TEs ), colnames(ovary.Foxl2))
joint.bcs.Foxl2_het_young_1  <- intersect(colnames(Foxl2_het_young_1_TEs ), colnames(ovary.Foxl2))
joint.bcs.Foxl2_het_young_2  <- intersect(colnames(Foxl2_het_young_2_TEs ), colnames(ovary.Foxl2))
joint.bcs.Foxl2_het  <- intersect(colnames(Foxl2_het_TEs ), colnames(ovary.Foxl2))
joint.bcs.Foxl2_wt_midage_1  <- intersect(colnames(Foxl2_wt_midage_1_TEs ), colnames(ovary.Foxl2))
joint.bcs.Foxl2_wt_midage_2  <- intersect(colnames(Foxl2_wt_midage_2_TEs ), colnames(ovary.Foxl2))
joint.bcs.Foxl2_wt_young_1  <- intersect(colnames(Foxl2_wt_young_1_TEs ), colnames(ovary.Foxl2))
joint.bcs.Foxl2_wt_young_2  <- intersect(colnames(Foxl2_wt_young_2_TEs ), colnames(ovary.Foxl2))
joint.bcs.Foxl2_wt  <- intersect(colnames(Foxl2_wt_TEs ), colnames(ovary.Foxl2))

# Select overlapping barcodes from scTE for merging
ovary.Foxl2_het_midage_1 <- Foxl2_het_midage_1_TEs [, joint.bcs.Foxl2_het_midage_1 ]
ovary.Foxl2_het_midage_2 <- Foxl2_het_midage_2_TEs [, joint.bcs.Foxl2_het_midage_2 ]
ovary.Foxl2_het_young_1 <- Foxl2_het_young_1_TEs [, joint.bcs.Foxl2_het_young_1 ]
ovary.Foxl2_het_young_2 <- Foxl2_het_young_2_TEs [, joint.bcs.Foxl2_het_young_2 ]
ovary.Foxl2_het <- Foxl2_het_TEs [, joint.bcs.Foxl2_het ]
ovary.Foxl2_wt_midage_1 <- Foxl2_wt_midage_1_TEs [, joint.bcs.Foxl2_wt_midage_1 ]
ovary.Foxl2_wt_midage_2 <- Foxl2_wt_midage_2_TEs [, joint.bcs.Foxl2_wt_midage_2 ]
ovary.Foxl2_wt_young_1 <- Foxl2_wt_young_1_TEs [, joint.bcs.Foxl2_wt_young_1 ]
ovary.Foxl2_wt_young_2 <- Foxl2_wt_young_2_TEs [, joint.bcs.Foxl2_wt_young_2 ]
ovary.Foxl2_wt <- Foxl2_wt_TEs [, joint.bcs.Foxl2_wt ]

# Check if row names of all objects are of same order before merging
all_same_order_1 <- identical(rownames(ovary.Foxl2_het_midage_1), rownames(ovary.Foxl2_het_midage_2)) &&
  identical(rownames(ovary.Foxl2_het_midage_1), rownames(ovary.Foxl2_het_young_1)) &&
  identical(rownames(ovary.Foxl2_het_midage_1), rownames(ovary.Foxl2_het_young_1))
print(all_same_order_1)     # [1] TRUE

all_same_order_2 <- identical(rownames(ovary.Foxl2_het_midage_1), rownames(ovary.Foxl2_het)) &&
  identical(rownames(ovary.Foxl2_het_midage_1), rownames(ovary.Foxl2_wt_midage_1)) &&
  identical(rownames(ovary.Foxl2_het_midage_1), rownames(ovary.Foxl2_wt_midage_2))
print(all_same_order_2)     # [1] TRUE

all_same_order_3 <- identical(rownames(ovary.Foxl2_het_midage_1), rownames(ovary.Foxl2_wt_young_1)) &&
  identical(rownames(ovary.Foxl2_het_midage_1), rownames(ovary.Foxl2_wt_young_2)) &&
  identical(rownames(ovary.Foxl2_het_midage_1), rownames(ovary.Foxl2_wt))
print(all_same_order_3)     # [1] TRUE

# Combine scTE results
ovary.Foxl2.scTE <- cbind(ovary.Foxl2_wt_young_1, ovary.Foxl2_wt_young_2, ovary.Foxl2_het_young_1, ovary.Foxl2_het_young_2,
                          ovary.Foxl2_wt_midage_1, ovary.Foxl2_wt_midage_2, ovary.Foxl2_wt,
                          ovary.Foxl2_het_midage_1, ovary.Foxl2_het_midage_2, ovary.Foxl2_het)

# Filter unique TEs (exclude feature names found in the original Seurat object)
scTE.unique <- setdiff(rownames(ovary.Foxl2.scTE), rownames(ovary.Foxl2.RNA.counts))
length(scTE.unique)     # [1] 34498
ovary.Foxl2.scTE.unique <- ovary.Foxl2.scTE[scTE.unique,]

# Add prefix to TE names
rownames(ovary.Foxl2.scTE.unique) <- paste0("scTE_",rownames(ovary.Foxl2.scTE.unique))

# Combine ovary.Foxl2 and scTE data
# Check that colnames are in the same order
if (!identical(colnames(ovary.Foxl2.RNA.counts), colnames(ovary.Foxl2.scTE.unique))) {
  # Reorder the TE counts to match the RNA counts
  print("Not identical.")
}

# Filter common cell IDs
common.cell.ids <- intersect(colnames(ovary.Foxl2.RNA.counts), colnames(ovary.Foxl2.scTE.unique))
length(common.cell.ids)     # [1] 21301

ovary.Foxl2.scTE.unique <- ovary.Foxl2.scTE.unique[, common.cell.ids]
ovary.Foxl2.RNA.counts <- ovary.Foxl2.RNA.counts[, common.cell.ids]

all_same_order <- identical(colnames(ovary.Foxl2.RNA.counts), colnames(ovary.Foxl2.scTE.unique))
print(all_same_order)     # [1] TRUE

# Combine datasets
ovary.Foxl2.genes.TEs.combined <- rbind(ovary.Foxl2.RNA.counts, ovary.Foxl2.scTE.unique)

common.cell.ids <- intersect(colnames(ovary.Foxl2), colnames(ovary.Foxl2.genes.TEs.combined))

ovary.Foxl2.genes.TEs.combined <- ovary.Foxl2.genes.TEs.combined[, common.cell.ids]

ovary.Foxl2.genes.TEs.combined.seurat <- CreateSeuratObject(counts = ovary.Foxl2.genes.TEs.combined)

ovary.Foxl2.genes.TEs.combined.seurat
# An object of class Seurat 
# 56877 features across 21301 samples within 1 assay 
# Active assay: RNA (56877 features, 0 variable features)

# Transfer metadata from my.ovary.Foxl2 - Age, Library, celltype.level1, celltype.level2

ovary.Foxl2.metadata <- ovary.Foxl2@meta.data
ovary.Foxl2.metadata <- ovary.Foxl2.metadata[common.cell.ids,]

all_same_order <- identical(rownames(ovary.Foxl2.metadata), colnames(ovary.Foxl2.genes.TEs.combined.seurat))
print(all_same_order)     # [1] TRUE

ovary.Foxl2.genes.TEs.combined.seurat$Age <- ovary.Foxl2.metadata$Age
ovary.Foxl2.genes.TEs.combined.seurat$Treatment <- ovary.Foxl2.metadata$Genotype
ovary.Foxl2.genes.TEs.combined.seurat$Library <- ovary.Foxl2.metadata$Library
ovary.Foxl2.genes.TEs.combined.seurat$Batch <- ovary.Foxl2.metadata$Batch
ovary.Foxl2.genes.TEs.combined.seurat$celltype.level1 <- ovary.Foxl2.metadata$celltype.level1
ovary.Foxl2.genes.TEs.combined.seurat$celltype.level2 <- ovary.Foxl2.metadata$celltype.level2
ovary.Foxl2.genes.TEs.combined.seurat$decontX_contamination <- ovary.Foxl2.metadata$decontX_contamination

###############################################################
# 2. Process and QC combined Seurat object
###############################################################

############### QC - mitochondrial genes & TE proportions ###############
ovary.Foxl2.genes.TEs.combined.seurat[["percent.mito"]] <- PercentageFeatureSet(ovary.Foxl2.genes.TEs.combined.seurat, pattern = "^mt-")

# Save data
save(ovary.Foxl2.genes.TEs.combined.seurat, file = paste(Sys.Date(),"Foxl2_Seurat_object_scTE_combined.RData",sep = "_"))

# Normalize data & SCTransform
ovary.Foxl2.genes.TEs.combined.seurat <- NormalizeData(ovary.Foxl2.genes.TEs.combined.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.Foxl2.genes.TEs.combined.seurat <- SCTransform(object = ovary.Foxl2.genes.TEs.combined.seurat, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library", "Batch"))

save(ovary.Foxl2.genes.TEs.combined.seurat, file = paste0(Sys.Date(),"_Foxl2_Seurat_object_scTE_combined_postSCT.RData"))

ovary.Foxl2.genes.TEs.combined.seurat <- RunPCA(ovary.Foxl2.genes.TEs.combined.seurat, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_Foxl2_scTE_combined_ElbowPlot.pdf"))
ElbowPlot(ovary.Foxl2.genes.TEs.combined.seurat, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.Foxl2.genes.TEs.combined.seurat[["pca"]]@stdev / sum(ovary.Foxl2.genes.TEs.combined.seurat[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 41

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 14

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 14

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_Foxl2_scTE_combined_elbowplot_threshmidage_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.Foxl2.genes.TEs.combined.seurat <- RunUMAP(ovary.Foxl2.genes.TEs.combined.seurat, dims = 1:pcs)
ovary.Foxl2.genes.TEs.combined.seurat <- FindNeighbors(ovary.Foxl2.genes.TEs.combined.seurat, dims = 1:pcs)
ovary.Foxl2.genes.TEs.combined.seurat <- FindClusters(object = ovary.Foxl2.genes.TEs.combined.seurat, resolution = 2)

save(ovary.Foxl2.genes.TEs.combined.seurat, file = paste0(Sys.Date(),"_Foxl2_Seurat_object_scTE_combined_FINAL.RData"))

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Foxl2_scTE_combined_Seurat_QC_and_filter_processing_session_Info.txt", sep =""))
sessionInfo()
sink()
