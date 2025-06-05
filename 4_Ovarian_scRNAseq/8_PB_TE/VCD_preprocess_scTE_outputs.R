setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/8_PB_TE")

library(Seurat)
library(SeuratDisk)
library(Matrix)
library(ggplot2)

###############################################################
# Preprocess scTE output for TE analysis
# VCD model dataset
###############################################################

###############################################################
# 1. Import data files and combine data
###############################################################

# Load annotated Seurat object
load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/8_PB_TE/Input/2024-10-24_10x_ovary_Benayoun_lab_VCD_Seurat_object_with_final_annotation.RData")

# Extract RNA counts
ovary.VCD.RNA.counts <- my.ovary.VCD@assays$RNA@counts

# Read scTE files

CTL_3m_30d_1_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_CTL_3m_30d_1.h5seurat")
CTL_3m_90d_1_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_CTL_3m_90d_1.h5seurat")
VCD_3m_30d_1_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_VCD_3m_30d_1.h5seurat")
VCD_3m_90d_1_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_VCD_3m_90d_1.h5seurat")
CTL_10m_30d_1_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_CTL_10m_30d_1.h5seurat")
CTL_10m_90d_1_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_CTL_10m_90d_1.h5seurat")
VCD_10m_30d_1_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_VCD_10m_30d_1.h5seurat")
VCD_10m_90d_1_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_VCD_10m_90d_1.h5seurat")
CTL_3m_30d_2_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_CTL_3m_30d_2.h5seurat")
CTL_3m_90d_2_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_CTL_3m_90d_2.h5seurat")
VCD_3m_30d_2_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_VCD_3m_30d_2.h5seurat")
VCD_3m_90d_2_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_VCD_3m_90d_2.h5seurat")
CTL_10m_30d_2_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_CTL_10m_30d_2.h5seurat")
CTL_10m_90d_2_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_CTL_10m_90d_2.h5seurat")
VCD_10m_30d_2_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_VCD_10m_30d_2.h5seurat")
VCD_10m_90d_2_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_output/scTE_out_VCD_10m_90d_2.h5seurat")

CTL_3m_30d_1_TEs <- CTL_3m_30d_1_TEs.h5@assays$RNA@data
CTL_3m_90d_1_TEs <- CTL_3m_90d_1_TEs.h5@assays$RNA@data
CTL_10m_30d_1_TEs <- CTL_10m_30d_1_TEs.h5@assays$RNA@data
CTL_10m_90d_1_TEs <- CTL_10m_90d_1_TEs.h5@assays$RNA@data
VCD_3m_30d_1_TEs <- VCD_3m_30d_1_TEs.h5@assays$RNA@data
VCD_3m_90d_1_TEs <- VCD_3m_90d_1_TEs.h5@assays$RNA@data
VCD_10m_30d_1_TEs <- VCD_10m_30d_1_TEs.h5@assays$RNA@data
VCD_10m_90d_1_TEs <- VCD_10m_90d_1_TEs.h5@assays$RNA@data
CTL_3m_30d_2_TEs <- CTL_3m_30d_2_TEs.h5@assays$RNA@data
CTL_3m_90d_2_TEs <- CTL_3m_90d_2_TEs.h5@assays$RNA@data
CTL_10m_30d_2_TEs <- CTL_10m_30d_2_TEs.h5@assays$RNA@data
CTL_10m_90d_2_TEs <- CTL_10m_90d_2_TEs.h5@assays$RNA@data
VCD_3m_30d_2_TEs <- VCD_3m_30d_2_TEs.h5@assays$RNA@data
VCD_3m_90d_2_TEs <- VCD_3m_90d_2_TEs.h5@assays$RNA@data
VCD_10m_30d_2_TEs <- VCD_10m_30d_2_TEs.h5@assays$RNA@data
VCD_10m_90d_2_TEs <- VCD_10m_90d_2_TEs.h5@assays$RNA@data

dim(CTL_3m_30d_1_TEs)   # [1] 56422  2646
dim(CTL_3m_90d_1_TEs)   # [1] 56422  3285
dim(CTL_10m_30d_1_TEs)   # [1] 56422  2982
dim(CTL_10m_90d_1_TEs)   # [1] 56422  1804
dim(VCD_3m_30d_1_TEs)   # [1] 56422  2386
dim(VCD_3m_90d_1_TEs)   # [1] 56422  2013
dim(VCD_10m_30d_1_TEs)   # [1] 56422  2371
dim(VCD_10m_90d_1_TEs)   # [1] 56422  2343
dim(CTL_3m_30d_2_TEs)   # [1] 56422  8075
dim(CTL_3m_90d_2_TEs)   # [1] 56422  7602
dim(CTL_10m_30d_2_TEs)   # [1] 56422  9474
dim(CTL_10m_90d_2_TEs)   # [1] 56422  9286
dim(VCD_3m_30d_2_TEs)   # [1] 56422 10000
dim(VCD_3m_90d_2_TEs)   # [1] 56422 10000
dim(VCD_10m_30d_2_TEs)   # [1] 56422 10000
dim(VCD_10m_90d_2_TEs)   # [1] 56422 10000

# Rename colnames to match Seurat object
colnames(CTL_3m_30d_1_TEs) <- paste0("CTL_3m_30d_1_",colnames(CTL_3m_30d_1_TEs),"-1")
colnames(CTL_3m_90d_1_TEs) <- paste0("CTL_3m_90d_1_",colnames(CTL_3m_90d_1_TEs),"-1")
colnames(CTL_10m_30d_1_TEs) <- paste0("CTL_10m_30d_1_",colnames(CTL_10m_30d_1_TEs),"-1")
colnames(CTL_10m_90d_1_TEs) <- paste0("CTL_10m_90d_1_",colnames(CTL_10m_90d_1_TEs),"-1")
colnames(VCD_3m_30d_1_TEs) <- paste0("VCD_3m_30d_1_",colnames(VCD_3m_30d_1_TEs),"-1")
colnames(VCD_3m_90d_1_TEs) <- paste0("VCD_3m_90d_1_",colnames(VCD_3m_90d_1_TEs),"-1")
colnames(VCD_10m_30d_1_TEs) <- paste0("VCD_10m_30d_1_",colnames(VCD_10m_30d_1_TEs),"-1")
colnames(VCD_10m_90d_1_TEs) <- paste0("VCD_10m_90d_1_",colnames(VCD_10m_90d_1_TEs),"-1")
colnames(CTL_3m_30d_2_TEs) <- paste0("CTL_3m_30d_2_",colnames(CTL_3m_30d_2_TEs),"-1")
colnames(CTL_3m_90d_2_TEs) <- paste0("CTL_3m_90d_2_",colnames(CTL_3m_90d_2_TEs),"-1")
colnames(CTL_10m_30d_2_TEs) <- paste0("CTL_10m_30d_2_",colnames(CTL_10m_30d_2_TEs),"-1")
colnames(CTL_10m_90d_2_TEs) <- paste0("CTL_10m_90d_2_",colnames(CTL_10m_90d_2_TEs),"-1")
colnames(VCD_3m_30d_2_TEs) <- paste0("VCD_3m_30d_2_",colnames(VCD_3m_30d_2_TEs),"-1")
colnames(VCD_3m_90d_2_TEs) <- paste0("VCD_3m_90d_2_",colnames(VCD_3m_90d_2_TEs),"-1")
colnames(VCD_10m_30d_2_TEs) <- paste0("VCD_10m_30d_2_",colnames(VCD_10m_30d_2_TEs),"-1")
colnames(VCD_10m_90d_2_TEs) <- paste0("VCD_10m_90d_2_",colnames(VCD_10m_90d_2_TEs),"-1")

# Select cell barcodes detected by both RNA and scTE
joint.bcs.CTL_3m_30d_1 <- intersect(colnames(CTL_3m_30d_1_TEs), colnames(my.ovary.VCD))
joint.bcs.CTL_3m_90d_1 <- intersect(colnames(CTL_3m_90d_1_TEs), colnames(my.ovary.VCD))
joint.bcs.CTL_10m_30d_1 <- intersect(colnames(CTL_10m_30d_1_TEs), colnames(my.ovary.VCD))
joint.bcs.CTL_10m_90d_1 <- intersect(colnames(CTL_10m_90d_1_TEs), colnames(my.ovary.VCD))
joint.bcs.VCD_3m_30d_1 <- intersect(colnames(VCD_3m_30d_1_TEs), colnames(my.ovary.VCD))
joint.bcs.VCD_3m_90d_1 <- intersect(colnames(VCD_3m_90d_1_TEs), colnames(my.ovary.VCD))
joint.bcs.VCD_10m_30d_1 <- intersect(colnames(VCD_10m_30d_1_TEs), colnames(my.ovary.VCD))
joint.bcs.VCD_10m_90d_1 <- intersect(colnames(VCD_10m_90d_1_TEs), colnames(my.ovary.VCD))
joint.bcs.CTL_3m_30d_2 <- intersect(colnames(CTL_3m_30d_2_TEs), colnames(my.ovary.VCD))
joint.bcs.CTL_3m_90d_2 <- intersect(colnames(CTL_3m_90d_2_TEs), colnames(my.ovary.VCD))
joint.bcs.CTL_10m_30d_2 <- intersect(colnames(CTL_10m_30d_2_TEs), colnames(my.ovary.VCD))
joint.bcs.CTL_10m_90d_2 <- intersect(colnames(CTL_10m_90d_2_TEs), colnames(my.ovary.VCD))
joint.bcs.VCD_3m_30d_2 <- intersect(colnames(VCD_3m_30d_2_TEs), colnames(my.ovary.VCD))
joint.bcs.VCD_3m_90d_2 <- intersect(colnames(VCD_3m_90d_2_TEs), colnames(my.ovary.VCD))
joint.bcs.VCD_10m_30d_2 <- intersect(colnames(VCD_10m_30d_2_TEs), colnames(my.ovary.VCD))
joint.bcs.VCD_10m_90d_2 <- intersect(colnames(VCD_10m_90d_2_TEs), colnames(my.ovary.VCD))

# Select overlapping barcodes from scTE for merging
my.ovary.VCD.CTL_3m_30d_1           <- CTL_3m_30d_1_TEs[, joint.bcs.CTL_3m_30d_1]
my.ovary.VCD.CTL_3m_90d_1           <- CTL_3m_90d_1_TEs[, joint.bcs.CTL_3m_90d_1]
my.ovary.VCD.CTL_10m_30d_1           <- CTL_10m_30d_1_TEs[, joint.bcs.CTL_10m_30d_1]
my.ovary.VCD.CTL_10m_90d_1           <- CTL_10m_90d_1_TEs[, joint.bcs.CTL_10m_90d_1]
my.ovary.VCD.VCD_3m_30d_1           <- VCD_3m_30d_1_TEs[, joint.bcs.VCD_3m_30d_1]
my.ovary.VCD.VCD_3m_90d_1           <- VCD_3m_90d_1_TEs[, joint.bcs.VCD_3m_90d_1]
my.ovary.VCD.VCD_10m_30d_1           <- VCD_10m_30d_1_TEs[, joint.bcs.VCD_10m_30d_1]
my.ovary.VCD.VCD_10m_90d_1           <- VCD_10m_90d_1_TEs[, joint.bcs.VCD_10m_90d_1]
my.ovary.VCD.CTL_3m_30d_2           <- CTL_3m_30d_2_TEs[, joint.bcs.CTL_3m_30d_2]
my.ovary.VCD.CTL_3m_90d_2           <- CTL_3m_90d_2_TEs[, joint.bcs.CTL_3m_90d_2]
my.ovary.VCD.CTL_10m_30d_2           <- CTL_10m_30d_2_TEs[, joint.bcs.CTL_10m_30d_2]
my.ovary.VCD.CTL_10m_90d_2           <- CTL_10m_90d_2_TEs[, joint.bcs.CTL_10m_90d_2]
my.ovary.VCD.VCD_3m_30d_2           <- VCD_3m_30d_2_TEs[, joint.bcs.VCD_3m_30d_2]
my.ovary.VCD.VCD_3m_90d_2           <- VCD_3m_90d_2_TEs[, joint.bcs.VCD_3m_90d_2]
my.ovary.VCD.VCD_10m_30d_2           <- VCD_10m_30d_2_TEs[, joint.bcs.VCD_10m_30d_2]
my.ovary.VCD.VCD_10m_90d_2           <- VCD_10m_90d_2_TEs[, joint.bcs.VCD_10m_90d_2]

# Check if row names of all objects are of same order before merging
all_same_order_1 <- identical(rownames(my.ovary.VCD.CTL_3m_30d_1), rownames(my.ovary.VCD.CTL_3m_90d_1)) &&
  identical(rownames(my.ovary.VCD.CTL_3m_30d_1), rownames(my.ovary.VCD.CTL_10m_30d_1)) &&
  identical(rownames(my.ovary.VCD.CTL_3m_30d_1), rownames(my.ovary.VCD.CTL_10m_90d_1))
print(all_same_order_1)     # [1] TRUE

all_same_order_2 <- identical(rownames(my.ovary.VCD.VCD_3m_30d_1), rownames(my.ovary.VCD.VCD_3m_90d_1)) &&
  identical(rownames(my.ovary.VCD.VCD_3m_30d_1), rownames(my.ovary.VCD.VCD_10m_30d_1)) &&
  identical(rownames(my.ovary.VCD.VCD_3m_30d_1), rownames(my.ovary.VCD.VCD_10m_90d_1))
print(all_same_order_2)     # [1] TRUE

all_same_order_3 <- identical(rownames(my.ovary.VCD.CTL_3m_30d_2), rownames(my.ovary.VCD.CTL_3m_90d_2)) &&
  identical(rownames(my.ovary.VCD.CTL_3m_30d_2), rownames(my.ovary.VCD.CTL_10m_30d_2)) &&
  identical(rownames(my.ovary.VCD.CTL_3m_30d_2), rownames(my.ovary.VCD.CTL_10m_90d_2))
print(all_same_order_3)     # [1] TRUE

all_same_order_4 <- identical(rownames(my.ovary.VCD.VCD_3m_30d_2), rownames(my.ovary.VCD.VCD_3m_90d_2)) &&
  identical(rownames(my.ovary.VCD.VCD_3m_30d_2), rownames(my.ovary.VCD.VCD_10m_30d_2)) &&
  identical(rownames(my.ovary.VCD.VCD_3m_30d_2), rownames(my.ovary.VCD.VCD_10m_90d_2))
print(all_same_order_4)     # [1] TRUE

all_same_order_5 <- identical(rownames(my.ovary.VCD.CTL_3m_30d_1), rownames(my.ovary.VCD.VCD_3m_30d_1)) &&
  identical(rownames(my.ovary.VCD.CTL_3m_30d_1), rownames(my.ovary.VCD.CTL_3m_30d_2)) &&
  identical(rownames(my.ovary.VCD.CTL_3m_30d_1), rownames(my.ovary.VCD.VCD_3m_30d_2))
print(all_same_order_5)     # [1] TRUE

# Combine scTE results
my.ovary.VCD.scTE <- cbind(my.ovary.VCD.CTL_3m_30d_1, my.ovary.VCD.CTL_3m_90d_1, my.ovary.VCD.CTL_10m_30d_1, my.ovary.VCD.CTL_10m_90d_1,
                           my.ovary.VCD.VCD_3m_30d_1, my.ovary.VCD.VCD_3m_90d_1, my.ovary.VCD.VCD_10m_30d_1, my.ovary.VCD.VCD_10m_90d_1,
                           my.ovary.VCD.CTL_3m_30d_2, my.ovary.VCD.CTL_3m_90d_2, my.ovary.VCD.CTL_10m_30d_2, my.ovary.VCD.CTL_10m_90d_2,
                           my.ovary.VCD.VCD_3m_30d_2, my.ovary.VCD.VCD_3m_90d_2, my.ovary.VCD.VCD_10m_30d_2, my.ovary.VCD.VCD_10m_90d_2)

# Filter unique TEs (exclude feature names found in the original Seurat object)
my.scTE.unique <- setdiff(rownames(my.ovary.VCD.scTE), rownames(ovary.VCD.RNA.counts))
length(my.scTE.unique)     # [1] 32484
my.ovary.VCD.scTE.unique <- my.ovary.VCD.scTE[my.scTE.unique,]

# Add prefix to TE names
rownames(my.ovary.VCD.scTE.unique) <- paste0("TE_",rownames(my.ovary.VCD.scTE.unique))

# Combine my.ovary.VCD and scTE data
# Check that colnames are in the same order
if (!identical(colnames(ovary.VCD.RNA.counts), colnames(my.ovary.VCD.scTE.unique))) {
  # Reorder the TE counts to match the RNA counts
  print("Not identical.\n")
}
# Filter common cell IDs
common.cell.ids <- intersect(colnames(ovary.VCD.RNA.counts), colnames(my.ovary.VCD.scTE.unique))
length(common.cell.ids)     # [1] 51424

my.ovary.VCD.scTE.unique <- my.ovary.VCD.scTE.unique[, common.cell.ids]
ovary.VCD.RNA.counts <- ovary.VCD.RNA.counts[, common.cell.ids]

all_same_order <- identical(colnames(ovary.VCD.RNA.counts), colnames(my.ovary.VCD.scTE.unique))
print(all_same_order)     # [1] TRUE

# Combine datasets
my.ovary.VCD.genes.TEs.combined <- rbind(ovary.VCD.RNA.counts, my.ovary.VCD.scTE.unique)

common.cell.ids <- intersect(colnames(my.ovary.VCD), colnames(my.ovary.VCD.genes.TEs.combined))

my.ovary.VCD.genes.TEs.combined <- my.ovary.VCD.genes.TEs.combined[, common.cell.ids]

my.ovary.VCD.genes.TEs.combined.seurat <- CreateSeuratObject(counts = my.ovary.VCD.genes.TEs.combined)

my.ovary.VCD.genes.TEs.combined.seurat
# An object of class Seurat 
# 56973 features across 51424 samples within 1 assay 
# Active assay: RNA (56973 features, 0 variable features)

# Transfer metadata from my.ovary.VCD - Age, Library, celltype.level1, celltype.level2

ovary.VCD.metadata <- my.ovary.VCD@meta.data
ovary.VCD.metadata <- ovary.VCD.metadata[common.cell.ids,]

all_same_order <- identical(rownames(ovary.VCD.metadata), colnames(my.ovary.VCD.genes.TEs.combined.seurat))
print(all_same_order)     # [1] TRUE

my.ovary.VCD.genes.TEs.combined.seurat$Age <- ovary.VCD.metadata$Age
my.ovary.VCD.genes.TEs.combined.seurat$Treatment <- ovary.VCD.metadata$Treatment
my.ovary.VCD.genes.TEs.combined.seurat$Duration <- ovary.VCD.metadata$Duration
my.ovary.VCD.genes.TEs.combined.seurat$Library <- ovary.VCD.metadata$Library
my.ovary.VCD.genes.TEs.combined.seurat$Group <- ovary.VCD.metadata$Group
my.ovary.VCD.genes.TEs.combined.seurat$Batch <- ovary.VCD.metadata$Batch
my.ovary.VCD.genes.TEs.combined.seurat$celltype.level1 <- ovary.VCD.metadata$celltype.level1
my.ovary.VCD.genes.TEs.combined.seurat$celltype.level2 <- ovary.VCD.metadata$celltype.level2
my.ovary.VCD.genes.TEs.combined.seurat$decontX_contamination <- ovary.VCD.metadata$decontX_contamination

###############################################################
# 2. Process and QC combined Seurat object
###############################################################

############### QC - mitochondrial genes & TE proportions ###############
my.ovary.VCD.genes.TEs.combined.seurat[["percent.mito"]] <- PercentageFeatureSet(my.ovary.VCD.genes.TEs.combined.seurat, pattern = "^mt-")

# Save data
save(my.ovary.VCD.genes.TEs.combined.seurat, file = paste(Sys.Date(),"VCD_Seurat_object_scTE_combined.RData",sep = "_"))

# Normalize data & SCTransform
my.ovary.VCD.genes.TEs.combined.seurat <- NormalizeData(my.ovary.VCD.genes.TEs.combined.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
my.ovary.VCD.genes.TEs.combined.seurat <- SCTransform(object = my.ovary.VCD.genes.TEs.combined.seurat, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library", "Batch"))

save(my.ovary.VCD.genes.TEs.combined.seurat, file = paste0(Sys.Date(),"_VCD_Seurat_object_scTE_combined_postSCT.RData"))

my.ovary.VCD.genes.TEs.combined.seurat <- RunPCA(my.ovary.VCD.genes.TEs.combined.seurat, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_VCD_scTE_combined_ElbowPlot.pdf"))
ElbowPlot(my.ovary.VCD.genes.TEs.combined.seurat, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- my.ovary.VCD.genes.TEs.combined.seurat[["pca"]]@stdev / sum(my.ovary.VCD.genes.TEs.combined.seurat[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 42

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 18

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 18

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_VCD_scTE_combined_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
my.ovary.VCD.genes.TEs.combined.seurat <- RunUMAP(my.ovary.VCD.genes.TEs.combined.seurat, dims = 1:pcs)
my.ovary.VCD.genes.TEs.combined.seurat <- FindNeighbors(my.ovary.VCD.genes.TEs.combined.seurat, dims = 1:pcs)
my.ovary.VCD.genes.TEs.combined.seurat <- FindClusters(object = my.ovary.VCD.genes.TEs.combined.seurat, resolution = 2)

Idents(my.ovary.VCD.genes.TEs.combined.seurat) <- "celltype.level2"

save(my.ovary.VCD.genes.TEs.combined.seurat, file = paste0(Sys.Date(),"_VCD_Seurat_object_scTE_combined_FINAL.RData"))

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_VCD_scTE_combined_Seurat_QC_and_filter_processing_session_Info.txt", sep =""))
sessionInfo()
sink()
