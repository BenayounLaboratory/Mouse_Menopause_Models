setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/8_PB_TE")

library(Seurat)
library(SeuratDisk)
library(Matrix)

###############################################################
# Preprocess scTE output for TE analysis
# Aging model dataset
###############################################################

###############################################################
# 1. Import data files and combine data
###############################################################

# Load annotated Seurat object
load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/8_PB_TE")

# Extract RNA counts
ovary.AC.RNA.counts <- ovary.AC@assays$RNA@counts

# Read scTE files
SeuratDisk::Convert("../scTE_out_YF_1.h5ad", "../scTE_out_YF_1.h5seurat")
SeuratDisk::Convert("../scTE_out_YF_2.h5ad", "../scTE_out_YF_2.h5seurat")
SeuratDisk::Convert("../scTE_out_OF_1.h5ad", "../scTE_out_OF_1.h5seurat")
SeuratDisk::Convert("../scTE_out_OF_2.h5ad", "../scTE_out_OF_2.h5seurat")

YF1_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_out_YF_1.h5seurat")
YF2_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_out_YF_2.h5seurat")
OF1_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_out_OF_1.h5seurat")
OF2_TEs.h5 <- SeuratDisk::LoadH5Seurat("../scTE_out_OF_2.h5seurat")

YF1_TEs <- YF1_TEs.h5@assays$RNA@data
YF2_TEs <- YF2_TEs.h5@assays$RNA@data
OF1_TEs <- OF1_TEs.h5@assays$RNA@data
OF2_TEs <- OF2_TEs.h5@assays$RNA@data

dim(YF1_TEs)     # [1] 56422  4919
dim(YF2_TEs)     # [1] 56422  4795
dim(OF1_TEs)     # [1] 56422  4079
dim(OF2_TEs)     # [1] 56422  1899

# Rename colnames to match Seurat object
colnames(YF1_TEs) <- paste0("YF_1_",colnames(YF1_TEs),"-1")
colnames(YF2_TEs) <- paste0("YF_2_",colnames(YF2_TEs),"-1")
colnames(OF1_TEs) <- paste0("OF_1_",colnames(OF1_TEs),"-1")
colnames(OF2_TEs) <- paste0("OF_2_",colnames(OF2_TEs),"-1")

# Select cell barcodes detected by both RNA and scTE
joint.bcs.YF1 <- intersect(colnames(YF1_TEs), colnames(ovary.AC))
joint.bcs.YF2 <- intersect(colnames(YF2_TEs), colnames(ovary.AC))
joint.bcs.OF1 <- intersect(colnames(OF1_TEs), colnames(ovary.AC))
joint.bcs.OF2 <- intersect(colnames(OF2_TEs), colnames(ovary.AC))

# Select overlapping barcodes from scTE for merging
ovary.AC.YF1           <- YF1_TEs[, joint.bcs.YF1]
ovary.AC.YF2           <- YF2_TEs[, joint.bcs.YF2]
ovary.AC.OF1           <- OF1_TEs[, joint.bcs.OF1]
ovary.AC.OF2           <- OF2_TEs[, joint.bcs.OF2]

# Check if row names of all objects are of same order before merging
all_same_order <- identical(rownames(ovary.AC.YF1), rownames(ovary.AC.YF2)) &&
                  identical(rownames(ovary.AC.YF1), rownames(ovary.AC.OF1)) &&
                  identical(rownames(ovary.AC.YF1), rownames(ovary.AC.OF2))
print(all_same_order)     # [1] TRUE

# Combine scTE results
ovary.AC.scTE <- cbind(ovary.AC.YF1,ovary.AC.YF2,ovary.AC.OF1,ovary.AC.OF2)

# Filter unique TEs (exclude feature names found in the original Seurat object)
my.scTE.unique <- setdiff(rownames(ovary.AC.scTE), rownames(ovary.AC.RNA.counts))
length(my.scTE.unique)     # [1] 36177
my.ovary.AC.scTE.unique <- ovary.AC.scTE[my.scTE.unique,]

# Add prefix to TE names
rownames(my.ovary.AC.scTE.unique) <- paste0("TE_",rownames(my.ovary.AC.scTE.unique))

# Combine my.ovary.AC and scTE data
# Check that colnames are in the same order
if (!identical(colnames(ovary.AC.RNA.counts), colnames(my.ovary.AC.scTE.unique))) {
  # Reorder the TE counts to match the RNA counts
  my.ovary.AC.scTE.unique <- my.ovary.AC.scTE.unique[, colnames(ovary.AC.RNA.counts)]
}

all_same_order <- identical(colnames(ovary.AC.RNA.counts), colnames(my.ovary.AC.scTE.unique))
print(all_same_order)     # [1] TRUE

# Combine datasets
ovary.AC.genes.TEs.combined <- rbind(ovary.AC.RNA.counts, my.ovary.AC.scTE.unique)

ovary.AC.genes.TEs.combined.seurat <- CreateSeuratObject(counts = ovary.AC.genes.TEs.combined)

ovary.AC.genes.TEs.combined.seurat
# An object of class Seurat 
# 56804 features across 9388 samples within 1 assay 
# Active assay: RNA (56804 features, 0 variable features)

# Transfer metadata from ovary.AC - Age, Library, celltype.level1, celltype.level2

all_same_order <- identical(colnames(ovary.AC), colnames(ovary.AC.genes.TEs.combined.seurat))
print(all_same_order)     # [1] TRUE

ovary.AC.genes.TEs.combined.seurat$Age <- ovary.AC$Age
ovary.AC.genes.TEs.combined.seurat$Library <- ovary.AC$Library
ovary.AC.genes.TEs.combined.seurat$celltype.level1 <- ovary.AC$celltype.level1
ovary.AC.genes.TEs.combined.seurat$celltype.level2 <- ovary.AC$celltype.level2
ovary.AC.genes.TEs.combined.seurat$decontX_contamination <- ovary.AC$decontX_contamination

ovary.AC.genes.TEs.combined.seurat$Age <- factor(ovary.AC.genes.TEs.combined.seurat$Age, levels = c("YF","OF"))

###############################################################
# 2. Process and QC combined Seurat object
###############################################################

############### QC ###############
ovary.AC.genes.TEs.combined.seurat[["percent.mito"]] <- PercentageFeatureSet(ovary.AC.genes.TEs.combined.seurat, pattern = "^mt-")

# Save data
save(ovary.AC.genes.TEs.combined.seurat, file = paste(Sys.Date(),"Aging_Seurat_object_scTE_combined.RData",sep = "_"))

# Normalize data & SCTransform
ovary.AC.genes.TEs.combined.seurat <- NormalizeData(ovary.AC.genes.TEs.combined.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.AC.genes.TEs.combined.seurat <- SCTransform(object = ovary.AC.genes.TEs.combined.seurat, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))

save(ovary.AC.genes.TEs.combined.seurat, file = paste0(Sys.Date(),"_Aging_Seurat_object_scTE_combined_postSCT.RData"))

ovary.AC.genes.TEs.combined.seurat <- RunPCA(ovary.AC.genes.TEs.combined.seurat, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_Aging_scTE_combined_ElbowPlot.pdf"))
ElbowPlot(ovary.AC.genes.TEs.combined.seurat, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.AC.genes.TEs.combined.seurat[["pca"]]@stdev / sum(ovary.AC.genes.TEs.combined.seurat[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 41

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
pdf(paste0(Sys.Date(), "_Aging_scTE_combined_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.AC.genes.TEs.combined.seurat <- RunUMAP(ovary.AC.genes.TEs.combined.seurat, dims = 1:pcs)
ovary.AC.genes.TEs.combined.seurat <- FindNeighbors(ovary.AC.genes.TEs.combined.seurat, dims = 1:pcs)
ovary.AC.genes.TEs.combined.seurat <- FindClusters(object = ovary.AC.genes.TEs.combined.seurat, resolution = 2)

Idents(ovary.AC.genes.TEs.combined.seurat) <- "celltype.level2"

save(ovary.AC.genes.TEs.combined.seurat, file = paste0(Sys.Date(),"_Aging_Seurat_object_scTE_combined_FINAL.RData"))

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Aging_scTE_combined_Seurat_QC_and_filter_processing_session_Info.txt", sep =""))
sessionInfo()
sink()
