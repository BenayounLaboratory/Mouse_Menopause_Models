setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/1_Preprocessing")
options(stringsAsFactors = FALSE)

library('celda')
library("singleCellTK")
library('Seurat')
library(bitops)
library(sctransform)
library(clustree)
library(scales)
library(dplyr)

library(DoubletFinder)

# Add cxds_bcds_hybrid related packages
library(scds)
library(scater)
library(bitops)

# rm(list = ls())

################################################################################
# 10x ovary Benayoun lab Aging - middle ages data
# Pre-processing - DecontX, Singlet filtering
################################################################################

################################################################################
# 1. Detect ambient RNA - DecontX
# https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html
################################################################################

# Read counts from CellRanger output
counts.13weeks <- Read10X("/Volumes/OIProject_II/10XGenomics_data/Benayoun_lab_data_backup/2_Cellranger/AC_middle_ages/AC_13weeks/outs/filtered_feature_bc_matrix/")
counts.76weeks <- Read10X("/Volumes/OIProject_II/10XGenomics_data/Benayoun_lab_data_backup/2_Cellranger/AC_middle_ages/AC_76weeks/outs/filtered_feature_bc_matrix/")
counts.86weeks <- Read10X("/Volumes/OIProject_II/10XGenomics_data/Benayoun_lab_data_backup/2_Cellranger/AC_middle_ages/AC_86weeks/outs/filtered_feature_bc_matrix/")

counts.raw.13weeks <- Read10X("/Volumes/OIProject_II/10XGenomics_data/Benayoun_lab_data_backup/2_Cellranger/AC_middle_ages/AC_13weeks/outs/raw_feature_bc_matrix/")
counts.raw.76weeks <- Read10X("/Volumes/OIProject_II/10XGenomics_data/Benayoun_lab_data_backup/2_Cellranger/AC_middle_ages/AC_76weeks/outs/raw_feature_bc_matrix/")
counts.raw.86weeks <- Read10X("/Volumes/OIProject_II/10XGenomics_data/Benayoun_lab_data_backup/2_Cellranger/AC_middle_ages/AC_86weeks/outs/raw_feature_bc_matrix/")

# Create a SingleCellExperiment object and run decontX
sce.13weeks <- SingleCellExperiment(list(counts = counts.13weeks))
sce.76weeks <- SingleCellExperiment(list(counts = counts.76weeks))
sce.86weeks <- SingleCellExperiment(list(counts = counts.86weeks))

sce.raw.13weeks <- SingleCellExperiment(list(counts = counts.raw.13weeks))
sce.raw.76weeks <- SingleCellExperiment(list(counts = counts.raw.76weeks))
sce.raw.86weeks <- SingleCellExperiment(list(counts = counts.raw.86weeks))

sce.13weeks <- decontX(sce.13weeks, background = sce.raw.13weeks)
sce.76weeks <- decontX(sce.76weeks, background = sce.raw.76weeks)
sce.86weeks <- decontX(sce.86weeks, background = sce.raw.86weeks)

# Save DecontX result
save(sce.13weeks, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_DecontX_SCE_object_sce.13weeks.RData"))
save(sce.76weeks, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_DecontX_SCE_object_sce.76weeks.RData"))
save(sce.86weeks, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_DecontX_SCE_object_sce.86weeks.RData"))

# Plot UMAP
UMAP.13weeks <- reducedDim(sce.13weeks, "decontX_UMAP")
UMAP.76weeks <- reducedDim(sce.76weeks, "decontX_UMAP")
UMAP.86weeks <- reducedDim(sce.86weeks, "decontX_UMAP")

plotDecontXContamination(sce.13weeks)
plotDecontXContamination(sce.76weeks)
plotDecontXContamination(sce.86weeks)

# Plot example box plots of contamination scores
boxplot(sce.13weeks$decontX_contamination)
boxplot(sce.76weeks$decontX_contamination)
boxplot(sce.86weeks$decontX_contamination)

# Create Seurat objects
seurat.13weeks <- CreateSeuratObject(round(decontXcounts(sce.13weeks)), meta.data=as.data.frame(colData(sce.13weeks)))
seurat.76weeks <- CreateSeuratObject(round(decontXcounts(sce.76weeks)), meta.data=as.data.frame(colData(sce.76weeks)))
seurat.86weeks <- CreateSeuratObject(round(decontXcounts(sce.86weeks)), meta.data=as.data.frame(colData(sce.86weeks)))

ovary.AC.MA <- merge(seurat.13weeks, 
                     y =  c(seurat.76weeks,
                            seurat.86weeks), 
                     add.cell.ids = c("13weeks",
                                      "76weeks",
                                      "86weeks"), 
                     project = "10x_ovary_AC_MA")

ovary.AC.MA
# An object of class Seurat 
# 32285 features across 24797 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)

save(ovary.AC.MA, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_middle_ages_post-DecontX_Seurat_object.RData"))

# Remove intermediate files
rm(counts.13weeks, counts.76weeks, counts.86weeks,
   counts.raw.13weeks, counts.raw.76weeks, counts.raw.86weeks,
   sce.13weeks, sce.76weeks, sce.86weeks,
   sce.raw.13weeks, sce.raw.76weeks, sce.raw.86weeks,
   seurat.13weeks, seurat.76weeks, seurat.86weeks)

################################################################################
# 2. Add metadata to Seurat object
################################################################################

# create Library label
Library <- rep("NA", length(colnames(ovary.AC.MA@assays$RNA)))
Library[grep("13weeks", colnames(ovary.AC.MA@assays$RNA))] <- "13weeks"
Library[grep("76weeks", colnames(ovary.AC.MA@assays$RNA))] <- "76weeks"
Library[grep("86weeks", colnames(ovary.AC.MA@assays$RNA))] <- "86weeks"

Library <- data.frame(Library)
rownames(Library) <- colnames(ovary.AC.MA@assays$RNA)

ovary.AC.MA <- AddMetaData(object = ovary.AC.MA, metadata = as.vector(Library), col.name = "Library")

################################################################################
# 3. Seurat QC
# Filter parameters: nFeature_RNA > 500
#                    500 < nCount_RNA < 1e5
#                    percent.mito < 15
#                    decontX_contamination < 0.25
# https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
################################################################################

ovary.AC.MA <- SetIdent(ovary.AC.MA, value = "Library")

########## QC - remove low/null genes ########## 
# https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
min.value = 0
min.cells = 20
genes.use <- rownames(ovary.AC.MA@assays$RNA)
num.cells <- Matrix::rowSums(ovary.AC.MA@assays$RNA@counts > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])

# remove low/null genes
ovary.AC.MA.filt <- subset(ovary.AC.MA, features = genes.use)
ovary.AC.MA.filt
# An object of class Seurat 
# 17092 features across 24797 samples within 1 assay 
# Active assay: RNA (17092 features, 0 variable features)

############### QC - mitochondrial genes ###############

ovary.AC.MA.filt[["percent.mito"]] <- PercentageFeatureSet(ovary.AC.MA.filt, pattern = "^mt-")
head(ovary.AC.MA.filt@meta.data)
#                               orig.ident nCount_RNA nFeature_RNA decontX_contamination decontX_clusters Library percent.mito
# 13weeks_AAACCCAAGCGTTACT-1 SeuratProject      11280         4093            0.38467538                1 13weeks     2.952128
# 13weeks_AAACCCACAGCATGCC-1 SeuratProject       4093         1832            0.20553968                2 13weeks     7.720498
# 13weeks_AAACCCACAGCCTTCT-1 SeuratProject       1877          963            0.17239475                2 13weeks     3.036761
# 13weeks_AAACCCAGTACTCAAC-1 SeuratProject       5658         2430            0.24012192                2 13weeks     1.908802
# 13weeks_AAACCCAGTAGAGATT-1 SeuratProject        249          172            0.89687978                2 13weeks     2.409639
# 13weeks_AAACCCAGTTCTCACC-1 SeuratProject       1987          925            0.02099857                2 13weeks     1.207851

pdf(paste(Sys.Date(),"10x_ovary_Benayoun_lab_AC_middle_ages_violinPlots_QC_gene_UMI_mito_decontX.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.AC.MA.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# plot important QC metrics
plot1 <- FeatureScatter(ovary.AC.MA.filt, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(ovary.AC.MA.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(ovary.AC.MA.filt, feature1 = "nCount_RNA", feature2 = "decontX_contamination")

pdf(paste(Sys.Date(),"10x_ovary_Benayoun_lab_AC_middle_ages_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
CombinePlots(plots = list(plot1, plot2, plot3))
dev.off()

# get clean cells
ovary.AC.MA.filt <- subset(ovary.AC.MA.filt, subset = nFeature_RNA > 500 & nCount_RNA > 500 & nCount_RNA < 1e5 & percent.mito < 15 & decontX_contamination < 0.25)
ovary.AC.MA.filt
# An object of class Seurat 
# 17092 features across 16886 samples within 1 assay 
# Active assay: RNA (17092 features, 0 variable features)

pdf(paste(Sys.Date(),"10x_ovary_Benayoun_lab_AC_middle_ages_violinPlots_QC_gene_UMI_mito_decontX_by_sample_POST_FILTER.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = ovary.AC.MA.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "decontX_contamination"), ncol = 2, pt.size=0)
dev.off()

# Save data
save(ovary.AC.MA.filt, file = paste(Sys.Date(),"10x_ovary_Benayoun_lab_AC_middle_ages_Seurat_object_post_QC.RData",sep = "_"))

# Normalize data & SCTransform
ovary.AC.MA.filt <- NormalizeData(ovary.AC.MA.filt, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.AC.MA.filt <- SCTransform(object = ovary.AC.MA.filt, vars.to.regress = c("nFeature_RNA", "percent.mito", "Library"))

save(ovary.AC.MA.filt, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_middle_ages_Seurat_object_postSCT.RData"))

ovary.AC.MA.filt <- RunPCA(ovary.AC.MA.filt, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_AC_middle_ages_ElbowPlot.pdf"))
ElbowPlot(ovary.AC.MA.filt, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.AC.MA.filt[["pca"]]@stdev / sum(ovary.AC.MA.filt[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 41

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 19

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 19

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_AC_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.AC.MA.filt <- RunUMAP(ovary.AC.MA.filt, dims = 1:pcs)
ovary.AC.MA.filt <- FindNeighbors(ovary.AC.MA.filt, dims = 1:pcs)
ovary.AC.MA.filt <- FindClusters(object = ovary.AC.MA.filt, resolution = 2)

################################################################################
# 4. Doublet detection part 1: DoubletFinder
################################################################################

########## 0. Doublet rate calculation ########## 
##################################################
#### 0. Assumed doublet information/to calculate %age for prediction
# Targeted Cell Recovery  # of Cells Loaded	Barcodes Detected	Singlets Multiplets	Multiplet Rate
#  3,000	           4,950           ~3,000        ~2,900	   ~80	    ~2.4%
#  4,000	           6,600           ~3,900	       ~3,800	   ~140     ~3.2%
#  5,000	           8,250           ~4,800        ~4,600	   ~210     ~4.0%
#  6,000            9,900           ~5,700	       ~5,400	   ~300	    ~4.8%
#  7,000            11,550	         ~6,600	       ~6,200	   ~400	    ~5.6%
#  8,000            13,200	         ~7,500	       ~7,000	   ~510	    ~6.4%
#  9,000            14,850	         ~8,400	       ~7,700	   ~640	    ~7.2%
# 10,000            16,500	         ~9,200	       ~8,400	   ~780	    ~8.0%
# 12,000            19,800	         ~10,900	     ~9,800	   ~1,100   ~9.6%
# 14,000            23,100	         ~12,500	     ~11,000	 ~1,500   ~11.2%
# 16,000            26,400	         ~14,000	     ~12,100	 ~1,900   ~12.8%
# 18,000            29,700	         ~15,500	     ~13,100	 ~2,300   ~14.4%
# 20,000            33,000	         ~16,900	     ~14,100	 ~2,800   ~16.0%

pred.10x.dblt <- data.frame( "cell_number" = c(3000,4000,5000,6000,7000,8000,9000, 10000, 12000, 14000, 16000, 18000, 20000),
                             "dblt_rate"   = c(2.4 ,3.2 ,4.0 ,4.8 ,5.6 ,6.4 ,7.2 , 8.0  , 9.6, 11.2, 12.8, 14.4, 16.0))

pred_dblt_lm <- lm(dblt_rate ~ cell_number, data = pred.10x.dblt)

pdf(paste0(Sys.Date(), "_10x_cell_number_vs_doublet_rate.pdf"))
plot(dblt_rate ~ cell_number, data = pred.10x.dblt)
abline(pred_dblt_lm, col = "red", lty = "dashed")
dev.off()

# Use distinct recovered cell count number to estimate doublets in dataset
#predict(pred_dblt_lm, data.frame("cell_number" = 6000))

########################################################
#### need to split by 10x sample to make sure to identify real doublets
# will run on one object at a time
ovary.AC.MA.filt.list <- SplitObject(ovary.AC.MA.filt, split.by = "Library")

## Assume doublet rate based on 10x information (+ 10% cells to be safe)
pred.dblt.rate <- predict(pred_dblt_lm, data.frame("cell_number" = 1.1 * unlist(lapply(ovary.AC.MA.filt.list, ncol))))/100
pred.dblt.rate
#   13weeks   76weeks   86weeks 
# 0.0489104 0.0513216 0.0483648 

############### Run DoubletFinder ###############
# loop over sample
for (i in 1:length(ovary.AC.MA.filt.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list.ovary.AC.MA <- paramSweep_v3(ovary.AC.MA.filt.list[[i]], PCs = 1:18, sct = TRUE, num.cores = 4)
  sweep.stats.ovary.AC.MA    <- summarizeSweep(sweep.res.list.ovary.AC.MA, GT = FALSE)
  bcmvn.ovary.AC.MA          <- find.pK(sweep.stats.ovary.AC.MA)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yield row number
  pk.ovary.AC.MA <- as.numeric(as.character(bcmvn.ovary.AC.MA[as.numeric(bcmvn.ovary.AC.MA$pK[bcmvn.ovary.AC.MA$BCmetric == max(bcmvn.ovary.AC.MA$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(ovary.AC.MA.filt.list[[i]]@meta.data$seurat_clusters)     ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi       <- round((pred.dblt.rate[[i]]+ 0.025)                                        ## use predicted doublet rate from 10x + 2.5% of possibly non dissociated cells (estimating)
                          *length(ovary.AC.MA.filt.list[[i]]@meta.data$orig.ident))          ## Assuming 10x provided doublet formation rate, based on observed cell yield
  
  ## Run DoubletFinder with varying classification stringencies
  ovary.AC.MA.filt.list[[i]] <- doubletFinder_v3(ovary.AC.MA.filt.list[[i]], PCs = 1:16, pN = 0.25, pK = pk.ovary.AC.MA, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(ovary.AC.MA.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.AC.MA.filt.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(ovary.AC.MA.filt.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(ovary.AC.MA.filt.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(ovary.AC.MA.filt.list)) {
  
  pdf(paste(Sys.Date(),"ovary_AC_MA",names(ovary.AC.MA.filt.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(ovary.AC.MA.filt.list[[i]], reduction = "umap", group.by = "DoubletFinder"))
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
# https://satijalab.org/seurat/articles/merge_vignette.html
ovary.AC.MA.DFsinglets <- merge(ovary.AC.MA.filt.list[[1]],
                                y = c(ovary.AC.MA.filt.list[[2]], ovary.AC.MA.filt.list[[3]]),
                                project = "ovary.AC.MA")
ovary.AC.MA.DFsinglets
# An object of class Seurat 
# 34184 features across 16886 samples within 2 assays 
# Active assay: SCT (17092 features, 0 variable features)
# 1 other assay present: RNA

# remove pANN columns that are 10xGenomics library lane specific
ovary.AC.MA.DFsinglets@meta.data <- ovary.AC.MA.DFsinglets@meta.data[,-grep("pANN",colnames(ovary.AC.MA.DFsinglets@meta.data))]

################################################################################
# 5. Doublet detection part 2: scds
################################################################################

############### Run scds:single cell doublet scoring (hybrid method) ###############
# https://www.bioconductor.org/packages/release/bioc/vignettes/scds/inst/doc/scds.html
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches                            

# create scds working object - convert list to SingleCellExperiment
ovary.AC.MA.filt.list.scds <- lapply(ovary.AC.MA.filt.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(ovary.AC.MA.filt.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  ovary.AC.MA.filt.list.scds[[i]] <- cxds_bcds_hybrid(ovary.AC.MA.filt.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.dblt.rate[i])*ncol(ovary.AC.MA.filt.list.scds[[i]]))                         ## Assume doublets based on nuclei isolation protocol performance
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(ovary.AC.MA.filt.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  ovary.AC.MA.filt.list.scds[[i]]$scds <- "Singlet"
  ovary.AC.MA.filt.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots
for (i in 1:length(ovary.AC.MA.filt.list.scds)) {
  
  p <- plotReducedDim(ovary.AC.MA.filt.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"ovary_AC_MA", names(ovary.AC.MA.filt.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to DoubletFinder annotated Seurat object
ovary.AC.MA.DFsinglets@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(ovary.AC.MA.filt.list.scds)) {
  
  # for each object compare and move doublet annotations over
  ovary.AC.MA.DFsinglets@meta.data[colnames(ovary.AC.MA.filt.list.scds[[i]]), ]$scds_hybrid <- ovary.AC.MA.filt.list.scds[[i]]$scds
  
}

################################################################################
# 6. Filter singlets
################################################################################

table(ovary.AC.MA.DFsinglets@meta.data$DoubletFinder, ovary.AC.MA.DFsinglets@meta.data$scds_hybrid)
#         Doublet Singlet
# Doublet      61    1198
# Singlet     776   14851

ovary.AC.MA.DFsinglets@meta.data$DoubletCall <- ifelse(bitOr(ovary.AC.MA.DFsinglets@meta.data$DoubletFinder == "Doublet", ovary.AC.MA.DFsinglets@meta.data$scds_hybrid == "Doublet") > 0, 
                                                       "Doublet", "Singlet")
table(ovary.AC.MA.DFsinglets@meta.data$DoubletCall)
#  Doublet Singlet 
#     2035   14851 

# re-run dimensionality reduction for plotting purposes
ovary.AC.MA.DFsinglets <- SCTransform(object = ovary.AC.MA.DFsinglets, vars.to.regress = c("nFeature_RNA", "percent.mito", "Library", "decontX_contamination"))
ovary.AC.MA.DFsinglets <- RunPCA(ovary.AC.MA.DFsinglets, npcs = 50)
ovary.AC.MA.DFsinglets <- RunUMAP(ovary.AC.MA.DFsinglets, dims = 1:50)

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_MA_DoubletCall_UMAP.pdf"))
DimPlot(ovary.AC.MA.DFsinglets, reduction = "umap", group.by = "DoubletCall")
dev.off()

# save annotated object
save(ovary.AC.MA.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_MA_Seurat_object_with_AnnotatedDoublets.RData"))

### extract/subset only singlets
# save data for singlets df
ovary.AC.MA.DFsinglets   <- subset(ovary.AC.MA.DFsinglets, subset = DoubletCall %in% "Singlet")  # only keep singlets
ovary.AC.MA.DFsinglets
# An object of class Seurat 
# 41239 features across 9559 samples within 2 assays 
# Active assay: SCT (20617 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

table(ovary.AC.MA.DFsinglets@meta.data$Library)
#  OF1  OF2  YF1  YF2 
# 2554 1125 2944 2936 

# save filtered/annotated object
save(ovary.AC.MA.DFsinglets, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_MA_Seurat_object_SINGLETS.RData"))

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Aging_MA_Seurat_QC_and_filter_processing_session_Info.txt", sep =""))
sessionInfo()
sink()
