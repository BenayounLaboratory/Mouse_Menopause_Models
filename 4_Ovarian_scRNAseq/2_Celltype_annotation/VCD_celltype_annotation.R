setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/2_Celltype_annotation")
options(stringsAsFactors = F)

library(Seurat)
library(SeuratWrappers)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scGate)
library(scProportionTest)
library(scran)
library(SingleR)
library(scSorter)
library(tidyverse)

# rm(list = ls())

################################################################################
# 10x ovary celltype annotation:
#                                1. Ptprc positive vs. Ptprc negative cells
#                                2. Ptprc positive sub-celltype
#                                3. Ptprc negative sub-celltype
# VCD model data
################################################################################

################################################################################
# 1. Load data
################################################################################

load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/2_Celltype_annotation/Input/2024-10-23_10x_ovary_Benayoun_lab_VCD_integrated_harmony_post_clustering.RData")
ovary.VCD.harmony
# An object of class Seurat 
# 48841 features across 47386 samples within 2 assays 
# Active assay: SCT (24356 features, 2000 variable features)
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, harmony, umap

ovary.VCD.harmony <- FindNeighbors(object = ovary.VCD.harmony, reduction = "harmony")
ovary.VCD.harmony <- FindClusters(ovary.VCD.harmony, resolution = 4)

annotation <- table(ovary.VCD.harmony@meta.data$seurat_clusters, ovary.VCD.harmony@meta.data$Batch)
# tail(annotation, 20)

################################################################################
# 2. Gate CD45+ cells (Ptprc expressing cells)
################################################################################

#Manually define a simple scGate gating model to purify immune cells using a positive marker Ptprc
my.scGate.model.immune <- gating_model(name = "immune", signature = c("Ptprc"))

# scGate it!
# Gate individual batches separately: https://github.com/carmonalab/scGate/issues/12
ovary.VCD.Ptprc.pos <- scGate(data = ovary.VCD.harmony, model = my.scGate.model.immune, assay = "RNA", verbose = TRUE)
# Running scGate on level 1...
# scGate: Detected 41333 non-pure cells at level 1
### Detected a total of 41333 non-pure cells for selected signatures - 87.23% cells marked for removal (active.ident)

t.pure <- as.matrix(table(ovary.VCD.Ptprc.pos$is.pure, ovary.VCD.Ptprc.pos$seurat_clusters))

# Calculate the total cells in each cluster
total <- t.pure[1,] + t.pure[2,]

# Get cluster numbers where Pure cells are more than 90% of total cells in each cluster
high_purity_clusters <- which(t.pure[1,] / total > 0.9)
high_purity_cluster_numbers <- as.numeric(colnames(t.pure)[high_purity_clusters])

# Print cluster numbers
print(high_purity_cluster_numbers)
# [1] 19 22 29 31 36 40 49 52 53 58 68

# Assign cell type to each seurat cluster: Ptprc (CD45) positive cells clusters: 19 22 29 31 36 40 49 52 53 58 68

ovary.VCD.Ptprc.pos <- SetIdent(ovary.VCD.Ptprc.pos, value = "seurat_clusters")
ovary.VCD <- ovary.VCD.Ptprc.pos

ovary.VCD@meta.data$celltype.level1 <- rep("Ptprc.neg",dim(ovary.VCD@meta.data)[1])
ovary.VCD@meta.data$celltype.level1[ovary.VCD@meta.data$seurat_clusters %in% high_purity_cluster_numbers] <- "Ptprc.pos"

save(ovary.VCD, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_VCD_Ptprc_gated_Seurat_object.RData"))

################################################################################
# 2.Ptprc positive sub-celltype annotation: 
#                                           1. Manual marker annotation
#                                           2. SingleR / scType / scSorter
################################################################################

################################
# Subset Ptprc positive cells and re-cluster
################################

ovary.VCD.Ptprc.pos <- subset(ovary.VCD, subset = celltype.level1 == "Ptprc.pos")

ovary.VCD.Ptprc.pos
# An object of class Seurat 
# 48841 features across 5931 samples within 2 assays 
# Active assay: RNA (24485 features, 0 variable features)
# 1 other assay present: SCT
# 3 dimensional reductions calculated: pca, harmony, umap

DefaultAssay(ovary.VCD.Ptprc.pos) <- "RNA"

ovary.VCD.Ptprc.pos <- DietSeurat(ovary.VCD.Ptprc.pos,
                                  counts = TRUE, 
                                  data = TRUE,
                                  scale.data = FALSE,
                                  features = NULL,
                                  assays = NULL,
                                  dimreducs = NULL,
                                  graphs = NULL)

ovary.VCD.Ptprc.pos[['nCount_SCT']] <- NULL
ovary.VCD.Ptprc.pos[['nFeature_SCT']] <- NULL
ovary.VCD.Ptprc.pos[['SCT_snn_res.2']] <- NULL
ovary.VCD.Ptprc.pos[['seurat_clusters']] <- NULL
ovary.VCD.Ptprc.pos[['DoubletFinder']] <- NULL
ovary.VCD.Ptprc.pos[['scds_hybrid']] <- NULL

########## Normalize, SCT, dimreduc data ##########

ovary.VCD.Ptprc.pos <- NormalizeData(ovary.VCD.Ptprc.pos, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.VCD.Ptprc.pos <- SCTransform(object = ovary.VCD.Ptprc.pos, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))

save(ovary.VCD.Ptprc.pos, file = paste0(Sys.Date(),"_10x_ovary_VCD_Ptprc_pos_cells_Seurat_object_clean_postSCT.RData"))

ovary.VCD.Ptprc.pos <- RunPCA(ovary.VCD.Ptprc.pos, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_VCD_Ptprc_pos_cells_ElbowPlot.pdf"))
ElbowPlot(ovary.VCD.Ptprc.pos, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.VCD.Ptprc.pos[["pca"]]@stdev / sum(ovary.VCD.Ptprc.pos[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 41

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 11

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 11

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_10x_ovary_VCD_Ptprc_pos_cells_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.VCD.Ptprc.pos <- RunUMAP(ovary.VCD.Ptprc.pos, dims = 1:pcs)
ovary.VCD.Ptprc.pos <- FindNeighbors(ovary.VCD.Ptprc.pos, dims = 1:pcs)
ovary.VCD.Ptprc.pos <- FindClusters(object = ovary.VCD.Ptprc.pos)
ovary.VCD.Ptprc.pos <- FindClusters(object = ovary.VCD.Ptprc.pos, resolution = 3)

################################
# SingleR - Immgen to annotate Ptprc positive cells
################################

DefaultAssay(ovary.VCD.Ptprc.pos) <- "RNA"

my.SingleCellExperiment.object <- as.SingleCellExperiment(ovary.VCD.Ptprc.pos)

immgen <- ImmGenData(ensembl = FALSE, cell.ont = c("all", "nonna", "none"))
immgen
# class: SummarizedExperiment 
# dim: 22134 830 
# metadata(0):
#        assays(1): logcounts
# rownames(22134): Zglp1 Vmn2r65 ... Tiparp Kdm1a
# rowData names(0):
#        colnames(830): GSM1136119_EA07068_260297_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_1.CEL GSM1136120_EA07068_260298_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_2.CEL
# ... GSM920654_EA07068_201214_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_1.CEL GSM920655_EA07068_201215_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_2.CEL
# colData names(3): label.main label.fine label.ont

table(immgen$label.main)
#   B cells      B cells, pro         Basophils                DC Endothelial cells       Eosinophils  Epithelial cells       Fibroblasts               ILC       Macrophages        Mast cells         Microglia 
#        79                 1                 6                88                20                 4                25                21                23                79                20                 3 
# Monocytes       Neutrophils          NK cells               NKT        Stem cells     Stromal cells           T cells               Tgd 
#        33                23                38                22                36                 7               231                71

# Drop cell types with less than 10 cells
immgen <- immgen[,immgen$label.main != 'B cells, pro']
immgen <- immgen[,immgen$label.main != 'Basophils']
immgen <- immgen[,immgen$label.main != 'Eosinophils']
immgen <- immgen[,immgen$label.main != 'Microglia']
immgen <- immgen[,immgen$label.main != 'Stromal cells']

my.singler.immgen = SingleR(test = my.SingleCellExperiment.object,
                            ref  = immgen,
                            assay.type.test = 1,
                            labels = immgen$label.main)

# Add Sample IDs to metadata
my.singler.immgen$meta.data$Library       =    ovary.VCD.Ptprc.pos@meta.data$Library                       # Library IDs

save(my.singler.immgen, file = paste(Sys.Date(),"10x_ovary_Benayoun_lab_VCD_Ptprc_pos_cells_SingleR_object_Immgen.RData",sep = "_"))

###########################
# Transfer cell annotations to Seurat object

# Transfer SingleR annotations to Seurat Object
ovary.VCD.Ptprc.pos[["SingleR.Immgen"]] <- my.singler.immgen$labels

table.SingleR.annotation <- table(ovary.VCD.Ptprc.pos@meta.data$seurat_clusters, ovary.VCD.Ptprc.pos@meta.data$SingleR.Immgen)

################################
# scSorter to annotate Ptprc positive cells
# https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html
################################

# Major cell types: Neutrophils, Macrophages, Monocytes, DC, NK cells, ILC, NKT, T cells (CD8 & CD4), B cells
# Use Immgen dataset to curate marker gene list

immgen
# class: SummarizedExperiment 
# dim: 22134 830 ?SingleCellExperiment
# metadata(0):
#   assays(1): logcounts
# rownames(22134): Zglp1 Vmn2r65 ... Tiparp Kdm1a
# rowData names(0):
#   colnames(830): GSM1136119_EA07068_260297_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_1.CEL GSM1136120_EA07068_260298_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_2.CEL ...
# GSM920654_EA07068_201214_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_1.CEL GSM920655_EA07068_201215_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_2.CEL
# colData names(3): label.main label.fine label.ont

se.immgen <- as(immgen, "SingleCellExperiment")
seurat.immgen <- as.Seurat(se.immgen, counts = NULL, data = "logcounts")
seurat.immgen <- SetIdent(seurat.immgen, value = "label.main")

seurat.immgen <- FindVariableFeatures(seurat.immgen, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat.immgen)
seurat.immgen <- ScaleData(seurat.immgen, features = all.genes)

Neutrophil.markers <- FindMarkers(seurat.immgen, ident.1 = "Neutrophils", min.pct = 0.25, only.pos = TRUE)
Macrophages.markers <- FindMarkers(seurat.immgen, ident.1 = "Macrophages", min.pct = 0.25, only.pos = TRUE)
Monocytes.markers <- FindMarkers(seurat.immgen, ident.1 = "Monocytes", min.pct = 0.25, only.pos = TRUE)
DC.markers <- FindMarkers(seurat.immgen, ident.1 = "DC", min.pct = 0.25, only.pos = TRUE)
NK.markers <- FindMarkers(seurat.immgen, ident.1 = "NK cells", min.pct = 0.25, only.pos = TRUE)
ILC.markers <- FindMarkers(seurat.immgen, ident.1 = "ILC", min.pct = 0.25, only.pos = TRUE)
NKTcell.markers <- FindMarkers(seurat.immgen, ident.1 = "NKT", min.pct = 0.25, only.pos = TRUE)
Tcell.markers <- FindMarkers(seurat.immgen, ident.1 = "T cells", min.pct = 0.25, only.pos = TRUE)
Bcell.markers <- FindMarkers(seurat.immgen, ident.1 = "B cells", min.pct = 0.25, only.pos = TRUE)

# rownames(head(Neutrophil.markers, n = 20))

# Curate list of canonical markers for each major cell type
# https://panglaodb.se/markers.html?cell_type=%27NK%20cells%27

my.immune.marker.list <- vector(mode = "list", 9)
names(my.immune.marker.list) <- c("Neutrophils", "Macrophages", "Monocytes", "DC", "NK cells", "ILC", "NKT", "T cells", "B cells")

my.immune.marker.list$Neutrophils <- unique(c(rownames(head(Neutrophil.markers, n = 20)), "Csf3r", "Trem1", "Ptgs2"))     #Ly6g, S100a8 missing from dataset
my.immune.marker.list$Macrophages <- unique(c(rownames(head(Macrophages.markers, n = 20)), "Cd68", "Fcgr1", "Naaa"))
my.immune.marker.list$Monocytes <- unique(c(rownames(head(Monocytes.markers, n = 20)), "Cfp", "Cd14", "Rhoc"))
my.immune.marker.list$DC <- unique(c(rownames(head(DC.markers, n = 20)), "Zbtb46", "Itgax", "Lamp3"))
my.immune.marker.list$'NK cells' <- unique(c(rownames(head(NK.markers, n = 20)), "Ly6c2", "Trdc", "Nkg7", "Klrd1"))
my.immune.marker.list$ILC <- unique(c(rownames(head(ILC.markers, n = 20)), "Il7r", "Gata3"))
my.immune.marker.list$NKT <- unique(c(rownames(head(NKTcell.markers, n = 20)), "Ncam1", "Il2rb", "Cd44"))
my.immune.marker.list$'T cells' <- unique(c(rownames(head(Tcell.markers, n = 20)), "Cd3e", "Trbc2", "Cd3d", "Cd3g"))
my.immune.marker.list$'B cells' <- unique(c(rownames(head(Bcell.markers, n = 20)), "Pxk", "Ms4a1", "Cd75"))

# Filter genes that do not exist in dataset

my.immune.marker.list$Neutrophils <- intersect(rownames(ovary.VCD.Ptprc.pos), my.immune.marker.list$Neutrophils)
my.immune.marker.list$Macrophages <- intersect(rownames(ovary.VCD.Ptprc.pos), my.immune.marker.list$Macrophages)
my.immune.marker.list$Monocytes <- intersect(rownames(ovary.VCD.Ptprc.pos), my.immune.marker.list$Monocytes)
my.immune.marker.list$DC <- intersect(rownames(ovary.VCD.Ptprc.pos), my.immune.marker.list$DC)
my.immune.marker.list$'NK cells' <- intersect(rownames(ovary.VCD.Ptprc.pos), my.immune.marker.list$'NK cells')
my.immune.marker.list$ILC <- intersect(rownames(ovary.VCD.Ptprc.pos), my.immune.marker.list$ILC)
my.immune.marker.list$NKT <- intersect(rownames(ovary.VCD.Ptprc.pos), my.immune.marker.list$NKT)
my.immune.marker.list$'T cells' <- intersect(rownames(ovary.VCD.Ptprc.pos), my.immune.marker.list$'T cells')
my.immune.marker.list$'B cells' <- intersect(rownames(ovary.VCD.Ptprc.pos), my.immune.marker.list$'B cells')

### Generate anno table
anno<- data.frame(cbind(rep(names(my.immune.marker.list)[1],length(my.immune.marker.list[[1]])), my.immune.marker.list[[1]], rep(2,length(my.immune.marker.list[[1]]))))
colnames(anno) <- c("Type", "Marker", "Weight")

for (i in 2:length(my.immune.marker.list)) {
  anno <- rbind(anno,
                cbind("Type" = rep(names(my.immune.marker.list)[i],length(my.immune.marker.list[[i]])),
                      "Marker" = my.immune.marker.list[[i]],
                      "Weight" = rep(2,length(my.immune.marker.list[[i]]))))
}

# Manually increase weight for canonical markers for well defined cell types
anno[anno$Marker == "Csf3r",]$Weight <- 6 # Neutrophils
anno[anno$Marker == "Trem1",]$Weight <- 6 # Neutrophils

anno[anno$Marker == "Cd68",]$Weight <- 6 # Macrophages

anno[anno$Marker == "Cfp",]$Weight <- 6 # Monocytes
anno[anno$Marker == "Cd14",]$Weight <- 6 # Monocytes

anno[anno$Marker == "Zbtb46",]$Weight <- 6 # DC

anno[anno$Marker == "Ly6c2",]$Weight <- 6 # NK cells
anno[anno$Marker == "Trdc",]$Weight <- 6 # NK cells
anno[anno$Marker == "Nkg7",]$Weight <- 6 # NK cells

anno[anno$Marker == "Il7r",]$Weight <- 6 # ILC
anno[anno$Marker == "Gata3",]$Weight <- 6 # ILC

anno[anno$Marker == "Cd3e",]$Weight <- 6 # NKT
anno[anno$Marker == "Ly6c2",]$Weight <- 6 # NKT

anno[anno$Marker == "Cd3e",]$Weight <- 6 # T cells
anno[anno$Marker == "Trbc2",]$Weight <- 6 # T cells

anno[anno$Marker == "Cd19",]$Weight <- 6 # B cells
anno[anno$Marker == "Pxk",]$Weight <- 6 # B cells

# Pre-process data for scSorter

DefaultAssay(ovary.VCD.Ptprc.pos) <- "SCT"

topgenes <- head(VariableFeatures(ovary.VCD.Ptprc.pos), 2500)
expr <- GetAssayData(ovary.VCD.Ptprc.pos)
topgene_filter <- rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.05
topgenes <- topgenes[topgene_filter]

# At last, we subset the preprocessed expression data. Now, we are ready to run scSorter.
picked_genes = unique(c(anno$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

# Run scSorter
rts <- scSorter(expr, anno)

# Transfer scSorter annotations to Seurat Object
ovary.VCD.Ptprc.pos[["scSorter.Immgen"]] <- rts$Pred_Type

table.scSorter.annotation <- table(ovary.VCD.Ptprc.pos@meta.data$seurat_clusters, ovary.VCD.Ptprc.pos@meta.data$scSorter.Immgen)

################################
# sc-type to annotate Ptprc positive cells
# https://github.com/IanevskiAleksandr/sc-type
################################

# Use same marker list as for scSorter
my.immune.marker.list
names(my.immune.marker.list)
# [1] "Neutrophils" "Macrophages" "Monocytes"   "DC"          "NK cells"    "ILC"         "NKT"         "T cells"     "B cells"   

cell_markers.tmp = data.frame("cellName" = names(my.immune.marker.list),
                              "geneSymbolmore1" = rep("",length(my.immune.marker.list)),
                              "geneSymbolmore2" = rep("",length(my.immune.marker.list)))

# grab positive markers automatically
for (i in 1:length(my.immune.marker.list)) {
  cell_markers.tmp$geneSymbolmore1[i] <- paste0(my.immune.marker.list[[cell_markers.tmp$cellName[i]]], collapse = ",")
}

# manually edit negative markers
neg.markers <- vector(mode = "list", length = length(my.immune.marker.list))
names(neg.markers) <- names(my.immune.marker.list)

neg.markers[names(neg.markers) != "Neutrophils" ] <- lapply(neg.markers[names(neg.markers) != "Neutrophils" ], function (x){append(x,"Mpx")}) 
neg.markers[names(neg.markers) != "T cells"     ] <- lapply(neg.markers[names(neg.markers) != "T cells"     ], function (x){append(x,"Zap70")}) # Tang, JEM 2017
neg.markers[names(neg.markers) != "T cells"     ] <- lapply(neg.markers[names(neg.markers) != "T cells"     ], function (x){append(x,"Cd4")}) 
neg.markers[names(neg.markers) != "T cells"     ] <- lapply(neg.markers[names(neg.markers) != "T cells"     ], function (x){append(x,"Cd8a")})
neg.markers[names(neg.markers) != "T cells"     ] <- lapply(neg.markers[names(neg.markers) != "T cells"     ], function (x){append(x,"Cd8b1")})
neg.markers[names(neg.markers) != "B cells"     ] <- lapply(neg.markers[names(neg.markers) != "B cells"     ], function (x){append(x,"Cd79a")}) 
neg.markers[names(neg.markers) != "Macrophages" ] <- lapply(neg.markers[names(neg.markers) != "Macrophages" ], function (x){append(x,"Csf1r")})
neg.markers[names(neg.markers) != c("Neutrophils", "Macrophages")] <- lapply(neg.markers[names(neg.markers) != c("Neutrophils", "Macrophages")], function (x){append(x, "Lyz1")}) #lyz cannot be outside of macrophages and neutrophils

# populate negative markers
for (i in 1:length(my.immune.marker.list)) {
  cell_markers.tmp$geneSymbolmore2[i] <- paste0(neg.markers[[cell_markers.tmp$cellName[i]]], collapse = ",")
}

###################################################################################################
######  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   #####
##### Modify code from "gene_sets_prepare" of scType to marker database in R
# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021

source("/Users/minhookim/Programs/sc-type-1.0/R/sctype_score_.R")

gene_sets_prepare_mOvary <- function(cell_markers.tmp){
  cell_markers = cell_markers.tmp
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}

# prepare gene sets
gs_list = gene_sets_prepare_mOvary(cell_markers.tmp)

# Finally, let's assign cell types to each cluster:
# get cell-type by cell matrix
es.max  = sctype_score(scRNAseqData = ovary.VCD.Ptprc.pos[["SCT"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(ovary.VCD.Ptprc.pos@meta.data$SCT_snn_res.4), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(ovary.VCD.Ptprc.pos@meta.data[ovary.VCD.Ptprc.pos@meta.data$SCT_snn_res.4==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(ovary.VCD.Ptprc.pos@meta.data$SCT_snn_res.4==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
# print(sctype_scores[,1:3])

ovary.VCD.Ptprc.pos@meta.data$scType_raw = ""

for(i in 1:ncol(es.max)){
  ovary.VCD.Ptprc.pos@meta.data$scType_raw[rownames(ovary.VCD.Ptprc.pos@meta.data) %in% colnames(es.max)[i] ] = as.character(rownames(es.max)[which.max(es.max[,i])])
}

table(ovary.VCD.Ptprc.pos@meta.data$scType_raw)
# B cells          DC         ILC Macrophages   Monocytes Neutrophils    NK cells         NKT     T cells 
#     530         351         739         635         215         266         869         700        1626 

table.scType.annotation <- table(ovary.VCD.Ptprc.pos@meta.data$seurat_clusters, ovary.VCD.Ptprc.pos@meta.data$scType_raw)

save(ovary.VCD.Ptprc.pos, file = paste(Sys.Date(),"10x_ovary_Benayoun_lab_VCD_Ptprc_positive_Seurat_object_with_cell_type_annotations_intermediate.RData",sep = "_"))

################################
# Assess cell type annotation by SingleR, scSorter and scType
################################

# SingleR

table.SingleR.annotation <- as.data.frame.matrix(table.SingleR.annotation)

# Manipulate data structure

df_long <- table.SingleR.annotation %>% 
  mutate(cluster = row_number()) %>%  # create a cluster column based on row numbers
  pivot_longer(cols = -cluster, names_to = "cell_type", values_to = "count")

df_count <- df_long %>%
  group_by(cluster, cell_type) %>%
  summarise(count = sum(count), .groups = "keep")

df_max.SingleR <- df_count %>%
  group_by(cluster) %>%
  slice(which.max(count))

# scSorter

table.scSorter.annotation <- as.data.frame.matrix(table.scSorter.annotation)

# Manipulate data structure

df_long <- table.scSorter.annotation %>% 
  mutate(cluster = row_number()) %>%  # create a cluster column based on row numbers
  pivot_longer(cols = -cluster, names_to = "cell_type", values_to = "count")

df_count <- df_long %>%
  group_by(cluster, cell_type) %>%
  summarise(count = sum(count), .groups = "keep")

df_max.scSorter <- df_count %>%
  group_by(cluster) %>%
  slice(which.max(count))

# scType

table.scType.annotation <- as.data.frame.matrix(table.scType.annotation)

# Manipulate data structure

df_long <- table.scType.annotation %>% 
  mutate(cluster = row_number()) %>%  # create a cluster column based on row numbers
  pivot_longer(cols = -cluster, names_to = "cell_type", values_to = "count")

df_count <- df_long %>%
  group_by(cluster, cell_type) %>%
  summarise(count = sum(count), .groups = "keep")

df_max.scType <- df_count %>%
  group_by(cluster) %>%
  slice(which.max(count))

# Combine results

# Check if clusters are in the same order
check_order <- identical(df_max.SingleR$cluster, df_max.scSorter$cluster) && 
  identical(df_max.SingleR$cluster, df_max.scType$cluster)

print(check_order)    # TRUE

annotation.results <- cbind(df_max.SingleR, df_max.scSorter, df_max.scType)
annotation.results <- annotation.results[,c(1,2,3,5,6,8,9)]
colnames(annotation.results) <- c("Cluster", "SingleR_Celltype", "SingleR_Cellcount", "scSorter_Celltype", "scSorter_Cellcount", "scType_Celltype", "scType_Cellcount")

save(annotation.results, file = paste(Sys.Date(),"10x_ovary_Benayoun_lab_VCD_Ptprc_pos_cells_annotation_results.RData",sep = "_"))
write.table(annotation.results, file = paste(Sys.Date(),"10x_ovary_Benayoun_lab_VCD_Ptprc_pos_cells_annotation_results.txt",sep = "_"), quote = FALSE)

# Manual marker gene analysis

ovary.VCD.Ptprc.pos <- SetIdent(ovary.VCD.Ptprc.pos, value = "seurat_clusters")

pdf(paste0(Sys.Date(), "_10x_ovary_VCD_Ptprc_positive_cells_Vlnplot_marker_genes_clusters_1_to_7.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.VCD.Ptprc.pos, subset = seurat_clusters %in% c(0:7)), 
        features = c("S100a8", "Adgre1", "C1qa", "Flt3", 
                     "Il2ra", "Ccdc184", "Klrb1c", "Fam184b", 
                     "Ncr1", "Gzma", "Gzmb", "Nkg7", 
                     "Cd3e", "Cd8a", "Cd4", "Cd19"))
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_VCD_Ptprc_positive_cells_Vlnplot_marker_genes_clusters_8_to_14.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.VCD.Ptprc.pos, subset = seurat_clusters %in% c(8:14)), 
        features = c("S100a8", "Adgre1", "C1qa", "Flt3", 
                     "Il2ra", "Ccdc184", "Klrb1c", "Fam184b", 
                     "Ncr1", "Gzma", "Gzmb", "Nkg7", 
                     "Cd3e", "Cd8a", "Cd4", "Cd19"))
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_VCD_Ptprc_positive_cells_Vlnplot_marker_genes_clusters_15_to_21.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.VCD.Ptprc.pos, subset = seurat_clusters %in% c(15:21)), 
        features = c("S100a8", "Adgre1", "C1qa", "Flt3", 
                     "Il2ra", "Ccdc184", "Klrb1c", "Fam184b", 
                     "Ncr1", "Gzma", "Gzmb", "Nkg7", 
                     "Cd3e", "Cd8a", "Cd4", "Cd19"))
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_VCD_Ptprc_positive_cells_Vlnplot_marker_genes_clusters_22_to_28.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.VCD.Ptprc.pos, subset = seurat_clusters %in% c(22:28)), 
        features = c("S100a8", "Adgre1", "C1qa", "Flt3", 
                     "Il2ra", "Ccdc184", "Klrb1c", "Fam184b", 
                     "Ncr1", "Gzma", "Gzmb", "Nkg7", 
                     "Cd3e", "Cd8a", "Cd4", "Cd19"))
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_VCD_Ptprc_positive_cells_Vlnplot_marker_genes_clusters_29_to_31.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.VCD.Ptprc.pos, subset = seurat_clusters %in% c(29:31)), 
        features = c("S100a8", "Ly6c2", "C1qa", "Flt3", 
                     "Il2ra", "Ccdc184", "Klrb1c", "Fam184b", 
                     "Ncr1", "Gzma", "Gzmb", "Nkg7", 
                     "Cd3e", "Cd8a", "Cd4", "Cd19"))
dev.off()

################################
# Transfer final annotations
################################

ovary.VCD.Ptprc.pos <- SetIdent(ovary.VCD.Ptprc.pos, value = "seurat_clusters")

VlnPlot(ovary.VCD.Ptprc.pos, features = c("nCount_RNA"))
VlnPlot(ovary.VCD.Ptprc.pos, features = c("nFeature_RNA"))
VlnPlot(ovary.VCD.Ptprc.pos, features = c("percent.mito"))
VlnPlot(ovary.VCD.Ptprc.pos, features = c("decontX_contamination"))

ovary.VCD.Ptprc.pos@meta.data$celltype.level2 <- rep("NA",dim(ovary.VCD.Ptprc.pos@meta.data)[1])
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 0] <- "CD8 NKT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 1] <- "DNT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 2] <- "DNT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 3] <- "CD8 NKT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 4] <- "NK"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 5] <- "Monocyte"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 6] <- "NK"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 7] <- "DPT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 8] <- "CD8 NKT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 9] <- "B"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 10] <- "NKT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 11] <- "DNT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 12] <- "B"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 13] <- "NK"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 14] <- "DC"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 15] <- "NKT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 16] <- "NKT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 17] <- "CD4 T"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 18] <- "NK"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 19] <- "DNT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 20] <- "ILC"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 21] <- "DC"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 22] <- "DNT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 23] <- "CD8 NKT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 24] <- "NKT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 25] <- "Neutrophil"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 26] <- "Macrophage"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 27] <- "Macrophage"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 28] <- "CD8 NKT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 29] <- "DNT"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 30] <- "B"
ovary.VCD.Ptprc.pos@meta.data$celltype.level2[ovary.VCD.Ptprc.pos@meta.data$seurat_clusters %in% 31] <- "Macrophage"

ovary.VCD.Ptprc.pos <- SetIdent(ovary.VCD.Ptprc.pos, value = "celltype.level2")
levels(ovary.VCD.Ptprc.pos) <- c("Neutrophil", "Monocyte", "Macrophage" , "DC", "ILC", "NK", "NKT", "CD8 NKT", "CD4 T", "DNT", "DPT", "B")

save(ovary.VCD.Ptprc.pos, file = paste(Sys.Date(),"10x_ovary_Benayoun_lab_VCD_Ptprc_positive_cells_Seurat_object_with_final_annotation.RData",sep = "_"))

################################################################################
# 3.Ptprc negative sub-celltype annotation: 
#                                           1. Manual marker annotation
#                                           2. SingleR / scType / scSorter
################################################################################

################################
# Subset Ptprc positive cells and re-cluster
################################

ovary.VCD.Ptprc.neg <- subset(ovary.VCD, subset = celltype.level1 == "Ptprc.neg")

ovary.VCD.Ptprc.neg
# An object of class Seurat 
# 48841 features across 41455 samples within 2 assays 
# Active assay: RNA (24485 features, 0 variable features)
# 1 other assay present: SCT
# 3 dimensional reductions calculated: pca, harmony, umap

DefaultAssay(ovary.VCD.Ptprc.pos) <- "RNA"

ovary.VCD.Ptprc.neg <- DietSeurat(ovary.VCD.Ptprc.neg,
                                  counts = TRUE, 
                                  data = TRUE,
                                  scale.data = FALSE,
                                  features = NULL,
                                  assays = NULL,
                                  dimreducs = NULL,
                                  graphs = NULL)

ovary.VCD.Ptprc.neg[['nCount_SCT']] <- NULL
ovary.VCD.Ptprc.neg[['nFeature_SCT']] <- NULL
ovary.VCD.Ptprc.neg[['SCT_snn_res.2']] <- NULL
ovary.VCD.Ptprc.neg[['seurat_clusters']] <- NULL
ovary.VCD.Ptprc.neg[['DoubletFinder']] <- NULL
ovary.VCD.Ptprc.neg[['scds_hybrid']] <- NULL

########## Normalize, SCT, dimreduc data ##########

ovary.VCD.Ptprc.neg <- NormalizeData(ovary.VCD.Ptprc.neg, normalization.method = "LogNormalize", scale.factor = 10000)
ovary.VCD.Ptprc.neg <- SCTransform(object = ovary.VCD.Ptprc.neg, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))

save(ovary.VCD.Ptprc.neg, file = paste0(Sys.Date(),"_10x_ovary_VCD_Ptprc_neg_cells_Seurat_object_clean_postSCT.RData"))

ovary.VCD.Ptprc.neg <- RunPCA(ovary.VCD.Ptprc.neg, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_10x_ovary_VCD_Ptprc_neg_cells_ElbowPlot.pdf"))
ElbowPlot(ovary.VCD.Ptprc.neg, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- ovary.VCD.Ptprc.neg[["pca"]]@stdev / sum(ovary.VCD.Ptprc.neg[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 42

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 16

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 16

# Based on these metrics, first 16 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_10x_ovary_VCD_Ptprc_pos_cells_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

###############################################################################
# run dimensionality reduction algorithm
ovary.VCD.Ptprc.neg <- RunUMAP(ovary.VCD.Ptprc.neg, dims = 1:pcs)
ovary.VCD.Ptprc.neg <- FindNeighbors(ovary.VCD.Ptprc.neg, dims = 1:pcs)
ovary.VCD.Ptprc.neg <- FindClusters(object = ovary.VCD.Ptprc.neg)
ovary.VCD.Ptprc.neg <- FindClusters(object = ovary.VCD.Ptprc.neg, resolution = 3)

ovary.VCD.Ptprc.neg <- SetIdent(ovary.VCD.Ptprc.neg, value = "seurat_clusters")

################################
# SingleR - Public datasets to annotate ovarian cells - PMID:36205477, PMID:38200272
################################

###############################
# PMID:36205477
###############################

# Prep reference file

cycling.ovary.seurat <- readRDS("/Volumes/OIProject_II/1_R/SingleR_Ref/ovary_0.rds")
cycling.ovary.seurat
# An object of class Seurat 
# 23657 features across 34712 samples within 1 assay 
# Active assay: RNA (23657 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

cycling.ovary.seurat@meta.data$celltype <- rep('NA')
cycling.ovary.seurat@meta.data$celltype[grep("GC_", cycling.ovary.seurat@meta.data$Level1)] <- "Granulosa"
cycling.ovary.seurat@meta.data$celltype[grep("M_", cycling.ovary.seurat@meta.data$Level1)] <- "Mesenchyme"
cycling.ovary.seurat@meta.data$celltype[grep("Stroma", cycling.ovary.seurat@meta.data$Level1)] <- "Stroma"
cycling.ovary.seurat@meta.data$celltype[grep("Theca", cycling.ovary.seurat@meta.data$Level1)] <- "Theca"
cycling.ovary.seurat@meta.data$celltype[grep("Smooth muscle", cycling.ovary.seurat@meta.data$Level1)] <- "Smooth muscle"
cycling.ovary.seurat@meta.data$celltype[grep("I_", cycling.ovary.seurat@meta.data$Level1)] <- "Immune"
cycling.ovary.seurat@meta.data$celltype[grep("EN_B", cycling.ovary.seurat@meta.data$Level1)] <- "BEC"
cycling.ovary.seurat@meta.data$celltype[grep("EN_L", cycling.ovary.seurat@meta.data$Level1)] <- "LEC"
cycling.ovary.seurat@meta.data$celltype[grep("Epithelium", cycling.ovary.seurat@meta.data$Level1)] <- "Epithelial"
cycling.ovary.seurat@meta.data$celltype[grep("Oocyte", cycling.ovary.seurat@meta.data$Level1)] <- "Oocyte"

cycling.ovary.seurat.2 <- subset(cycling.ovary.seurat, subset = celltype == "Immune", invert = TRUE)

cycle.data.sce <- as.SingleCellExperiment(cycling.ovary.seurat.2)
dim(cycle.data.sce)
# [1] 23657 33063

# Prep test data
DefaultAssay(ovary.VCD.Ptprc.neg) <- "RNA"
my.SingleCellExperiment.object <- as.SingleCellExperiment(ovary.VCD.Ptprc.neg)

# Find genes in common to the data sets and limit analysis to these
commonGenes <- intersect(rownames(cycle.data.sce), rownames(my.SingleCellExperiment.object))
# length(commonGenes)

my.cycling.ovary.ref <- cycle.data.sce[commonGenes,]
# dim(my.cycling.ovary.ref)

my.SingleCellExperiment.object <- my.SingleCellExperiment.object[commonGenes,]
# dim(my.SingleCellExperiment.object)

my.trained.cycling <- trainSingleR(my.cycling.ovary.ref, my.cycling.ovary.ref$celltype)

my.singleR.cycling <- classifySingleR(test = my.SingleCellExperiment.object, trained = my.trained.cycling)

save(my.trained.cycling, file = paste(Sys.Date(),"10X_ovary_Benayoun_VCD_Ptprc_negative_cells_SingleR_cycing_data_trained_classifier.RData",sep = "_"))
save(my.singleR.cycling, file = paste(Sys.Date(),"10X_ovary_Benayoun_VCD_Ptprc_negative_cells_SingleR_cycing_data_trained_result.RData",sep = "_"))

# Transfer cell type annotation to Seurat
ovary.VCD.Ptprc.neg[["SingleR.cycling"]] <- my.singleR.cycling$labels

table.SingleR.cycling.annotation <- table(ovary.VCD.Ptprc.neg@meta.data$seurat_clusters, ovary.VCD.Ptprc.neg@meta.data$SingleR.cycling)

###############################
# PMID:38200272
###############################

GSE232309.ovary.seurat <- readRDS("/Volumes/OIProject_II/1_R/SingleR_Ref/GSE232309_age.combined.RDS")
GSE232309.ovary.seurat
# An object of class Seurat 
# 20264 features across 14349 samples within 1 assay 
# Active assay: RNA (20264 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

GSE232309.ovary.seurat@meta.data$celltype <- rep('NA')
GSE232309.ovary.seurat@meta.data$celltype[grep("Granulosa", GSE232309.ovary.seurat@meta.data$cluster.names)] <- "Granulosa"
GSE232309.ovary.seurat@meta.data$celltype[grep("Luteal", GSE232309.ovary.seurat@meta.data$cluster.names)] <- "Granulosa"
GSE232309.ovary.seurat@meta.data$celltype[grep("Stroma", GSE232309.ovary.seurat@meta.data$cluster.names)] <- "Stroma"
GSE232309.ovary.seurat@meta.data$celltype[grep("Theca", GSE232309.ovary.seurat@meta.data$cluster.names)] <- "Theca"
GSE232309.ovary.seurat@meta.data$celltype[grep("Epithelium", GSE232309.ovary.seurat@meta.data$cluster.names)] <- "Epithelial"
GSE232309.ovary.seurat@meta.data$celltype[grep("Endothelium", GSE232309.ovary.seurat@meta.data$cluster.names)] <- "Endothelial"
GSE232309.ovary.seurat@meta.data$celltype[grep("Phagocytes", GSE232309.ovary.seurat@meta.data$cluster.names)] <- "Immune"
GSE232309.ovary.seurat@meta.data$celltype[grep("Lymphocytes", GSE232309.ovary.seurat@meta.data$cluster.names)] <- "Immune"
GSE232309.ovary.seurat@meta.data$celltype[grep("Oocyte", GSE232309.ovary.seurat@meta.data$cluster.names)] <- "Oocyte"

GSE232309.ovary.seurat.2 <- subset(GSE232309.ovary.seurat, subset = celltype %in% c("Immune"), invert = TRUE)

GSE232309.data.sce <- as.SingleCellExperiment(GSE232309.ovary.seurat.2)
dim(GSE232309.data.sce)
# [1] 20264 12320

# Prep test data
DefaultAssay(ovary.VCD.Ptprc.neg) <- "RNA"
my.SingleCellExperiment.object <- as.SingleCellExperiment(ovary.VCD.Ptprc.neg)

# Find genes in common to the data sets and limit analysis to these
commonGenes <- intersect(rownames(GSE232309.data.sce), rownames(my.SingleCellExperiment.object))
# length(commonGenes)

my.GSE232309.ovary.ref <- GSE232309.data.sce[commonGenes,]
# dim(my.cycling.ovary.ref)

my.SingleCellExperiment.object <- my.SingleCellExperiment.object[commonGenes,]
# dim(my.SingleCellExperiment.object)

my.trained.GSE232309 <- trainSingleR(my.GSE232309.ovary.ref, my.GSE232309.ovary.ref$celltype)

my.singleR.GSE232309 <- classifySingleR(test = my.SingleCellExperiment.object, trained = my.trained.GSE232309)

save(my.trained.GSE232309, file = paste(Sys.Date(),"10X_ovary_Benayoun_VCD_Ptprc_negative_cells_SingleR_cycing_data_trained_classifier.RData",sep = "_"))
save(my.singleR.GSE232309, file = paste(Sys.Date(),"10X_ovary_Benayoun_VCD_Ptprc_negative_cells_SingleR_cycing_data_trained_result.RData",sep = "_"))

# Transfer cell type annotation to Seurat
ovary.VCD.Ptprc.neg[["SingleR.GSE232309"]] <- my.singleR.GSE232309$labels

table.SingleR.GSE232309.annotation <- table(ovary.VCD.Ptprc.neg@meta.data$seurat_clusters, ovary.VCD.Ptprc.neg@meta.data$SingleR.GSE232309)

################################
# scSorter to annotate Ptprc negative cells
# https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html
################################

# Use curated marker gene list
# Marker genes derived from cycling data + manually added from literature

cycling.ovary.seurat.2 <- NormalizeData(cycling.ovary.seurat.2, normalization.method = "LogNormalize", scale.factor = 10000)
cycling.ovary.seurat.2 <- SCTransform(object = cycling.ovary.seurat.2, vars.to.regress = c("nGene", "fraction.mito", "mouse"))
cycling.ovary.seurat.2 <- RunPCA(cycling.ovary.seurat.2, npcs = 50)
cycling.ovary.seurat.2 <- RunUMAP(cycling.ovary.seurat.2, dims = 1:20)
cycling.ovary.seurat.2 <- FindNeighbors(cycling.ovary.seurat.2, dims = 1:20)
cycling.ovary.seurat.2 <- FindClusters(object = cycling.ovary.seurat.2)

cycling.ovary.seurat.2
# An object of class Seurat 
# 46196 features across 33063 samples within 2 assays 
# Active assay: SCT (22539 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

cycling.ovary.seurat.2 <- SetIdent(cycling.ovary.seurat.2, value = "celltype")
cycling.ovary.markers <- FindAllMarkers(cycling.ovary.seurat.2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cycling.ovary.markers.top20 <- cycling.ovary.markers %>% 
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

my.ovarian.marker.list <- vector(mode = "list", 9)
names(my.ovarian.marker.list) <- c("Oocyte",
                                   "Granulosa",
                                   "Theca",
                                   "Stroma",
                                   "Mesenchyme",
                                   "BEC",
                                   "LEC",
                                   "Epithelial")

# Curate marker gene list (top 20 marker genes from cycling data + additional markers)
my.ovarian.marker.list$Oocyte <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "Oocyte"], "Ddx4", "DazI"))
my.ovarian.marker.list$Granulosa <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "Granulosa"], "Amh", "Fshr", "Cyp19a1", "Foxl2", "Bex1", "Cnmd", "Inha", "Gja1", "Amhr2", "KitI", "Cdh2", "Star"))
my.ovarian.marker.list$Theca <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "Theca"], "Top2a", "Lhcgr", "Dcn", "Mdk", "Star", "Cyp11a1", "Cyp17a1", "Col1a2", "Col3a1", "Ogn", "Fbln5"))
my.ovarian.marker.list$Stroma <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "Stroma"],  "Mdk", "Nr2f2", "Col6a3", "Lum", "Col1a1", "Pdgfra", "Nr2f2"))
my.ovarian.marker.list$Mesenchyme <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "Mesenchyme"],  "Acta2", "Mcam", "Rgs5", "Notch3"))
my.ovarian.marker.list$BEC <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "BEC"],  "Flt1", "Vwf"))
my.ovarian.marker.list$LEC <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "LEC"],  "Flt4", "Prox1"))
my.ovarian.marker.list$Epithelial <- unique(c(cycling.ovary.markers.top20$gene[cycling.ovary.markers.top20$cluster == "Epithelial"],  "Cldn1", "Pax8"))

save(my.ovarian.marker.list, file = paste(Sys.Date(),"Ovarian_non-immune_cell_markers_for_scSorter.RData",sep = "_"))


# Filter genes that do not exist in dataset

my.ovarian.marker.list$Oocyte <- intersect(rownames(ovary.VCD.Ptprc.neg), my.ovarian.marker.list$Oocyte)
my.ovarian.marker.list$Granulosa <- intersect(rownames(ovary.VCD.Ptprc.neg), my.ovarian.marker.list$Granulosa)
my.ovarian.marker.list$Theca <- intersect(rownames(ovary.VCD.Ptprc.neg), my.ovarian.marker.list$Theca)
my.ovarian.marker.list$Stroma <- intersect(rownames(ovary.VCD.Ptprc.neg), my.ovarian.marker.list$Stroma)
my.ovarian.marker.list$Mesenchyme <- intersect(rownames(ovary.VCD.Ptprc.neg), my.ovarian.marker.list$Mesenchyme)
my.ovarian.marker.list$BEC <- intersect(rownames(ovary.VCD.Ptprc.neg), my.ovarian.marker.list$BEC)
my.ovarian.marker.list$LEC <- intersect(rownames(ovary.VCD.Ptprc.neg), my.ovarian.marker.list$LEC)
my.ovarian.marker.list$Epithelial <- intersect(rownames(ovary.VCD.Ptprc.neg), my.ovarian.marker.list$Epithelial)

### Generate anno table
anno<- data.frame(cbind(rep(names(my.ovarian.marker.list)[1],length(my.ovarian.marker.list[[1]])), my.ovarian.marker.list[[1]], rep(2,length(my.ovarian.marker.list[[1]]))))
colnames(anno) <- c("Type", "Marker", "Weight")

for (i in 2:length(my.ovarian.marker.list)) {
  anno <- rbind(anno,
                cbind("Type" = rep(names(my.ovarian.marker.list)[i],length(my.ovarian.marker.list[[i]])),
                      "Marker" = my.ovarian.marker.list[[i]],
                      "Weight" = rep(2,length(my.ovarian.marker.list[[i]]))))
}

# Manually increase weight for canonical markers for well defined cell types
anno[anno$Marker == "Ddx4",]$Weight <- 6
anno[anno$Marker == "Cyp19a1",]$Weight <- 6
anno[anno$Marker == "Amh",]$Weight <- 6
anno[anno$Marker == "Fshr",]$Weight <- 6
anno[anno$Marker == "Foxl2",]$Weight <- 6
anno[anno$Marker == "Cyp11a1",]$Weight <- 6
anno[anno$Marker == "Cyp17a1",]$Weight <- 6
anno[anno$Marker == "Nr2f2",]$Weight <- 6
anno[anno$Marker == "Col1a1",]$Weight <- 6
anno[anno$Marker == "Rgs5",]$Weight <- 6
anno[anno$Marker == "Flt1",]$Weight <- 6
anno[anno$Marker == "Flt4",]$Weight <- 6

# Pre-process data for scSorter
topgenes <- head(VariableFeatures(ovary.VCD.Ptprc.neg), 2500)
expr <- GetAssayData(ovary.VCD.Ptprc.neg)
topgene_filter <- rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.05
topgenes <- topgenes[topgene_filter]

# At last, we subset the preprocessed expression data. Now, we are ready to run scSorter.
picked_genes = unique(c(anno$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

# Run scSorter
rts <- scSorter(expr, anno)

# Transfer scSorter annotations to Seurat Object
ovary.VCD.Ptprc.neg[["scSorter.labels"]] <- rts$Pred_Type

# table(ovary.VCD.Ptprc.neg@meta.data$seurat_clusters, ovary.VCD.Ptprc.neg@meta.data$scSorter.labels)

table.scSorter.annotation <- table(ovary.VCD.Ptprc.neg@meta.data$seurat_clusters, ovary.VCD.Ptprc.neg@meta.data$scSorter.labels)

ovary.VCD.Ptprc.neg <- SetIdent(ovary.VCD.Ptprc.neg, value = "scSorter.labels")
levels(ovary.VCD.Ptprc.neg) <- c("Oocyte",
                                 "Granulosa",
                                 "Theca",
                                 "Stroma",
                                 "Mesenchyme",
                                 "Pericyte",
                                 "BEC",
                                 "LEC",
                                 "Epithelial",
                                 "Unknown")

################################
# sc-type to annotate Ptprc negative cells
# https://github.com/IanevskiAleksandr/sc-type
################################

# Use same marker list as for scSorter
my.ovarian.marker.list
names(my.ovarian.marker.list)
# [1] "Oocyte"     "Granulosa"  "Theca"      "Stroma"     "Mesenchyme" "BEC"        "LEC"        "Epithelial"

cell_markers.tmp = data.frame("cellName" = names(my.ovarian.marker.list),
                              "geneSymbolmore1" = rep("",length(my.ovarian.marker.list)),
                              "geneSymbolmore2" = rep("",length(my.ovarian.marker.list)))

# grab positive markers automatically
for (i in 1:length(my.ovarian.marker.list)) {
  cell_markers.tmp$geneSymbolmore1[i] <- paste0(my.ovarian.marker.list[[cell_markers.tmp$cellName[i]]], collapse = ",")
}

# manually edit negative markers
neg.markers <- vector(mode = "list", length = length(my.ovarian.marker.list))
names(neg.markers) <- names(my.ovarian.marker.list)

# populate negative markers
for (i in 1:length(my.ovarian.marker.list)) {
  cell_markers.tmp$geneSymbolmore2[i] <- paste0(neg.markers[[cell_markers.tmp$cellName[i]]], collapse = ",")
}

###################################################################################################
######  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   #####
##### Modify code from "gene_sets_prepare" of scType to marker database in R
# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021

source("/Users/minhookim/Programs/sc-type-1.0/R/sctype_score_.R")

gene_sets_prepare_mOvary <- function(cell_markers.tmp){
  cell_markers = cell_markers.tmp
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}

# prepare gene sets
gs_list = gene_sets_prepare_mOvary(cell_markers.tmp)

# Finally, let's assign cell types to each cluster:
# get cell-type by cell matrix
es.max  = sctype_score(scRNAseqData = ovary.VCD.Ptprc.neg[["SCT"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(ovary.VCD.Ptprc.neg@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(ovary.VCD.Ptprc.neg@meta.data[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(ovary.VCD.Ptprc.neg@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
# print(sctype_scores[,1:3])

ovary.VCD.Ptprc.neg@meta.data$scType_raw = ""

for(i in 1:ncol(es.max)){
  ovary.VCD.Ptprc.neg@meta.data$scType_raw[rownames(ovary.VCD.Ptprc.neg@meta.data) %in% colnames(es.max)[i] ] = as.character(rownames(es.max)[which.max(es.max[,i])])
}

table.scType.annotation <- table(ovary.VCD.Ptprc.neg@meta.data$seurat_clusters, ovary.VCD.Ptprc.neg@meta.data$scType_raw)

save(ovary.VCD.Ptprc.neg, file = paste(Sys.Date(),"10x_ovary_Benayoun_lab_VCD_Ptprc_negativetive_Seurat_object_with_cell_type_annotations_intermediate.RData",sep = "_"))

################################
# Assess cell type annotation by SingleR, scSorter and scType
################################

# SingleR - cycling

table.SingleR.cycling.annotation <- as.data.frame.matrix(table.SingleR.cycling.annotation)

# Manipulate data structure

df_long <- table.SingleR.cycling.annotation %>% 
  mutate(cluster = row_number()) %>%  # create a cluster column based on row numbers
  pivot_longer(cols = -cluster, names_to = "cell_type", values_to = "count")

df_count <- df_long %>%
  group_by(cluster, cell_type) %>%
  summarise(count = sum(count), .groups = "keep")

df_max.SingleR.cycling <- df_count %>%
  group_by(cluster) %>%
  slice(which.max(count))

# SingleR - GSE232309

table.SingleR.GSE232309.annotation <- as.data.frame.matrix(table.SingleR.GSE232309.annotation)

# Manipulate data structure

df_long <- table.SingleR.GSE232309.annotation %>% 
  mutate(cluster = row_number()) %>%  # create a cluster column based on row numbers
  pivot_longer(cols = -cluster, names_to = "cell_type", values_to = "count")

df_count <- df_long %>%
  group_by(cluster, cell_type) %>%
  summarise(count = sum(count), .groups = "keep")

df_max.SingleR.GSE232309 <- df_count %>%
  group_by(cluster) %>%
  slice(which.max(count))

# scSorter

table.scSorter.annotation <- as.data.frame.matrix(table.scSorter.annotation)

# Manipulate data structure

df_long <- table.scSorter.annotation %>% 
  mutate(cluster = row_number()) %>%  # create a cluster column based on row numbers
  pivot_longer(cols = -cluster, names_to = "cell_type", values_to = "count")

df_count <- df_long %>%
  group_by(cluster, cell_type) %>%
  summarise(count = sum(count), .groups = "keep")

df_max.scSorter <- df_count %>%
  group_by(cluster) %>%
  slice(which.max(count))

# scType

table.scType.annotation <- as.data.frame.matrix(table.scType.annotation)

# Manipulate data structure

df_long <- table.scType.annotation %>% 
  mutate(cluster = row_number()) %>%  # create a cluster column based on row numbers
  pivot_longer(cols = -cluster, names_to = "cell_type", values_to = "count")

df_count <- df_long %>%
  group_by(cluster, cell_type) %>%
  summarise(count = sum(count), .groups = "keep")

df_max.scType <- df_count %>%
  group_by(cluster) %>%
  slice(which.max(count))

# Combine results

# Check if clusters are in the same order
check_order <- identical(df_max.SingleR.cycling$cluster, df_max.SingleR.GSE232309$cluster) && 
  identical(df_max.SingleR.cycling$cluster, df_max.scSorter$cluster) &&
  identical(df_max.SingleR.cycling$cluster, df_max.scType$cluster)

print(check_order)    # TRUE

annotation.results <- cbind(df_max.SingleR.cycling, df_max.SingleR.GSE232309, df_max.scSorter, df_max.scType)
annotation.results <- annotation.results[,c(1,2,3,5,6,8,9,11,12)]
colnames(annotation.results) <- c("Cluster", "SingleR_cycling_Celltype", "SingleR_cycling_Cellcount", "SingleR_GSE232309_Celltype", "SingleR_GSE232309_Cellcount", "scSorter_Celltype", "scSorter_Cellcount", "scType_Celltype", "scType_Cellcount")

save(annotation.results, file = paste(Sys.Date(),"10x_ovary_Benayoun_lab_VCD_Ptprc_negative_cells_annotation_results.RData",sep = "_"))
write.table(annotation.results, file = paste(Sys.Date(),"10x_ovary_Benayoun_lab_VCD_Ptprc_negative_cells_annotation_results.txt",sep = "_"), quote = FALSE)

################################
# Manual marker gene analysis
################################
pdf(paste0(Sys.Date(), "_10x_ovary_VCD_Ptprc_negative_cells_Vlnplot_marker_genes_clusters_0_to_7.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.VCD.Ptprc.neg, subset = seurat_clusters %in% c(0:7)), 
        features = c("Ddx4", "Amh", "Foxl2", "Cyp19a1",
                     "Cyp17a1", "Cyp11a1", "Lhcgr", 
                     "Pdgfra", "Col1a1", 
                     "Rgs5", "Acta2", "Notch3", 
                     "Vwf", "Flt1", "Flt4", "Prox1", 
                     "Pax8", "Cdh1"))
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_VCD_Ptprc_negative_cells_Vlnplot_marker_genes_clusters_8_to_16.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.VCD.Ptprc.neg, subset = seurat_clusters %in% c(8:16)), 
        features = c("Ddx4", "Amh", "Foxl2", "Cyp19a1",
                     "Cyp17a1", "Cyp11a1", "Lhcgr", 
                     "Pdgfra", "Col1a1", 
                     "Rgs5", "Acta2", "Notch3", 
                     "Vwf", "Flt1", "Flt4", "Prox1", 
                     "Pax8", "Cdh1"))
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_VCD_Ptprc_negative_cells_Vlnplot_marker_genes_clusters_17_to24.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.VCD.Ptprc.neg, subset = seurat_clusters %in% c(17:24)), 
        features = c("Ddx4", "Amh", "Foxl2", "Cyp19a1",
                     "Cyp17a1", "Cyp11a1", "Lhcgr", 
                     "Pdgfra", "Col1a1", 
                     "Rgs5", "Acta2", "Notch3", 
                     "Vwf", "Flt1", "Flt4", "Prox1", 
                     "Pax8", "Cdh1"))
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_VCD_Ptprc_negative_cells_Vlnplot_marker_genes_clusters_24_to_32.pdf"), width = 15, height = 20)
VlnPlot(subset(ovary.VCD.Ptprc.neg, subset = seurat_clusters %in% c(24:32)), 
        features = c("Ddx4", "Amh", "Foxl2", "Cyp19a1",
                     "Cyp17a1", "Cyp11a1", "Lhcgr", 
                     "Pdgfra", "Col1a1", 
                     "Rgs5", "Acta2", "Notch3", 
                     "Vwf", "Flt1", "Flt4", "Prox1", 
                     "Pax8", "Cdh1"))
dev.off()

################################
# Transfer final annotations
################################

ovary.VCD.Ptprc.neg <- SetIdent(ovary.VCD.Ptprc.neg, value = "seurat_clusters")

VlnPlot(ovary.VCD.Ptprc.neg, features = c("nCount_RNA"))
VlnPlot(ovary.VCD.Ptprc.neg, features = c("nFeature_RNA"))
VlnPlot(ovary.VCD.Ptprc.neg, features = c("percent.mito"))
VlnPlot(ovary.VCD.Ptprc.neg, features = c("decontX_contamination"))

ovary.VCD.Ptprc.neg@meta.data$celltype.level2 <- rep("NA",dim(ovary.VCD.Ptprc.neg@meta.data)[1])
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 0] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 1] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 2] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 3] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 4] <- "LEC"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 5] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 6] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 7] <- "LEC"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 8] <- "BEC"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 9] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 10] <- "Theca"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 11] <- "BEC"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 12] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 13] <- "Epithelial"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 14] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 15] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 16] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 17] <- "Theca"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 18] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 19] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 20] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 21] <- "Epithelial"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 22] <- "Theca"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 23] <- "BEC"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 24] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 25] <- "BEC"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 26] <- "BEC"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 27] <- "Theca"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 28] <- "Epithelial"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 29] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 30] <- "Granulosa"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 31] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 32] <- "Epithelial"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 33] <- "Mesenchyme"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 34] <- "Granulosa"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 35] <- "Granulosa"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 36] <- "Granulosa"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 37] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 38] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 39] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 40] <- "Epithelial"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 41] <- "Pericyte"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 42] <- "Granulosa"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 43] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 44] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 45] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 46] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 47] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 48] <- "BEC"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 49] <- "BEC"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 50] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 51] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 52] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 53] <- "BEC"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 54] <- "BEC"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 55] <- "Granulosa"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 56] <- "LEC"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 57] <- "Stroma"
ovary.VCD.Ptprc.neg@meta.data$celltype.level2[ovary.VCD.Ptprc.neg@meta.data$seurat_clusters %in% 58] <- "Epithelial"

save(ovary.VCD.Ptprc.neg, file = paste(Sys.Date(),"10x_ovary_Benayoun_lab_VCD_Ptprc_negative_cells_Seurat_object_with_final_annotation.RData",sep = "_"))

################################################################################
# 4.Combine cell type annotation
################################################################################

# Extract celltype.level2 metadata
celltype.level2.pos <- ovary.VCD.Ptprc.pos@meta.data[, "celltype.level2", drop = FALSE]
celltype.level2.neg <- ovary.VCD.Ptprc.neg@meta.data[, "celltype.level2", drop = FALSE]

# Combine the metadata from both subsetted objects
combined.celltype.level2 <- rbind(celltype.level2.pos, celltype.level2.neg)

# Ensure that the cell IDs match
common_cell_ids <- intersect(rownames(ovary.VCD@meta.data), rownames(combined.celltype.level2))

# Update the celltype.level2 in ovary.VCD for the matching cells
ovary.VCD@meta.data[common_cell_ids, "celltype.level2"] <- combined.celltype.level2[common_cell_ids, "celltype.level2"]

# Verify that the metadata has been successfully updated
table(is.na(ovary.VCD@meta.data$celltype.level2))
# FALSE 
# 47386

save(ovary.VCD, file = paste(Sys.Date(),"10x_ovary_Benayoun_lab_VCD_Seurat_object_with_final_annotation.RData",sep = "_"))

################################################################################
# 5.Generate plots
################################################################################

ovary.VCD <- SetIdent(ovary.VCD, value = "celltype.level2")
levels(ovary.VCD) <- c("Granulosa",
                       "Theca",
                       "Stroma",
                       "Mesenchyme",
                       "Pericyte",
                       "BEC",
                       "LEC",
                       "Epithelial",
                       "Neutrophil", "Monocyte", "Macrophage" , "DC", "ILC", "NK", "NKT", "CD8 NKT", "CD4 T", "DNT", "DPT", "B")

pdf(paste(Sys.Date(),"10X_ovary_Benayoun_lab_VCD_UMAP_with_label.pdf",sep = "_"), height = 7, width = 10)
DimPlot(ovary.VCD, reduction = "umap", shuffle = TRUE, label = TRUE)
dev.off()

pdf(paste(Sys.Date(),"10X_ovary_Benayoun_lab_VCD_UMAP_without_label.pdf",sep = "_"), height = 7, width = 10)
DimPlot(ovary.VCD, reduction = "umap", shuffle = TRUE, label = FALSE)
dev.off()

pdf(paste(Sys.Date(),"10X_ovary_Benayoun_lab_VCD_Dotplot_marker_genes.pdf",sep = "_"), height = 10, width = 15)
DotPlot(ovary.VCD, features = c("Amh", "Cyp19a1",           # Granulosa cell marker
                                "Cyp17a1", "Cyp11a1",       # Theca cell marker 
                                "Pdgfra", "Col1a1",         # Stroma cell marker
                                "Acta2", "Notch3", "Rgs5",  # Mesenchymal cell marker (pericyte + smooth muscle cell)
                                "Flt1",  "Prox1",           # Endothelial cell marker (BEC + LEC)
                                "Upk1b",                    # Epithelial cell marker
                                "S100a8", "Adgre1", "C1qa","Flt3", "Il2ra", "Klrb1c", "Gzma", "Nkg7", "Cd3e", "Cd8b1", "Cd4", "Cd19"))                   
dev.off()

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_VCD_cell_type_annotation_session_Info.txt", sep =""))
sessionInfo()
sink()

