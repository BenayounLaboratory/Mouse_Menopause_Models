setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/1_Preprocessing")
options(stringsAsFactors = FALSE)

library(harmony)
library(Seurat)
library(SeuratWrappers)

library(ggplot2)

# rm(list = ls())

################################################################################
# Foxl2 haploinsufficiency data
# Data integration to account for batch effects - Harmony
# https://github.com/satijalab/seurat-wrappers/blob/master/docs/harmony.md
# https://hbctraining.github.io/scRNA-seq_online/lessons/06a_integration_harmony.html
# https://portals.broadinstitute.org/harmony/articles/quickstart.htmls
# Remove batch-specific clusters
################################################################################

##################
do_scatter <- function(umap_use, meta_data, label_name, no_guides = TRUE,
                       do_labels = TRUE, nice_names,
                       pt_size = 4, point_size = .5, base_size = 12, 
                       do_points = TRUE, do_density = FALSE, h = 6, w = 8) {
  umap_use <- umap_use[, 1:2]
  colnames(umap_use) <- c('X1', 'X2')
  plt_df <- umap_use %>% data.frame() %>% 
    cbind(meta_data) %>% 
    dplyr::sample_frac(1L) 
  plt_df$given_name <- plt_df[[label_name]]
  
  if (!missing(nice_names)) {
    plt_df %<>%
      dplyr::inner_join(nice_names, by = "given_name") %>% 
      subset(nice_name != "" & !is.na(nice_name))
    
    plt_df[[label_name]] <- plt_df$nice_name        
  }
  
  plt <- plt_df %>% 
    ggplot2::ggplot(aes_string("X1", "X2", col = label_name, fill = label_name)) + 
    theme_test(base_size = base_size) + 
    theme(panel.background = element_rect(fill = NA, color = "black")) + 
    guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1,
                                                    shape = 16, size = 4)), 
           alpha = FALSE) +
    theme(plot.title = element_text(hjust = .5)) + 
    labs(x = "PC 1", y = "PC 2") 
  
  if (do_points) 
    plt <- plt + geom_point(shape = '.')
  if (do_density) 
    plt <- plt + geom_density_2d()    
  
  
  if (no_guides)
    plt <- plt + guides(col = FALSE, fill = FALSE, alpha = FALSE)
  
  if (do_labels) {
    data_labels <- plt_df %>% 
      dplyr::group_by_(label_name) %>% 
      dplyr::summarise(X1 = mean(X1), X2 = mean(X2)) %>% 
      dplyr::ungroup()
    
    plt <- plt + geom_label(data = data_labels, label.size = NA,
                            aes_string(label = label_name), 
                            color = "white", size = pt_size, alpha = 1,
                            segment.size = 0) +
      guides(col = FALSE, fill = FALSE)
  }
  
  return(plt)
}

##################

# We library normalized the cells, 
# log transformed the counts, and 
# scaled the genes. 
# Then we performed PCA and kept the top PCs.

# load dataset

load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/1_Preprocessing/Input/2024-10-24_10x_ovary_Benayoun_lab_Foxl2_Seurat_object_SINGLETS.RData")
ovary.Foxl2.DFsinglets
# An object of class Seurat 
# 44755 features across 21301 samples within 2 assays 
# Active assay: SCT (22376 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

DefaultAssay(ovary.Foxl2.DFsinglets) <- "RNA"

# SCTransform should be performed separately
# https://satijalab.org/seurat/archive/v3.1/integration.html#sctransform
# https://github.com/immunogenomics/harmony/issues/41

ovary.Foxl2.list <- SplitObject(ovary.Foxl2.DFsinglets, split.by = "Batch")

ovary.Foxl2.list <- lapply(X = ovary.Foxl2.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- SCTransform(x, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))
})

ovary.features <- SelectIntegrationFeatures(object.list = ovary.Foxl2.list, verbose = FALSE, selection.method = "vst", nfeatures = 2000)

ovary.Foxl2.merged <- merge(ovary.Foxl2.list[[1]], 
                          y = ovary.Foxl2.list[2:length(ovary.Foxl2.list)],  
                          project = "ovary.Foxl2", 
                          merge.data = TRUE)

VariableFeatures(ovary.Foxl2.merged) <- ovary.features
ovary.Foxl2.harmony <- ScaleData(object = ovary.Foxl2.merged)
ovary.Foxl2.harmony <- RunPCA(object = ovary.Foxl2.harmony)

# Datasets integration via Harmony

ovary.Foxl2.harmony
# An object of class Seurat 
# 44731 features across 21301 samples within 2 assays 
# Active assay: SCT (22352 features, 2000 variable features)
# 1 other assay present: RNA
# 1 dimensional reduction calculated: pca

# Get cell embeddings and metadata from Seurat object
ovary.Foxl2.harmony.pc <- Embeddings(ovary.Foxl2.harmony, reduction = 'pca')
ovary.Foxl2.harmony.metadata <- ovary.Foxl2.harmony@meta.data

p1 <- do_scatter(ovary.Foxl2.harmony.pc, ovary.Foxl2.harmony.metadata, 'Library') + 
  labs(title = 'Colored by Library')

p2 <- do_scatter(ovary.Foxl2.harmony.pc, ovary.Foxl2.harmony.metadata, 'Batch') + 
  labs(title = 'Colored by Batch')

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_Foxl2_Harmony_integration_Scatterplots_PCA_before_integration.pdf"), width = 15, height = 5)
cowplot::plot_grid(p1, p2)
dev.off()

harmony_embeddings <- harmony::HarmonyMatrix(ovary.Foxl2.harmony.pc, 
                                             ovary.Foxl2.harmony.metadata, 
                                             c("Library", "Batch"), 
                                             do_pca = FALSE, verbose=FALSE)

p1 <- do_scatter(harmony_embeddings, ovary.Foxl2.harmony.metadata, 'Library') + 
  labs(title = 'Colored by Library')

p2 <- do_scatter(harmony_embeddings, ovary.Foxl2.harmony.metadata, 'Batch') + 
  labs(title = 'Colored by Batch')

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_Foxl2_Harmony_integration_Scatterplots_PCA_after_integration.pdf"), width = 15, height = 5)
cowplot::plot_grid(p1, p2)
dev.off()

save(harmony_embeddings, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_harmony_embeddings.RData"))

ovary.Foxl2.harmony[["harmony"]] <- CreateDimReducObject(embeddings = harmony_embeddings, key = "Harmony_", assay = DefaultAssay(ovary.Foxl2.harmony))
ovary.Foxl2.harmony <- RunUMAP(ovary.Foxl2.harmony, reduction = "harmony", dims = 1:30)

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_Foxl2_Harmony_integration_UMAP_Dimplot_splitby_batch.pdf"), width = 15, height = 5)
DimPlot(ovary.Foxl2.harmony, reduction = "umap", split.by = "Batch", ncol = 3, label = TRUE, shuffle = TRUE)
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_Foxl2_Harmony_integration_UMAP_Dimplot_groupby_batch.pdf"), width = 8, height = 5)
DimPlot(ovary.Foxl2.harmony, reduction = "umap", group.by = "Batch", shuffle = TRUE)
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_Foxl2_Harmony_integration_UMAP_Dimplot_groupby_library.pdf"), width = 8, height = 5)
DimPlot(ovary.Foxl2.harmony, reduction = "umap", group.by = "Library", shuffle = TRUE)
dev.off()

table(ovary.Foxl2.harmony@meta.data$seurat_clusters, ovary.Foxl2.harmony$Batch)

save(ovary.Foxl2.harmony, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_integrated_harmony.RData"))

########################################################################################################################################################
sink(file = paste0(Sys.Date(),"_Foxl2_data_integration_using_Harmony_session_Info.txt"))
sessionInfo()
sink()
