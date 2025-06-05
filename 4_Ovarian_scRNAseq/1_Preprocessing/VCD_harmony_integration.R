setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/1_Preprocessing")
options(stringsAsFactors = FALSE)

library(harmony)
library(Seurat)
library(SeuratWrappers)

library(ggplot2)

# rm(list = ls())

################################################################################
# VCD data
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

load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/1_Preprocessing/Input/2024-10-21_10x_ovary_Benayoun_lab_VCD_Seurat_object_SINGLETS.RData")
ovary.VCD.DFsinglets
# An object of class Seurat 
# 48969 features across 47357 samples within 2 assays 
# Active assay: SCT (24484 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

DefaultAssay(ovary.VCD.DFsinglets) <- "RNA"

# SCTransform should be performed separately
# https://satijalab.org/seurat/archive/v3.1/integration.html#sctransform
# https://github.com/immunogenomics/harmony/issues/41

ovary.VCD.list <- SplitObject(ovary.VCD.DFsinglets, split.by = "Batch")

ovary.VCD.list <- lapply(X = ovary.VCD.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- SCTransform(x, vars.to.regress = c("nFeature_RNA", "percent.mito", "decontX_contamination", "Library"))
})

ovary.features <- SelectIntegrationFeatures(object.list = ovary.VCD.list, verbose = FALSE, selection.method = "vst", nfeatures = 2000)

ovary.VCD.merged <- merge(ovary.VCD.list[[1]], 
                          y = ovary.VCD.list[2:length(ovary.VCD.list)],  
                          project = "ovary.VCD", 
                          merge.data = TRUE)

VariableFeatures(ovary.VCD.merged) <- ovary.features
ovary.VCD.harmony <- ScaleData(object = ovary.VCD.merged)
ovary.VCD.harmony <- RunPCA(object = ovary.VCD.harmony)

# Datasets integration via Harmony

# ovary.VCD.harmony

# Get cell embeddings and metadata from Seurat object
ovary.VCD.harmony.pc <- Embeddings(ovary.VCD.harmony, reduction = 'pca')
ovary.VCD.harmony.metadata <- ovary.VCD.harmony@meta.data

p1 <- do_scatter(ovary.VCD.harmony.pc, ovary.VCD.harmony.metadata, 'Library') + 
  labs(title = 'Colored by Library')

p2 <- do_scatter(ovary.VCD.harmony.pc, ovary.VCD.harmony.metadata, 'Batch') + 
  labs(title = 'Colored by Batch')

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_VCD_Harmony_integration_Scatterplots_PCA_before_integration.pdf"), width = 15, height = 5)
cowplot::plot_grid(p1, p2)
dev.off()

harmony_embeddings <- harmony::HarmonyMatrix(ovary.VCD.harmony.pc, 
                                             ovary.VCD.harmony.metadata, 
                                             c("Library", "Batch"), 
                                             do_pca = FALSE, verbose=FALSE)

p1 <- do_scatter(harmony_embeddings, ovary.VCD.harmony.metadata, 'Library') + 
  labs(title = 'Colored by Library')

p2 <- do_scatter(harmony_embeddings, ovary.VCD.harmony.metadata, 'Batch') + 
  labs(title = 'Colored by Batch')

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_VCD_Harmony_integration_Scatterplots_PCA_after_integration.pdf"), width = 15, height = 5)
cowplot::plot_grid(p1, p2)
dev.off()

save(harmony_embeddings, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_harmony_embeddings.RData"))

ovary.VCD.harmony[["harmony"]] <- CreateDimReducObject(embeddings = harmony_embeddings, key = "Harmony_", assay = DefaultAssay(ovary.VCD.harmony))
ovary.VCD.harmony <- RunUMAP(ovary.VCD.harmony, reduction = "harmony", dims = 1:30)

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_VCD_Harmony_integration_UMAP_Dimplot_splitby_batch.pdf"), width = 15, height = 5)
DimPlot(ovary.VCD.harmony, reduction = "umap", split.by = "Batch", ncol = 3, label = TRUE, shuffle = TRUE)
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_VCD_Harmony_integration_UMAP_Dimplot_groupby_batch.pdf"), width = 8, height = 5)
DimPlot(ovary.VCD.harmony, reduction = "umap", group.by = "Batch", shuffle = TRUE)
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_VCD_Harmony_integration_UMAP_Dimplot_groupby_library.pdf"), width = 8, height = 5)
DimPlot(ovary.VCD.harmony, reduction = "umap", group.by = "Library", shuffle = TRUE)
dev.off()

table(ovary.VCD.harmony@meta.data$seurat_clusters, ovary.VCD.harmony$Batch)

save(ovary.VCD.harmony, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_integrated_harmony.RData"))

########################################################################################################################################################
sink(file = paste0(Sys.Date(),"_VCD_data_integration_using_Harmony_session_Info.txt"))
sessionInfo()
sink()
