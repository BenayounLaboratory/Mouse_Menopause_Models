library(Seurat)
library(future)

rm(list = ls())

plan(sequential)
options(future.globals.maxSize = 80 * 1024^3)  # 80 GB

#################################
# Menopause-model project / revision
# Integrate mouse ovarian scRNA-seq datasets (AC, VCD, Foxl2) via Seurat
# SCTransform + RPCA anchoring
#################################

#################################
# 1. Load data
#################################

load("/project2/bbenayou_34/kim/Ovarian_aging_single_cell_project/Celltype_annotated_Seurat_objects/2025-11-05_10x_ovary_Benayoun_lab_AC_Seurat_object_with_final_annotation.RData")
load("/project2/bbenayou_34/kim/Ovarian_aging_single_cell_project/Celltype_annotated_Seurat_objects/2026-04-28_10x_ovary_Foxl2_Seurat_object_celltypes_annotated.RData")
load("/project2/bbenayou_34/kim/Ovarian_aging_single_cell_project/Celltype_annotated_Seurat_objects/2025-06-27_10x_ovary_Benayoun_lab_VCD_Seurat_object_with_final_annotation.RData")

#################################
# 2. Update old Seurat objects
#################################

ovary.AC    <- UpdateSeuratObject(ovary.AC)
ovary.VCD   <- UpdateSeuratObject(ovary.VCD)
ovary.Foxl2 <- UpdateSeuratObject(ovary.Foxl2)

#################################
# 3. Label each object by dataset
#################################

ovary.AC$Dataset    <- "AC"
ovary.VCD$Dataset   <- "VCD"
ovary.Foxl2$Dataset <- "Foxl2"

obj.list <- list(
  AC = ovary.AC,
  VCD = ovary.VCD,
  Foxl2 = ovary.Foxl2
)

#################################
# 4. SCTransform + PCA per dataset
#################################

obj.list <- lapply(obj.list, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(
    x,
    method = "glmGamPoi",
    return.only.var.genes = TRUE,
    conserve.memory = TRUE,
    verbose = FALSE
  )
  x <- RunPCA(x, npcs = 50, verbose = FALSE)
  x
})

save(obj.list, file = paste0(Sys.Date(), "_MeMo_object_list_post_SCTransform.RData"))

#################################
# 5. Select integration features + prep SCT integration
#################################

features <- SelectIntegrationFeatures(
  object.list = obj.list,
  nfeatures = 2000
)

obj.list <- PrepSCTIntegration(
  object.list = obj.list,
  anchor.features = features,
  verbose = FALSE
)

save(obj.list, file = paste0(Sys.Date(), "_MeMo_object_list_post_prepsctintegration.RData"))

#################################
# 6. Run PCA using the same integration features
#################################

obj.list <- lapply(obj.list, function(x) {
  x <- RunPCA(
    x,
    features = features,
    npcs = 50,
    verbose = FALSE
  )
  return(x)
})

save(obj.list, file = paste0(Sys.Date(), "_MeMo_object_list_post_PCA.RData"))

#################################
# 7. Find RPCA anchors + integrate data
#################################

anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "rpca",
  reference = c(1, 2, 3),  # AC, VCD, Foxl2
  dims = 1:30,
  k.anchor = 5,
  verbose = FALSE
)

save(anchors, file = paste0(Sys.Date(), "_MeMo_anchors_by_reference.RData"))

ovary.integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  dims = 1:30,
  k.weight = 50,
  preserve.order = TRUE,
  verbose = TRUE
)

#################################
# 8. Downstream analysis on the integrated object
#################################

DefaultAssay(ovary.integrated) <- "integrated"

ovary.integrated <- RunPCA(ovary.integrated, npcs = 50, verbose = FALSE)
ovary.integrated <- RunUMAP(ovary.integrated, dims = 1:30)
ovary.integrated <- FindNeighbors(ovary.integrated, dims = 1:30)
ovary.integrated <- FindClusters(ovary.integrated, resolution = 0.5)

save(ovary.integrated, file = paste0(Sys.Date(), "_MeMo_RPCA_integrated_object.RData"))

#################################
sink(file = paste0(Sys.Date(), "_Cross_model_integration_session_info.txt"))
sessionInfo()
sink()
