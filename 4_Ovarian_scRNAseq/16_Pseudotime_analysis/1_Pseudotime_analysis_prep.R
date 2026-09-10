library(Seurat)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)
library(tidyverse)

options(future.globals.maxSize = 8000 * 1024^2)

###############################################################################
# Menopause-model project / revision
# Prepare granulosa and theca CDS objects for interactive pseudotime root selection
###############################################################################

###############################################################################
# 1. Load integrated object, assign numeric age from Library name, subset to granulosa / theca
###############################################################################

load("./MeMo_RPCA_integrated_object.RData")

# Numeric age per cell, from its Library name
assign_age_numeric <- function(library_names) {
  age <- rep(NA_real_, length(library_names))

  age[grepl("^YF", library_names)] <- 4
  age[grepl("^OF", library_names)] <- 20

  age[grepl("3m_30d",  library_names)] <- 4
  age[grepl("3m_90d",  library_names)] <- 6
  age[grepl("10m_30d", library_names)] <- 11
  age[grepl("10m_90d", library_names)] <- 13

  age[grepl("Foxl2_wt_young",   library_names)] <- 3
  age[grepl("Foxl2_het_young",  library_names)] <- 3
  age[grepl("Foxl2_wt_supold",  library_names)] <- 17
  age[grepl("Foxl2_het_supold", library_names)] <- 17
  age[grepl("Foxl2_wt_old",     library_names)] <- 9
  age[grepl("Foxl2_het_old",    library_names)] <- 9

  age
}

ovary.integrated$Age_numeric <- assign_age_numeric(ovary.integrated$Library)

# Library -> Age_numeric mapping check 
lib_age_map <- unique(ovary.integrated@meta.data[, c("Library", "Age_numeric")])
lib_age_map <- lib_age_map[order(lib_age_map$Age_numeric, lib_age_map$Library), ]
print(lib_age_map, row.names = FALSE)

n_missing_age <- sum(is.na(ovary.integrated$Age_numeric))

granulosa.integrated <- subset(ovary.integrated, subset = celltype.level2 == "Granulosa")
theca.integrated     <- subset(ovary.integrated, subset = celltype.level2 == "Theca")

###############################################################################
# 2. Per cell type: RPCA re-integration, monocle3 CDS, trajectory graph, diagnostic plots
###############################################################################

cell_type_list <- list(
  list(obj = granulosa.integrated, label = "granulosa"),
  list(obj = theca.integrated,     label = "theca")
)

for (ct in cell_type_list) {

  obj   <- ct$obj
  label <- ct$label

  obj.list <- SplitObject(obj, split.by = "Batch")

  lib_sizes <- sapply(obj.list, ncol)

  min_cells <- 25
  dropped <- names(lib_sizes)[lib_sizes < min_cells]
  if (length(dropped) > 0) {
    cat("\nExcluding", length(dropped), "libraries with <", min_cells,
        "cells for", label, ":", paste(dropped, collapse = ", "), "\n")
    obj.list <- obj.list[!(names(obj.list) %in% dropped)]
    lib_sizes <- lib_sizes[!(names(lib_sizes) %in% dropped)]
  }

  safe_npcs <- min(30, min(lib_sizes) - 1)

  obj.list <- lapply(obj.list, function(x) {
    SCTransform(x, vst.flavor = "v2", verbose = FALSE)
  })

  max_npcs_per_lib <- sapply(obj.list, function(x) {
    min(nrow(x[["SCT"]]$scale.data), ncol(x)) - 1
  })
  safe_npcs <- min(safe_npcs, max_npcs_per_lib)

  obj.list <- lapply(obj.list, function(x) {
    RunPCA(x, npcs = safe_npcs, verbose = FALSE)
  })

  features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
  obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)

  anchors <- FindIntegrationAnchors(
    object.list     = obj.list,
    anchor.features = features,
    normalization.method = "SCT",
    reduction       = "rpca",
    dims            = 1:safe_npcs
  )

  obj.integrated <- IntegrateData(anchorset = anchors,
                                  normalization.method = "SCT",
                                  dims = 1:safe_npcs,
                                  k.weight = 20)

  DefaultAssay(obj.integrated) <- "integrated"
  obj.integrated <- RunPCA(obj.integrated, npcs = safe_npcs, verbose = FALSE)
  obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:safe_npcs)
  obj.integrated <- FindNeighbors(obj.integrated, reduction = "pca", dims = 1:safe_npcs)
  obj.integrated <- FindClusters(obj.integrated, resolution = 0.5)

  saveRDS(obj.integrated, paste0(Sys.Date(), "_", label, "_reintegrated.rds"))

  # Build monocle3 CDS, transfer UMAP, learn graph

  seurat_obj <- obj.integrated

  if ("RNA" %in% names(seurat_obj@assays) &&
      inherits(seurat_obj[["RNA"]], "Assay5")) {
    seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
  }

  counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  cell_meta  <- seurat_obj@meta.data
  gene_meta  <- data.frame(gene_short_name = rownames(counts_mat),
                           row.names       = rownames(counts_mat))

  cds <- new_cell_data_set(
    expression_data = counts_mat,
    cell_metadata   = cell_meta,
    gene_metadata   = gene_meta
  )

  cds <- estimate_size_factors(cds)

  umap_embed <- Embeddings(seurat_obj, "umap")
  reducedDims(cds)[["UMAP"]] <- umap_embed

  set.seed(1234)
  cds <- cluster_cells(cds, reduction_method = "UMAP")

  cluster_sizes <- table(clusters(cds))

  cds <- learn_graph(cds, use_partition = FALSE)

  saveRDS(cds, paste0(Sys.Date(), "_", label, "_cds_prerooted.rds"))

  # Diagnostic plots colored by Age / Age_numeric / Dataset / Library / Treatment / Genotype

  color_vars <- c("Age", "Age_numeric", "Dataset", "Library", "Treatment", "Genotype")
  color_vars <- color_vars[color_vars %in% colnames(colData(cds))]

  pdf_path <- paste0(Sys.Date(), "_", label, "_diagnostic_plots.pdf")
  pdf(pdf_path, width = 10, height = 8)

  for (v in color_vars) {
    p <- plot_cells(cds,
                    color_cells_by          = v,
                    label_groups_by_cluster = FALSE,
                    label_leaves            = FALSE,
                    label_branch_points     = FALSE,
                    trajectory_graph_color  = "grey40") +
      ggtitle(paste(label, "-", v))
    print(p)
  }

  if ("Library" %in% color_vars) {
    p_lib <- plot_cells(cds,
                        color_cells_by          = "Library",
                        label_groups_by_cluster = TRUE,
                        label_leaves            = FALSE,
                        label_branch_points     = FALSE,
                        trajectory_graph_color  = "grey40") +
      ggtitle(paste(label, "- Library (for root selection reference)"))
    print(p_lib)
  }

  dev.off()

  rm(obj.list, anchors, obj.integrated, seurat_obj, counts_mat, cell_meta, gene_meta, umap_embed, cds, p)
  gc()
}

###############################################################################
sink(file = paste0(Sys.Date(), "_1_Pseudotime_analysis_prep_session_info.txt"))
sessionInfo()
sink()
