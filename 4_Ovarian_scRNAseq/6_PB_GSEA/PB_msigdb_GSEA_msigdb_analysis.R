setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/6_PB_GSEA")
options(stringsAsFactors = FALSE)

# Load libraries
library('clusterProfiler')
library('msigdbr')
library('ggplot2')
library('scales')
library(stringr)
library('reshape2')
library(dplyr)

theme_set(theme_bw())
rm(list = ls())

################################################################################
# 1. Load MSigDB Gene Sets
################################################################################

# Retrieve MSigDB gene sets for mouse
msig_bp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
msig_reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")

# Combine gene sets into a list
all_sets <- list(
  BP = msig_bp,
  Reactome = msig_reactome
)

################################################################################
# 2. Define Functions
################################################################################

perform_gsea <- function(geneList, gene_sets, output_prefix, fdr_cutoff = 0.1) {
  GSEA(
    geneList,
    TERM2GENE = gene_sets[, c("gs_name", "gene_symbol")],
    pvalueCutoff = fdr_cutoff
  )
}

run_gsea_for_all_celltypes <- function(deseq_results_list, all_sets, fdr_cutoff = 0.1) {
  gsea_results_list <- list()
  
  for (cell_type in names(deseq_results_list)) {
    gsea_results_list[[cell_type]] <- list()
    
    for (db_name in names(all_sets)) {
      geneList <- sort(deseq_results_list[[cell_type]]$stat, decreasing = TRUE)
      names(geneList) <- rownames(deseq_results_list[[cell_type]])
      
      gsea_result <- tryCatch(
        perform_gsea(geneList, all_sets[[db_name]], paste(cell_type, db_name), fdr_cutoff),
        error = function(e) NULL
      )
      
      if (!is.null(gsea_result) && nrow(gsea_result@result) > 0) {
        gsea_results_list[[cell_type]][[db_name]] <- gsea_result@result
      } else {
        message("No significant results for: ", cell_type, " using ", db_name)
      }
    }
  }
  
  return(gsea_results_list)
}

extract_terms_in_multiple_celltypes <- function(gsea_results_list, fdr_cutoff = 0.1, min_celltypes = 4) {
  terms_across_celltypes <- list()
  
  for (db_name in unique(unlist(lapply(gsea_results_list, names)))) {
    terms_per_celltype <- lapply(gsea_results_list, function(celltype_results) {
      if (db_name %in% names(celltype_results)) {
        sig_terms <- celltype_results[[db_name]]
        sig_terms <- sig_terms[sig_terms$qvalues < fdr_cutoff, "Description"]
        return(sig_terms)
      } else {
        return(character(0))
      }
    })
    
    all_terms <- unlist(terms_per_celltype)
    term_counts <- table(all_terms)
    common_terms <- names(term_counts[term_counts >= min_celltypes])
    terms_across_celltypes[[db_name]] <- common_terms
  }
  
  return(terms_across_celltypes)
}

################################################################################
# 3. Define model paths and run
################################################################################

model_files <- list(
  Aging = "./Input/2025-01-14_10x_ovary_Benayoun_lab_AC_PB_DESeq2_object.RData",
  VCD = "./Input/2025-01-15_10x_ovary_Benayoun_lab_VCD_PB_DESeq2_object.RData",
  Foxl2 = "./Input/2025-01-15_10x_ovary_Benayoun_lab_Foxl2_NULL_PB_DESeq2_object.RData"
)

gsea_all_models <- list()
common_terms_all_models <- list()

for (model in names(model_files)) {

  # Load data
  load(model_files[[model]]) 
  
  # Run GSEA
  gsea.results.list <- run_gsea_for_all_celltypes(deseq.res.list, all_sets, fdr_cutoff = 0.1)
  gsea_all_models[[model]] <- gsea.results.list
  
  # Extract shared terms
  terms_in_multiple_celltypes <- extract_terms_in_multiple_celltypes(gsea.results.list)
  common_terms_all_models[[model]] <- terms_in_multiple_celltypes
  
  # Save results
  save(gsea.results.list, file = paste0(Sys.Date(), "_", model, "_PB_GSEA_list_object.RData"))
  save(terms_in_multiple_celltypes, file = paste0(Sys.Date(), "_", model, "_PB_GSEA_terms_found_in_multiple_celltypes_list.RData"))
}

################################################################################
# 4. Generate combined bubble plot
################################################################################

# Define functions
combine_gsea_results_all <- function(gsea_list, dataset_name, db_name) {
  res <- lapply(names(gsea_list), function(cell_type) {
    if (!is.null(gsea_list[[cell_type]][[db_name]])) {
      df <- gsea_list[[cell_type]][[db_name]]
      df$CellType <- cell_type
      df$Dataset <- dataset_name
      df$DB <- db_name
      return(df)
    }
  })
  do.call(rbind, res)
}

plot_bubble_shared_terms <- function(df, term_list, db_name, file_prefix = Sys.Date()) {
  bubble_data <- df %>%
    filter(DB == db_name, Description %in% term_list) %>%
    mutate(
      Pathway = str_to_title(str_remove(Description, "GOBP_")),
      MinusLog10Pval = -log10(p.adjust),  # Using adjusted p-value
      Group = interaction(CellType, Dataset, drop = TRUE)
    ) %>%
    select(Pathway, CellType, Dataset, Group, NES, MinusLog10Pval)
  
  # Define color palette
  my_colors <- c("darkblue", "dodgerblue4", "dodgerblue3", "dodgerblue1", 
                 "white", "lightcoral", "brown1", "firebrick2", "firebrick4")
  my_values <- rescale(c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3))
  
  # Plot
  pdf(paste0(file_prefix, "_Shared_GSEA_", db_name, "_BubblePlot.pdf"), height = 6, width = 15)
  print(
    ggplot(bubble_data, aes(x = Group, y = Pathway, size = MinusLog10Pval, color = NES)) +
      geom_point(alpha = 0.8) +
      scale_color_gradientn(colors = my_colors, values = my_values, limits = c(-3, 3), na.value = "grey70") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = paste("Shared GSEA Terms Across Models -", db_name),
        x = "Cell Type - Dataset", y = "Pathway",
        color = "NES", size = "-log10(FD)"
      )
  )
  dev.off()
}

# Combine results across all models and both databases
all_gsea_results <- bind_rows(
  combine_gsea_results_all(gsea_all_models$Aging, "Aging", "BP"),
  combine_gsea_results_all(gsea_all_models$VCD, "VCD", "BP"),
  combine_gsea_results_all(gsea_all_models$Foxl2, "Foxl2", "BP"),
  combine_gsea_results_all(gsea_all_models$Aging, "Aging", "Reactome"),
  combine_gsea_results_all(gsea_all_models$VCD, "VCD", "Reactome"),
  combine_gsea_results_all(gsea_all_models$Foxl2, "Foxl2", "Reactome")
)

# Extract shared terms per DB (more than 4 cell types in at least 1 model)

get_terms_from_list <- function(overlap_list, db = "BP") {
  if (!is.null(overlap_list[[db]])) return(unique(overlap_list[[db]]))
  return(character(0))
}

combined_terms <- list(
  BP = unique(unlist(lapply(common_terms_all_models, get_terms_from_list, db = "BP"))),
  Reactome = unique(unlist(lapply(common_terms_all_models, get_terms_from_list, db = "Reactome")))
)

# Generate Plots for BP and Reactome

plot_bubble_shared_terms(all_gsea_results, combined_terms$BP, "BP")
plot_bubble_shared_terms(all_gsea_results, combined_terms$Reactome, "Reactome")

################################################################################
sink(file = paste0(Sys.Date(), "_PB_msigdb_GSEA_session_Info.txt"))
sessionInfo()
sink()
