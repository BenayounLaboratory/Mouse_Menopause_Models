setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/3_Proportion_test")
options(stringsAsFactors = F)

# Load required libraries
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scProportionTest)
library(gridExtra)
library(beeswarm)

rm(list=ls())

################################################################################
# Perform cell type proportion tests
# scProportionTest
# Foxl2 dataset
################################################################################

################################################################################
# Define color scheme and custom order of celltypes
################################################################################

# Define colors for different ages

# Define the endpoint colors
colors_wt <- colorRampPalette(c("deeppink1", "deeppink4"))
colors_het <- colorRampPalette(c("#ADDFE4", "#508B8E"))

# Generate intermediate colors
month_wt_colors <- colors_wt(10)

# Assigning the colors to specific ages
wt_colors <- list(
  "3_months"  = month_wt_colors[1],
  "4_months"  = "deeppink1",      
  "5_months"  = month_wt_colors[3],
  "6_months"  = month_wt_colors[4],
  "7_months"  = month_wt_colors[5],
  "8_months"  = month_wt_colors[6],
  "10_months" = month_wt_colors[7],
  "12_months" = month_wt_colors[8],
  "14_months" = month_wt_colors[9],
  "20_months" = "deeppink4"
)

# [WT] Young: deeppink1, midage: #B10D66
# [HET] Young:#ADDFE4, midage: #508B8E
# Celltype order

custom_order_level1 <- c("Ptprc.neg", "Ptprc.pos")

custom_order_level2 <- c(
  "Granulosa",
  "Theca",
  "Stroma",
  "Mesenchyme",
  "Pericyte",
  "BEC",
  "LEC",
  "Epithelial",
  "Neutrophil", "Monocyte", "Macrophage" , "DC", "ILC", "NK", "NKT", "CD8 NKT", "DNT", "DPT", "B"
)

################################################################################
# Functions
################################################################################

run_prop_test <- function(data, group_by) {
  prop_test <- sc_utils(data)
  permutation_test(
    prop_test, cluster_identity = group_by,
    sample_1 = "wt", sample_2 = "het",
    sample_identity = "Genotype"
  )
}

generate_permutation_plot <- function(..., groups, title, custom_order, max_range = 8) {
  plot_list <- list(...)
  plot_data <- do.call(rbind, lapply(seq_along(plot_list), function(i) {
    df <- plot_list[[i]]@results$permutation
    df$groups <- rep(groups[i], nrow(df))
    df
  }))
  
  plot_data$groups <- factor(plot_data$groups, levels = c("Midage", "Young"))
  plot_data$clusters <- factor(plot_data$clusters, levels = custom_order)
  
  max_value <- max(plot_data$obs_log2FD[is.finite(plot_data$obs_log2FD)], na.rm = TRUE)
  min_value <- min(plot_data$obs_log2FD[is.finite(plot_data$obs_log2FD)], na.rm = TRUE)
  
  plot_data <- plot_data %>%
    mutate(
      obs_log2FD = case_when(
        obs_log2FD == Inf ~ max_value,
        obs_log2FD == -Inf ~ min_value,
        TRUE ~ obs_log2FD
      ),
      is_inf = ifelse(obs_log2FD %in% c(max_value, min_value), "*", NA),
      significance = ifelse(FDR < 0.05 & boot_mean_log2FD > 0, "Increased",
                            ifelse(FDR < 0.05 & boot_mean_log2FD < 0, "Decreased", "n.s.")),
      group_significance = paste(groups, significance, sep = "_")
    )
  
  custom_colors <- c(
    "Young_Increased" = "#ADDFE4", "Midage_Increased" = "#508B8E",
    "Young_Decreased" = "deeppink1", "Midage_Decreased" = "#B10D66",
    "Young_n.s." = "grey", "Midage_n.s." = "grey"
  )
  shape_mapping <- c("Young" = 16, "Midage" = 17)
  
  ggplot(plot_data, aes(x = clusters, y = obs_log2FD, group = interaction(groups, clusters))) +
    theme_bw() +
    geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = group_significance, shape = groups),
                    position = position_dodge2(width = 0.8, padding = 0.5)) +
    geom_text(aes(label = is_inf), nudge_y = 0.2, na.rm = TRUE, size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = custom_colors) +
    scale_shape_manual(values = shape_mapping) +
    coord_flip() + ylim(-max_range, max_range) +
    ggtitle(title) +
    theme(legend.position = "right", axis.text.y = element_text(size = 12))
}

export_frequency <- function(data, group_by, celltypes, file_prefix) {
  freq_table <- prop.table(table(data@meta.data[[group_by]], data@meta.data$Library), margin = 2)
  freq_table <- freq_table[match(celltypes, rownames(freq_table)), , drop = FALSE]
  freqs <- as.data.frame(as.matrix(t(freq_table)))
  write.table(freqs, file = paste0(Sys.Date(), "_", file_prefix, "_proportion_frequency.txt"), quote = FALSE, sep = "\t")
}

################################################################################
# Load data
################################################################################

# Load dataset
load("./Input/2024-10-24_10x_ovary_Benayoun_lab_Foxl2_Seurat_object_with_final_annotation.RData")

# Subset the Seurat object by Age
data_subsets <- list(
  "young" = subset(ovary.Foxl2, subset = Age == "Young"),
  "midage" = subset(ovary.Foxl2, subset = Age == "Midage")
)

################################################################################
# 2. Run proportion tests and generate plots
################################################################################

results_Level1 <- lapply(data_subsets, function(x) run_prop_test(x, "celltype.level1"))
results_Level2 <- lapply(data_subsets, function(x) run_prop_test(x, "celltype.level2"))

plot_Level1 <- generate_permutation_plot(
  results_Level1$young, results_Level1$midage,
  groups = c("Young", "Midage"),
  title = "Foxl2 Level1: WT vs. HET", custom_order = custom_order_level1
)

plot_Level2 <- generate_permutation_plot(
  results_Level2$young, results_Level2$midage,
  groups = c("Young", "Midage"),
  title = "Foxl2 Level2: WT vs. HET", custom_order = custom_order_level2
)

pdf(paste(Sys.Date(), "WT_vs_HET_Level1_Permutation_Plot.pdf", sep = "_"), height = 8, width = 8)
print(plot_Level1)
dev.off()

pdf(paste(Sys.Date(), "WT_vs_HET_Level2_Permutation_Plot.pdf", sep = "_"), height = 12, width = 8)
print(plot_Level2)
dev.off()

# Export frequency tables
export_frequency(ovary.Foxl2, "celltype.level1", custom_order_level1, "Foxl2_Level1")
export_frequency(ovary.Foxl2, "celltype.level2", custom_order_level2, "Foxl2_Level2")

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Foxl2_proportion_test_session_Info.txt", sep =""))
sessionInfo()
sink()
