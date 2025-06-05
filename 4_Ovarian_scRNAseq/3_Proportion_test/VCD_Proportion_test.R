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
# VCD model dataset
################################################################################

################################################################################
# Define color scheme and custom order of celltypes
################################################################################

# Colors

colors_CTL <- colorRampPalette(c("deeppink1", "deeppink4"))(10)
colors_VCD <- colorRampPalette(c("gold", "darkgoldenrod4"))(4)

CTL_colors <- list(
  "3_months"  = colors_CTL[1], "4_months" = "deeppink1", "5_months" = colors_CTL[3],
  "6_months"  = colors_CTL[4], "7_months" = colors_CTL[5], "8_months" = colors_CTL[6],
  "10_months" = colors_CTL[7], "12_months" = colors_CTL[8], "14_months" = colors_CTL[9],
  "20_months" = "deeppink4"
)

VCD_colors <- list(
  "5_months"  = colors_VCD[1], "7_months" = colors_VCD[2],
  "12_months" = colors_VCD[3], "14_months" = colors_VCD[4]
)

# Celltype order

custom_order_level1 <- c("Ptprc.neg", "Ptprc.pos")

custom_order_level2 <- c(
  "Granulosa", "Theca", "Stroma", "Mesenchyme", "Pericyte", "BEC", "LEC", "Epithelial",
  "Neutrophil", "Monocyte", "Macrophage" , "DC", "ILC", "NK", "NKT", "CD8 NKT", "CD4 T", "DNT", "DPT", "B"
)

################################################################################
# Functions
################################################################################

run_prop_test <- function(data, group_by) {
  prop_test <- sc_utils(data)
  permutation_test(prop_test, cluster_identity = group_by, sample_1 = "CTL", sample_2 = "VCD", sample_identity = "Treatment")
}

generate_permutation_plot <- function(..., groups, title, custom_order, max_range = 8) {
  plot_list <- list(...)
  plot_data <- do.call(rbind, lapply(seq_along(plot_list), function(i) {
    df <- plot_list[[i]]@results$permutation
    df$groups <- groups[i]
    df
  }))
  
  plot_data$groups <- factor(plot_data$groups, levels = rev(c("3m_30d", "3m_90d", "10m_30d", "10m_90d")))
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
      significance = case_when(
        FDR < 0.05 & boot_mean_log2FD < 0 ~ "Decreased",
        FDR < 0.05 & boot_mean_log2FD > 0 ~ "Increased",
        TRUE ~ "n.s."
      ),
      group_significance = paste(groups, significance, sep = "_")
    )
  
  custom_colors <- c(
    "3m_30d_Decreased" = "#FF1493", "3m_90d_Decreased" = "#D8107C", 
    "10m_30d_Decreased" = "#A40C5E", "10m_90d_Decreased" = "#970B57", 
    "3m_30d_Increased" = "#FFD700", "3m_90d_Increased" = "#D8B102", 
    "10m_30d_Increased" = "#B18B05", "10m_90d_Increased" = "#8B6508", 
    "3m_30d_n.s." = "grey", "3m_90d_n.s." = "grey", 
    "10m_30d_n.s." = "grey", "10m_90d_n.s." = "grey"
  )
  
  shape_mapping <- c("3m_30d" = 16, "3m_90d" = 17, "10m_30d" = 15, "10m_90d" = 18)
  
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

subset_data <- function(data) {
  list(
    "3m_30d" = subset(data, subset = Age == "3m" & Duration == "30d"),
    "3m_90d" = subset(data, subset = Age == "3m" & Duration == "90d"),
    "10m_30d" = subset(data, subset = Age == "10m" & Duration == "30d"),
    "10m_90d" = subset(data, subset = Age == "10m" & Duration == "90d")
  )
}

# Export frequency table by group (Library) and selected cell types
export_frequency <- function(data, group_by, celltypes, file_prefix) {
  freq_table <- prop.table(table(data@meta.data[[group_by]], data@meta.data$Library), margin = 2)
  freq_table <- freq_table[match(celltypes, rownames(freq_table)), , drop = FALSE]
  freqs <- as.data.frame(as.matrix(t(freq_table)))
  
  out_file <- paste0(Sys.Date(), "_", file_prefix, "_proportion_frequency.txt")
  write.table(freqs, file = out_file, quote = FALSE, sep = "\t")
}

################################################################################
# Load data
################################################################################

# Load dataset
load("./Input/2024-10-24_10x_ovary_Benayoun_lab_VCD_Seurat_object_with_final_annotation.RData")

# Subset the Seurat object by Age and Duration
data_subsets <- subset_data(ovary.VCD)

################################################################################
# Run proportion tests and generate plots
################################################################################

results_Level1 <- lapply(data_subsets, function(x) run_prop_test(x, "celltype.level1"))
results_Level2 <- lapply(data_subsets, function(x) run_prop_test(x, "celltype.level2"))

plot_Level1 <- generate_permutation_plot(
  results_Level1$`3m_30d`, results_Level1$`3m_90d`, results_Level1$`10m_30d`, results_Level1$`10m_90d`,
  groups = c("3m_30d", "3m_90d", "10m_30d", "10m_90d"), 
  title = "VCD cohort: Level1 - CTL vs. VCD", custom_order = custom_order_level1
)

plot_Level2 <- generate_permutation_plot(
  results_Level2$`3m_30d`, results_Level2$`3m_90d`, results_Level2$`10m_30d`, results_Level2$`10m_90d`,
  groups = c("3m_30d", "3m_90d", "10m_30d", "10m_90d"), 
  title = "VCD cohort: Level2 - CTL vs. VCD", custom_order = custom_order_level2
)

# Save plots
pdf(paste(Sys.Date(), "VCD_Level1_Permutation_Plot.pdf", sep = "_"), height = 8, width = 8)
print(plot_Level1)
dev.off()

pdf(paste(Sys.Date(), "VCD_Level2_Permutation_Plot.pdf", sep = "_"), height = 12, width = 8)
print(plot_Level2)
dev.off()

# Export frequency table
export_frequency(
  data = ovary.VCD,
  group_by = "celltype.level1",
  celltypes = custom_order_level1,
  file_prefix = "VCD"
)

export_frequency(
  data = ovary.VCD,
  group_by = "celltype.level2",
  celltypes = custom_order_level2,
  file_prefix = "VCD"
)

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_VCD_proportion_test_session_Info.txt", sep =""))
sessionInfo()
sink()
