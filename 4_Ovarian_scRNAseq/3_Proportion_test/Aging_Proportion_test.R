setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/3_Proportion_test")
options(stringsAsFactors = F)

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(beeswarm)
library(scProportionTest)
library(gridExtra)


################################################################################
# Perform cell type proportion tests
# scProportionTest
# Aging ("AC") model dataset
################################################################################

# Functions

## Generate beeswarm plot
export_frequency <- function(data, group_by, celltypes, level) {
  
  # Generate frequency table
  freq_table <- prop.table(table(data@meta.data[[group_by]], data@meta.data$Library), margin = 2)
  freq_table <- freq_table[match(celltypes, rownames(freq_table)), , drop = FALSE]
  freqs <- as.data.frame(as.matrix(t(freq_table)))
  
  write.table(freqs, file = paste0(Sys.Date(), "_Aging_proportion_frequency.txt"), quote = FALSE, sep = "\t")
}


## Run scProportionTest and generate plot
run_prop_test <- function(data, group_by, level, custom_order = NULL) {
  
  # Run scProportionTest
  prop_test <- sc_utils(data)
  results <- permutation_test(
    prop_test, cluster_identity = group_by, sample_1 = "YF", sample_2 = "OF", 
    sample_identity = "Age"
  )
  
  # Parse results
  plot_data <- results@results$permutation
  
  # Set custom cluster order if provided
  if (!is.null(custom_order)) {
    plot_data$clusters <- factor(plot_data$clusters, levels = custom_order)
  }
  
  # Define significance labels
  plot_data$significance <- factor(ifelse(
    plot_data$FDR < 0.05 & plot_data$boot_mean_log2FD < 0, "FDR<0.05 (Decreased)",
    ifelse(plot_data$FDR < 0.05 & plot_data$boot_mean_log2FD > 0, "FDR<0.05 (Increased)", "n.s.")
  ), levels = c("FDR<0.05 (Decreased)", "FDR<0.05 (Increased)", "n.s."))
  
  # Set filename
  filename <- paste(Sys.Date(), "10x_ovary_Benayoun_lab_AC_celltype_proportion_TEST_FREQUENCY", 
                    level, ".pdf", sep = "_")
  
  # Plot
  p <- ggplot(plot_data, aes(x = clusters, y = obs_log2FD)) + 
    theme_bw() +
    geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance),
                    position = position_dodge(width = -1)) +
    geom_hline(yintercept = 0) + 
    scale_color_manual(values = c(
      "FDR<0.05 (Decreased)" = "deeppink1", 
      "FDR<0.05 (Increased)" = "deeppink4", 
      "n.s." = "grey"
    )) +
    coord_flip() + 
    ylim(-8, 8) + 
    ggtitle(paste("Proportion Test -", level)) +
    theme(legend.title = element_blank())
  
  # Save the plot to PDF
  pdf(filename, height = 10, width = 7)
  print(p)
  dev.off()
}

################################################################################
# 1. Load data
################################################################################

# Load dataset
load("./Input/2024-10-15_10x_ovary_Benayoun_lab_AC_Seurat_object_with_final_annotation.RData")

# Custom order for celltype.level1
custom_order_level1 <- c("Ptprc.neg", "Ptprc.pos")

# Custom order for celltype.level2
custom_order_level2 <- c(
  "Granulosa", "Theca", "Stroma", "Mesenchyme", "Pericyte", "BEC", "LEC", 
  "Epithelial", "Neutrophil", "Monocyte", "Macrophage", "DC", "ILC", "NK", 
  "NKT", "CD8 NKT", "CD8 T", "CD4 T", "DNT", "DPT", "B"
)

################################################################################
# 2. Run proportion tests and generate plots
################################################################################

# Generate beeswarm plots
export_frequency(ovary.AC, "celltype.level1", custom_order_level1, "LEVEL1")
export_frequency(ovary.AC, "celltype.level2", custom_order_level2, "LEVEL2")

# Run scProportionTest for both levels
run_prop_test(ovary.AC, group_by = "celltype.level1", level = "LEVEL1", 
              custom_order = custom_order_level1)
run_prop_test(ovary.AC, group_by = "celltype.level2", level = "LEVEL2", 
              custom_order = custom_order_level2)

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Aging_proportion_test_session_Info.txt", sep =""))
sessionInfo()
sink()
