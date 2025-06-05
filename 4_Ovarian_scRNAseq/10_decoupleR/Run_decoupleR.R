setwd('/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/For_github/decoupleR')
options(stringsAsFactors = F)

library("decoupleR")
library("OmnipathR")

library(dplyr)
library(tidyr)
library(tibble)

library(ggplot2)
library(ComplexHeatmap)
library(circlize)

theme_set(theme_bw())   

rm(list = ls())

######################################################
# Run decouple R
# Use pseudobulked datasets for inference
# Aging, VCD and Foxl2 datasets
######################################################

######################################################
# 1. Load datasets
######################################################

# Aging
load("/Volumes/jinho01/Benayoun_lab/Projects/Menopause_model_project/Data/Ovarian_scRNAseq/PB_DESeq2/Without_scTE/AC/2025-01-14/DESeq2_results/2025-01-14_10x_ovary_Benayoun_lab_AC_PB_DESeq2_object.RData")
deseq.res.list.aging <- deseq.res.list
rm(deseq.res.list)

# VCD
load("/Volumes/jinho01/Benayoun_lab/Projects/Menopause_model_project/Data/Ovarian_scRNAseq/PB_DESeq2/Without_scTE/VCD/2025-01-15/DESeq2_results/2025-01-15_10x_ovary_Benayoun_lab_VCD_PB_DESeq2_object.RData")
deseq.res.list.vcd <- deseq.res.list
rm(deseq.res.list)

# Foxl2
load("/Volumes/jinho01/Benayoun_lab/Projects/Menopause_model_project/Data/Ovarian_scRNAseq/PB_DESeq2/Without_scTE/Foxl2/2025-01-15/DESeq2_results/2025-01-15_10x_ovary_Benayoun_lab_Foxl2_NULL_PB_DESeq2_object.RData")
deseq.res.list.foxl2 <- deseq.res.list
rm(deseq.res.list)

# Mouse CollecTRI network
mouse.net <- get_collectri(organism='mouse', split_complexes=FALSE)

######################################################
# 2. Run TF activity inference
######################################################

# Function to run decoupleR
run_decoupleR <- function(deseq_list, model_name) {
  
  result <- data.frame(matrix(0, 0, 8))
  colnames(result) <- c("statistic", "source", "condition", "score", "p_value", "Reg_Rank", "Cell_Type", "Model")
  
  for (i in seq_along(deseq_list)) {
    # Run fgsea
    contrast_acts <- run_fgsea(mat = deseq_list[[i]][, 'stat', drop=FALSE],
                               net = mouse.net,
                               minsize = 5)
    
    # Filter norm_fgsea results
    tf.scores <- contrast_acts[contrast_acts$statistic == 'norm_fgsea',]
    tf.sig <- tf.scores[tf.scores$p_value < 0.05,]
    tf.sig$Reg_Rank <- NA
    tf.sig$Cell_Type <- names(deseq_list)[i]
    tf.sig$Model <- model_name
    
    # Keep only expressed TFs
    tf.sig <- tf.sig[tf.sig$source %in% row.names(deseq_list[[i]]),]
    
    # Rank assignment
    msk <- tf.sig$score > 0
    tf.sig$Reg_Rank[msk]  <- round(rank(-tf.sig[msk, 'score']))
    tf.sig$Reg_Rank[!msk] <- round(rank(-abs(tf.sig[!msk, 'score'])))
    
    result <- rbind(result, tf.sig)
  }
  return(result)
}

Aging_TFs  <- run_decoupleR(deseq.res.list.aging, "Aging")
VCD_TFs    <- run_decoupleR(deseq.res.list.vcd, "VCD")
Foxl2_TFs  <- run_decoupleR(deseq.res.list.foxl2, "Foxl2")

combined_TF_res <- bind_rows(Aging_TFs, VCD_TFs, Foxl2_TFs)

# Save results
save(combined_TF_res, file = paste0(Sys.Date(), "_combined_decoupleR_TF_results.RData"))

######################################################
# 3. Assess recurrent TFs
######################################################

# Identify common cell types
Aging.unique.celltype <- unique(Aging.TFs.res$Cell_Type)
VCD.unique.celltype   <- unique(VCD.TFs.res$Cell_Type)
Foxl2.unique.celltype <- unique(Foxl2.TFs.res$Cell_Type)

common.celltype <- intersect(Aging.unique.celltype, intersect(VCD.unique.celltype, Foxl2.unique.celltype))
common.celltype
# [1] "Granulosa"  "Theca"      "Stroma"     "BEC"        "Epithelial" "DNT"       

#############################
# 3-1. Extract TFs
#############################

# Extract TFs for each cell type
extract_TFs_by_celltype <- function(TF_results, dataset_name) {
  TF_results %>%
    select(source, Cell_Type) %>%
    mutate(dataset = dataset_name)
}

# Process each dataset
Aging_TFs <- extract_TFs_by_celltype(Aging.TFs.res, "Aging")
VCD_TFs <- extract_TFs_by_celltype(VCD.TFs.res, "VCD")
Foxl2_TFs <- extract_TFs_by_celltype(Foxl2.TFs.res, "Foxl2")

# Combine datasets
all_TFs <- bind_rows(Aging_TFs, Foxl2_TFs, VCD_TFs)

filtered_TFs <- all_TFs %>%
  filter(Cell_Type %in% common.celltype)

# unique(filtered_TFs$Cell_Type)

# Get list of unique cell types
cell_types <- unique(filtered_TFs$Cell_Type)

#############################
# 3-2. Identify TFs - Found in at least two datasets
#############################

# Initialize lists
merged_TF_activity_list <- list()

for (cell_type in cell_types) {
  
  # Extract TFs for cell type
  TFs_Aging <- Aging.TFs.res %>% filter(Cell_Type == cell_type) %>% select(source, score) %>% rename(score_Aging = score)
  TFs_Foxl2 <- Foxl2.TFs.res %>% filter(Cell_Type == cell_type) %>% select(source, score) %>% rename(score_Foxl2 = score)
  TFs_VCD <- VCD.TFs.res %>% filter(Cell_Type == cell_type) %>% select(source, score) %>% rename(score_VCD = score)
  
  # Merge TF scores by TF name
  merged_TF_activity <- full_join(TFs_Aging, TFs_Foxl2, by = "source") %>%
    full_join(., TFs_VCD, by = "source")
  
  # Filter for TFs that exist in at least two datasets
  merged_TF_activity <- merged_TF_activity %>%
    filter(rowSums(!is.na(merged_TF_activity[, -1])) >= 2)
  
  # Rank TFs within each dataset
  merged_TF_activity <- merged_TF_activity %>%
    mutate(rank_Aging = rank(-score_Aging, na.last = "keep"),
           rank_Foxl2 = rank(-score_Foxl2, na.last = "keep"),
           rank_VCD = rank(-score_VCD, na.last = "keep"))
  
  # Store results
  merged_TF_activity_list[[cell_type]] <- merged_TF_activity
  
}

# Save results
save(merged_TF_activity_list, file = "common_TFs_By_CellType.RData")

write.table(do.call(rbind, merged_TF_activity_list), "Common_TFs_By_CellType.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#############################
# 3-3. Generate plots
#############################

# Define color palette
color_palette <- c("darkblue", "dodgerblue4", "dodgerblue3", "dodgerblue1",
                   "white", "lightcoral", "brown1", "firebrick2", "firebrick4")

# Set minimum and maximum limits for score values
score_limits <- c(-3, 3)

# Function to generate heatmap for a given TF dataset
generate_tf_heatmap <- function(tf_data, cell_type) {
  # Skip if there are no TFs
  if (nrow(tf_data) == 0) {
    message("Skipping heatmap for ", cell_type, " - No TFs found.")
    return(NULL)
  }
  
  # Sort TFs alphabetically
  tf_data <- tf_data %>% arrange(desc(source))
  
  # Convert to long format plotting
  heatmap_data <- tf_data %>%
    pivot_longer(cols = starts_with("score_"), names_to = "Dataset", values_to = "Score") %>%
    mutate(Dataset = factor(Dataset, levels = c("score_Aging", "score_VCD", "score_Foxl2")))
  
  # Generate heatmap
  pdf_filename <- paste0(Sys.Date(), "_", cell_type, "_TF_Activity_Heatmap.pdf")
  pdf(pdf_filename, width = 8, height = 12)
  print(ggplot(heatmap_data, aes(x = Dataset, y = factor(source, levels = tf_data$source), fill = Score)) +
    geom_tile(color = "white") + 
    scale_fill_gradientn(
      colours = color_palette, 
      limits = score_limits, 
      na.value = "gray"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10)
    ) +
    labs(
      x = "Dataset",
      y = "Transcription Factor",
      fill = "Activity Score"
    ))
  dev.off()
}

for (cell_type in names(merged_TF_activity_list)) {
  generate_tf_heatmap(merged_TF_activity_list[[cell_type]], cell_type)
}

######################################################
sink(file = paste(Sys.Date(), "_Run_decoupleR_session_info.txt"))
sessionInfo()
sink()
