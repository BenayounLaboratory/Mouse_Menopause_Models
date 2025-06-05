setwd("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/RNAscope/")

# Load libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(beeswarm)
library(outliers)
library(stringr)

################################################################################
# RNAscope data - plot counts
# Remove outliers
################################################################################

############################
# Functions
############################

remove_grubbs_outliers <- function(df, value_col = "Proportion", group_col = "Group", marker_col = "Marker") {
  df_clean <- df %>%
    group_by(across(all_of(c(group_col, marker_col)))) %>%
    group_modify(~ {
      values <- .x[[value_col]]
      while (length(values) >= 3) {
        g_test <- tryCatch(grubbs.test(values), error = function(e) return(NULL))
        if (is.null(g_test) || g_test$p.value >= 0.05) break
        outlier_val <- as.numeric(str_extract(g_test$alternative, "\\d+\\.*\\d*"))
        values <- values[values != outlier_val]
      }
      .x[.x[[value_col]] %in% values, , drop = FALSE]
    }) %>%
    ungroup()
  return(df_clean)
}

############################
# Aging & VCD
############################

plot_marker_proportions <- function(df, model_name, group_assign_fn, color_map, comparisons) {
  for (this_set in unique(df$Set)) {
    
    # Subset and preprocess
    data_set <- df %>% filter(Set == this_set)
    detections_df <- data_set %>%
      filter(Category == "Detections") %>%
      select(Sample_ID, Detections = Count)
    
    df_proportion <- data_set %>%
      filter(!Category %in% c("Detections", "Cd4")) %>%
      left_join(detections_df, by = "Sample_ID") %>%
      mutate(
        Proportion = Count / Detections,
        Marker = Category,
        Group = group_assign_fn(Sample_ID)
      ) %>%
      select(Sample_ID, Marker, Group, Proportion)
    
    df_proportion$Group <- factor(df_proportion$Group, levels = names(color_map))
    df_proportion$Marker <- factor(df_proportion$Marker, levels = unique(df_proportion$Marker))
    
    marker_list <- unique(df_proportion$Marker)
    n_col <- 4
    n_row <- ceiling(length(marker_list) / n_col)
    
    pdf(paste0(Sys.Date(), "_", model_name, "_Set", this_set, "_marker_proportion_boxplot.pdf"), width = 14, height = 7)
    par(mfrow = c(n_row, n_col), mar = c(4.5, 4.5, 3, 1))
    
    for (marker_name in marker_list) {
      marker_data <- df_proportion %>% filter(Marker == marker_name)
      y_max <- max(marker_data$Proportion, na.rm = TRUE) * 1.2
      
      boxplot(Proportion ~ Group, data = marker_data,
              col = color_map,
              ylim = c(0, y_max),
              outline = TRUE, las = 1,
              main = marker_name,
              ylab = "Proportion", xlab = "Group")
      
      beeswarm(Proportion ~ Group, data = marker_data,
               pch = 1, col = "black", add = TRUE, cex = 0.9)
      
      for (comp in comparisons) {
        g1 <- comp[1]
        g2 <- comp[2]
        test_data <- marker_data %>% filter(Group %in% c(g1, g2))
        p_val <- tryCatch(
          wilcox.test(Proportion ~ Group, data = test_data)$p.value,
          error = function(e) NA
        )
        mid_x <- (which(levels(marker_data$Group) == g1) + which(levels(marker_data$Group) == g2)) / 2
        text(mid_x, y_max * 0.95, paste0(g1, " vs ", g2, ": p=", format(p_val, digits = 3)), cex = 0.7)
      }
    }
    dev.off()
  }
}

# Aging model

Aging.model <- read.delim("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/RNAscope/Raw_count_data/Aging_model_combined.txt", sep = "\t", header = TRUE)

plot_marker_proportions(
  Aging.model,
  model_name = "Aging",
  group_assign_fn = function(id) ifelse(grepl("^YF", id), "YF", "OF"),
  color_map = c("YF" = "deeppink1", "OF" = "deeppink4"),
  comparisons = list(c("YF", "OF"))
)


VCD.model <- read.delim("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/RNAscope/Raw_count_data/VCD_model_combined.txt", sep = "\t", header = TRUE)

plot_marker_proportions(
  VCD.model,
  model_name = "VCD",
  group_assign_fn = function(id) ifelse(grepl("^CTL", id), "CTL", "VCD"),
  color_map = c("CTL" = "deeppink1", "VCD" = "gold"),
  comparisons = list(c("CTL", "VCD"))
)

############################
# Foxl2 haploinsufficiency
############################

Foxl2.model <- read.delim("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/RNAscope/Raw_count_data/Foxl2_model_combined.txt", sep = "\t", header = TRUE)

plot_marker_proportions_grouped <- function(df, model_name, color_map, comparisons) {
  for (this_set in unique(df$Set)) {
    
    data_set <- df %>% filter(Set == this_set)
    
    detections_df <- data_set %>%
      filter(Category == "Detections") %>%
      select(Sample_ID, Detections = Count)
    
    df_proportion <- data_set %>%
      filter(!Category %in% c("Detections", "Cd4")) %>%
      left_join(detections_df, by = "Sample_ID") %>%
      mutate(
        Proportion = Count / Detections,
        Marker = Category
      ) %>%
      select(Sample_ID, Marker, Group, Proportion)
    
    df_proportion$Group <- factor(df_proportion$Group, levels = names(color_map))
    df_proportion$Marker <- factor(df_proportion$Marker, levels = unique(df_proportion$Marker))
    
    marker_list <- unique(df_proportion$Marker)
    n_col <- 4
    n_row <- ceiling(length(marker_list) / n_col)
    
    pdf(paste0(Sys.Date(), "_", model_name, "_Set", this_set, "_marker_proportion_boxplot.pdf"), width = 14, height = 7)
    par(mfrow = c(n_row, n_col), mar = c(4.5, 4.5, 3, 1))
    
    for (marker_name in marker_list) {
      marker_data <- df_proportion %>% filter(Marker == marker_name)
      y_max <- max(marker_data$Proportion, na.rm = TRUE) * 1.2
      
      boxplot(Proportion ~ Group, data = marker_data,
              col = color_map,
              ylim = c(0, y_max),
              outline = TRUE, las = 1,
              main = marker_name,
              ylab = "Proportion", xlab = "Group")
      
      beeswarm(Proportion ~ Group, data = marker_data,
               pch = 18, col = "black", add = TRUE, cex = 0.9)
      
      # Pairwise Wilcoxon for defined comparisons
      for (comp in comparisons) {
        g1 <- comp[1]
        g2 <- comp[2]
        test_data <- marker_data %>% filter(Group %in% c(g1, g2))
        p_val <- tryCatch(
          wilcox.test(Proportion ~ Group, data = test_data)$p.value,
          error = function(e) NA
        )
        mid_x <- (which(levels(marker_data$Group) == g1) + which(levels(marker_data$Group) == g2)) / 2
        text(mid_x, y_max * 0.95, paste0(g1, " vs ", g2, ": p=", format(p_val, digits = 3)), cex = 0.7)
      }
    }
    
    dev.off()
  }
}

plot_marker_proportions_grouped(
  df = Foxl2.model,
  model_name = "Foxl2",
  color_map = c("Young_wt" = "deeppink1", "Young_het" = "springgreen",
                "Old_wt" = "deeppink1", "Old_het" = "springgreen"),
  comparisons = list(c("Young_wt", "Young_het"), c("Old_wt", "Old_het"))
)

################################################################################
sink(file = paste(Sys.Date(),"RNAscope_quantification_session_Info.txt",sep="_"))
sessionInfo()
sink()
