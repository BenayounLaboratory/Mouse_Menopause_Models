setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/11_Transcriptional_clock")

library(Seurat)
library(glmnet)
library(caret)
library(dplyr)

library(ggplot2)
library(gridExtra)

library(beeswarm)

rm(list=ls())

options(mc.cores = 4)
set.seed(123)       

###############################################################
# Transcriptional clock - Granulosa, Theca
###############################################################

###############################################################
# 1. Load and prep objects
###############################################################

load("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/11_Transcriptional_clock/Input/2024-10-15_10x_ovary_Benayoun_lab_AC_Seurat_object_with_final_annotation.RData")
load("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/11_Transcriptional_clock/Input/2025-02-13_10x_ovary_Benayoun_lab_AC_MA_Seurat_object_with_final_annotation.RData")
load("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/11_Transcriptional_clock/Input/2024-10-24_10x_ovary_Benayoun_lab_VCD_Seurat_object_with_final_annotation.RData")
load("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/11_Transcriptional_clock/Input/2024-10-24_10x_ovary_Benayoun_lab_Foxl2_Seurat_object_with_final_annotation.RData")

################################################################################
# 2. Subset Granulosa and Theca cells separately
################################################################################

# Subset for Granulosa
sub.ovary.AC.granulosa <- subset(ovary.AC, subset = celltype.level2 == "Granulosa")
sub.ovary.AC.MA.granulosa <- subset(ovary.AC.MA, subset = celltype.level2 == "Granulosa")
sub.ovary.VCD.granulosa <- subset(ovary.VCD, subset = celltype.level2 == "Granulosa")
sub.ovary.Foxl2.granulosa <- subset(ovary.Foxl2, subset = celltype.level2 == "Granulosa")

seurat.Granulosa.combined <- merge(sub.ovary.AC.granulosa,
                                   y = list(sub.ovary.AC.MA.granulosa,
                                            sub.ovary.VCD.granulosa,
                                            sub.ovary.Foxl2.granulosa),
                                   project = "Granulosa_Combined")

# Subset for Theca
sub.ovary.AC.theca <- subset(ovary.AC, subset = celltype.level2 == "Theca")
sub.ovary.AC.MA.theca <- subset(ovary.AC.MA, subset = celltype.level2 == "Theca")
sub.ovary.VCD.theca <- subset(ovary.VCD, subset = celltype.level2 == "Theca")
sub.ovary.Foxl2.theca <- subset(ovary.Foxl2, subset = celltype.level2 == "Theca")

seurat.Theca.combined <- merge(sub.ovary.AC.theca,
                               y = list(sub.ovary.AC.MA.theca,
                                        sub.ovary.VCD.theca,
                                        sub.ovary.Foxl2.theca),
                               project = "Theca_Combined")

################################################################################
# 3. Update metadata
################################################################################

# Function to update metadata
update_metadata <- function(seurat_obj) {
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    mutate(
      Treatment = NA_character_,
      Age = NA_character_,
      Duration = NA_character_,
      Genotype = NA_character_
    ) %>%
    mutate(
      Treatment = case_when(
        grepl("^CTL", Library) ~ "CTL",
        grepl("^VCD", Library) ~ "VCD",
        TRUE ~ NA_character_
      ),
      Age = case_when(
        grepl("^YF", Library) ~ "4m",
        grepl("^OF", Library) ~ "20m",
        grepl("CTL_3m_30d|VCD_3m_30d", Library) ~ "4m",
        grepl("CTL_3m_90d|VCD_3m_90d", Library) ~ "6m",
        grepl("CTL_10m_30d|VCD_10m_30d", Library) ~ "11m",
        grepl("CTL_10m_90d|VCD_10m_90d", Library) ~ "13m",
        grepl("Foxl2_wt_young|Foxl2_het_young", Library) ~ "3m",
        grepl("Foxl2_wt_midage|Foxl2_het_midage", Library) ~ "8m",
        TRUE ~ NA_character_
      ),
      Duration = case_when(
        grepl("VCD_.*_30d", Library) ~ "30d",
        grepl("VCD_.*_90d", Library) ~ "90d",
        TRUE ~ NA_character_
      ),
      Genotype = case_when(
        grepl("Foxl2_wt", Library) ~ "wt",
        grepl("Foxl2_het", Library) ~ "het",
        TRUE ~ NA_character_
      ),
      Age_num = case_when(
        grepl("^YF", Library) ~ 17,
        grepl("^OF", Library) ~ 87,
        grepl("^13weeks", Library) ~ 13,
        grepl("^76weeks", Library) ~ 76,
        grepl("^86weeks", Library) ~ 86,
        grepl("CTL_3m_30d|VCD_3m_30d", Library) ~ 20,
        grepl("CTL_3m_90d|VCD_3m_90d", Library) ~ 29,
        grepl("CTL_10m_30d|VCD_10m_30d", Library) ~ 48,
        grepl("CTL_10m_90d|VCD_10m_90d", Library) ~ 57,
        grepl("Foxl2_wt_young|Foxl2_het_young", Library) ~ 15,
        grepl("Foxl2_wt_midage|Foxl2_het_midage", Library) ~ 39,
        TRUE ~ NA_real_
      )
    )
  return(seurat_obj)
}

# Apply function
seurat.Granulosa.combined <- update_metadata(seurat.Granulosa.combined)
seurat.Theca.combined <- update_metadata(seurat.Theca.combined)

# Save updated objects
save(seurat.Granulosa.combined, file = paste0(Sys.Date(), "_10x_ovary_Granulosa_combined_Seurat_object.RData"))
save(seurat.Theca.combined, file = paste0(Sys.Date(), "_10x_ovary_Theca_combined_Seurat_object.RData"))

DefaultAssay(seurat.Granulosa.combined) <- "RNA"
DefaultAssay(seurat.Theca.combined) <- "RNA"

###############################################################
# 2. Function to train clock
###############################################################

train_lasso_model <- function(seurat_obj, cell_type) {
  
  # Subset dataset by controls only (YF, OF, weeks, wt, CTL)
  seurat_controls <- subset(seurat_obj,
                            cells = grep("YF|OF|weeks|wt|CTL", colnames(seurat_obj), value = TRUE))
  
  # Extract metadata
  metadata <- seurat_controls@meta.data
  metadata$Age_num <- as.numeric(metadata$Age_num)
  
  # Partition data
  train_indices <- createDataPartition(metadata$Age_num, p = 0.75, list = FALSE)
  
  # Get cell names for training and testing
  train_cells <- rownames(metadata)[train_indices]
  test_cells <- rownames(metadata)[-train_indices]
  
  # Subset Seurat objects for training and testing
  train_obj <- subset(seurat_controls, cells = train_cells)
  test_obj <- subset(seurat_controls, cells = test_cells)
  
  ###############################################################
  # 3. Train LASSO Model
  ###############################################################
  
  # Extract expression data
  train_data <- GetAssayData(train_obj, slot = "data", assay = "RNA")
  test_data <- GetAssayData(test_obj, slot = "data", assay = "RNA")
  
  train_age <- train_obj@meta.data$Age_num
  test_age <- test_obj@meta.data$Age_num
  
  # Set the number of folds for cross-validation
  k_folds <- 5
  folds <- createFolds(train_age, k = k_folds, list = TRUE)
  
  fold_metrics <- list()
  
  # Perform k-fold cross-validation
  for (fold in 1:k_folds) {
    train_indices <- unlist(folds[-fold])
    val_indices <- unlist(folds[fold])
    
    train_matrix <- as.matrix(t(train_data[, train_indices]))
    val_matrix <- as.matrix(t(train_data[, val_indices]))
    
    train_age_fold <- train_age[train_indices]
    val_age_fold <- train_age[val_indices]
    
    # Train LASSO model
    inner_cv_model <- cv.glmnet(
      x = train_matrix,
      y = train_age_fold,
      alpha = 1,
      family = "gaussian"
    )
    
    best_lambda <- inner_cv_model$lambda.min
    
    final_model <- glmnet(
      x = train_matrix,
      y = train_age_fold,
      alpha = 1,
      lambda = best_lambda
    )
    
    # Predict on validation set
    predicted_val_age <- predict(final_model, s = best_lambda, newx = val_matrix)
    
    # Compute performance metrics
    mae <- mean(abs(predicted_val_age - val_age_fold))
    mse <- mean((predicted_val_age - val_age_fold)^2)
    rmse <- sqrt(mse)
    r_squared <- 1 - sum((val_age_fold - predicted_val_age)^2) / sum((val_age_fold - mean(val_age_fold))^2)
    correlation <- cor(predicted_val_age, val_age_fold, method = "spearman")
    
    fold_metrics[[fold]] <- c(MAE = mae, MSE = mse, RMSE = rmse, R2 = r_squared, Correlation = correlation)
  }
  
  # Compute average cross-validation performance
  avg_metrics <- colMeans(do.call(rbind, fold_metrics), na.rm = TRUE)
  
  # Train final model on full training data
  final_train_matrix <- as.matrix(t(train_data))
  final_test_matrix <- as.matrix(t(test_data))
  
  final_cv_model <- cv.glmnet(
    x = final_train_matrix,
    y = train_age,
    alpha = 1,
    family = "gaussian"
  )
  
  final_lambda <- final_cv_model$lambda.min
  
  final_model <- glmnet(
    x = final_train_matrix,
    y = train_age,
    alpha = 1,
    lambda = final_lambda
  )
  
  # Predict on test set
  predicted_test_age <- predict(final_model, s = final_lambda, newx = final_test_matrix)
  
  # Compute test performance
  test_mae <- mean(abs(predicted_test_age - test_age))
  test_mse <- mean((predicted_test_age - test_age)^2)
  test_rmse <- sqrt(test_mse)
  test_r_squared <- 1 - sum((test_age - predicted_test_age)^2) / sum((test_age - mean(test_age))^2)
  
  # Compute test correlation 
  test_correlation <- cor.test(predicted_test_age, test_age, method = "spearman")
  
  # Extract features (genes with nonzero coefficients)
  feature_coefficients <- as.matrix(coef(final_model, s = final_lambda))
  lasso_features <- feature_coefficients[feature_coefficients[,1] != 0, , drop = FALSE]
  
  # Save features and coefficients
  write.table(lasso_features, 
              file = paste0(Sys.Date(), "_", cell_type, "_Lasso_nonzero_features.txt"), 
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
  
  # Store performance results
  performance_summary <- data.frame(
    CellType = cell_type,
    MAE_CV = avg_metrics["MAE"],
    MSE_CV = avg_metrics["MSE"],
    RMSE_CV = avg_metrics["RMSE"],
    R2_CV = avg_metrics["R2"],
    Correlation_CV = avg_metrics["Correlation"],
    MAE_Test = test_mae,
    MSE_Test = test_mse,
    RMSE_Test = test_rmse,
    R2_Test = test_r_squared,
    Correlation_Test_Rho = test_correlation$estimate, 
    Correlation_Test_pval = test_correlation$p.value
  )
  
  # Save performance summary 
  write.table(performance_summary, file = paste0(Sys.Date(), "_", cell_type, "_Test_Set_Performance_Summary.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  # Save predictions
  predictions_df <- data.frame(
    Cell_ID = colnames(test_data),
    Actual_Age = test_age,
    Predicted_Age = as.numeric(predicted_test_age),
    CellType = cell_type
  )
  
  write.table(predictions_df, file = paste0(Sys.Date(), "_", cell_type, "_Test_Set_Predictions.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  # Save model
  save(final_model, file = paste0(Sys.Date(), "_", cell_type, "_lasso_model.RData"))
  
  return(performance_summary)
}

###############################################################
# 4. Train Granulosa & Theca Transcriptional Clocks
###############################################################

granulosa_performance <- train_lasso_model(seurat.Granulosa.combined, "Granulosa")
theca_performance <- train_lasso_model(seurat.Theca.combined, "Theca")

###############################################################
# 5. Assess model performance - chronological vs. predicted age scatter plot
###############################################################

# Load predictions

granulosa.predictions <- read.table("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/11_Transcriptional_clock/Input/2025-03-19_Granulosa_Test_Set_Predictions.txt",
                                    sep = "\t", header = TRUE)
theca.predictions <- read.table("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/11_Transcriptional_clock/Input/2025-03-19_Theca_Test_Set_Predictions.txt",
                                sep = "\t", header = TRUE)

predictions.combined <- rbind(granulosa.predictions, theca.predictions)

# Function to assign model
assign_model <- function(cell_id) {
  if (grepl("YF|OF|weeks", cell_id)) {
    return("AC")
  } else if (grepl("VCD|CTL", cell_id)) {
    return("VCD")
  } else if (grepl("wt|het", cell_id)) {
    return("Foxl2")
  } else {
    return(NA) 
  }
}

# Model assignment
predictions.combined <- predictions.combined %>%
  mutate(Model = sapply(Cell_ID, assign_model))

# Shapes mapping
shape_mapping <- c("AC" = 16, "VCD" = 17, "Foxl2" = 18) 

# Generate scatter plots per cell type
plot_list <- list()

for (cell_type in unique(predictions.combined$CellType)) {
  
  plot <- ggplot(predictions.combined %>% filter(CellType == cell_type), 
                 aes(x = Actual_Age, y = Predicted_Age, shape = Model, color = Model)) +
    geom_point(size = 3, alpha = 0.5) +  
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") +
    scale_shape_manual(values = shape_mapping) +  
    coord_cartesian(xlim = c(0, 100), ylim = c(-20, 120)) + 
    theme_minimal() +
    labs(title = paste("Predicted vs. Actual Age -", cell_type),
         x = "Actual Age", y = "Predicted Age") +
    theme(legend.title = element_blank()) 
  
  plot_list[[cell_type]] <- plot
}

# Save plots
for (cell_type in names(plot_list)) {
  ggsave(filename = paste0(Sys.Date(), "_Predicted_vs_Actual_Age_", cell_type, ".pdf"), 
         plot = plot_list[[cell_type]], width = 7, height = 7)
}

###############################################################
# 6. Test model for test groups
###############################################################

# Load models
load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/11_Transcriptional_clock/Input/2025-03-19_Granulosa_lasso_model.RData")
granulosa.model <- final_model
rm(final_model)

load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/11_Transcriptional_clock/Input/2025-03-19_Theca_lasso_model.RData")
theca.model <- final_model
rm(final_model)

# Function to predict ages
predict_test_groups <- function(seurat_obj, trained_model, cell_type) {
  
  # Subset dataset by test groups only
  test_obj <- subset(seurat_obj,
                     cells = grep("het|VCD", colnames(seurat_obj), value = TRUE))
  
  # Extract test data
  test_data <- GetAssayData(test_obj, slot = "data")
  
  # Get metadata
  metadata_test <- test_obj@meta.data
  
  # Predict age
  predictions <- predict(trained_model, 
                         s = trained_model$lambda, 
                         newx = t(test_data)) 
  
  # Format predictions for merging
  predictions_df <- data.frame(
    Cell_ID = colnames(test_data),
    Predicted_Age = as.numeric(predictions)
  )
  
  metadata_df <- data.frame(
    Cell_ID = rownames(metadata_test),
    Actual_Age = metadata_test$Age_num
  )
  
  # Merge predictions with actual ages
  predictions_merged <- merge(predictions_df, metadata_df, by = "Cell_ID")
  
  # Save results
  write.table(predictions_merged, file = paste0(Sys.Date(), "_", cell_type, "_test_group_Predictions.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  return(predictions_merged)
}

# Apply function
Granulosa_test_predictions <- predict_test_groups(seurat.Granulosa.combined, 
                                                  granulosa.model, 
                                                  "Granulosa")
Theca_test_predictions <- predict_test_groups(seurat.Theca.combined, 
                                              theca.model, 
                                              "Theca")

Granulosa.predictions.combined <- rbind(granulosa.predictions[,c("Cell_ID", "Actual_Age", "Predicted_Age")], Granulosa_test_predictions)
Theca.predictions.combined <- rbind(theca.predictions[,c("Cell_ID", "Actual_Age", "Predicted_Age")], Theca_test_predictions)

# Assign models
Granulosa.predictions.combined <- Granulosa.predictions.combined %>%
  mutate(Model = sapply(Cell_ID, assign_model))

Theca.predictions.combined <- Theca.predictions.combined %>%
  mutate(Model = sapply(Cell_ID, assign_model))

# Save predictions
write.table(Granulosa.predictions.combined, file = paste0(Sys.Date(), "_MeMo_Granulosa_clock_all_predictions_combined.txt"),
            sep = "\t", quote = FALSE)
write.table(Theca.predictions.combined, file = paste0(Sys.Date(), "_MeMo_Theca_clock_all_predictions_combined.txt"),
            sep = "\t", quote = FALSE)

# Generate scatter plots per cell type

# VCD model

Granulosa.predictions.VCD <- Granulosa.predictions.combined[Granulosa.predictions.combined$Model == "VCD",]
Theca.predictions.VCD <- Theca.predictions.combined[Theca.predictions.combined$Model == "VCD",]

# Add offset & add color
Granulosa.predictions.VCD <- Granulosa.predictions.VCD %>%
  mutate(Color = ifelse(grepl("CTL", Cell_ID), "deeppink", "gold"),
         X_offset = ifelse(grepl("CTL", Cell_ID), -0.5, 0.5))
Theca.predictions.VCD <- Theca.predictions.VCD %>%
  mutate(Color = ifelse(grepl("CTL", Cell_ID), "deeppink", "gold"),
         X_offset = ifelse(grepl("CTL", Cell_ID), -0.5, 0.5))

pdf(file = paste0(Sys.Date(), "_MeMo_Granulosa_clock_predictions_VCD_scatter_plot.pdf"), width = 7, height = 7)
ggplot(Granulosa.predictions.VCD, 
       aes(x = Actual_Age + X_offset, y = Predicted_Age, color = Color)) +
  geom_point(shape = 17, size = 3, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") + 
  scale_color_identity() + 
  theme_minimal() +
  labs(title = "Predicted vs. Actual Age - Granulosa",
       x = "Actual Age", y = "Predicted Age") +
  coord_cartesian(xlim = c(0, 100), ylim = c(-20, 120)) + 
  theme(legend.position = "none")  
dev.off()

pdf(file = paste0(Sys.Date(), "_MeMo_Theca_clock_predictions_VCD_scatter_plot.pdf"), width = 7, height = 7)
ggplot(Theca.predictions.VCD, 
       aes(x = Actual_Age + X_offset, y = Predicted_Age, color = Color)) +
  geom_point(shape = 17, size = 3, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") + 
  scale_color_identity() + 
  theme_minimal() +
  labs(title = "Predicted vs. Actual Age - Theca",
       x = "Actual Age", y = "Predicted Age") +
  coord_cartesian(xlim = c(0, 100), ylim = c(-20, 120)) + 
  theme(legend.position = "none")  
dev.off()

# Foxl2 model

Granulosa.predictions.Foxl2 <- Granulosa.predictions.combined[Granulosa.predictions.combined$Model == "Foxl2",]
Theca.predictions.Foxl2 <- Theca.predictions.combined[Theca.predictions.combined$Model == "Foxl2",]

# Add offset & add color
Granulosa.predictions.Foxl2 <- Granulosa.predictions.Foxl2 %>%
  mutate(Color = ifelse(grepl("wt", Cell_ID), "deeppink", "springgreen"),
         X_offset = ifelse(grepl("wt", Cell_ID), -0.5, 0.5))
Theca.predictions.Foxl2 <- Theca.predictions.Foxl2 %>%
  mutate(Color = ifelse(grepl("wt", Cell_ID), "deeppink", "springgreen"),
         X_offset = ifelse(grepl("wt", Cell_ID), -0.5, 0.5))

pdf(file = paste0(Sys.Date(), "_MeMo_Granulosa_clock_predictions_Foxl2_scatter_plot.pdf"), width = 7, height = 7)
ggplot(Granulosa.predictions.Foxl2, 
       aes(x = Actual_Age + X_offset, y = Predicted_Age, color = Color)) +
  geom_point(shape = 17, size = 3, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") + 
  scale_color_identity() + 
  theme_minimal() +
  labs(title = "Predicted vs. Actual Age - Granulosa",
       x = "Actual Age", y = "Predicted Age") +
  coord_cartesian(xlim = c(0, 100), ylim = c(-20, 120)) + 
  theme(legend.position = "none")  
dev.off()

pdf(file = paste0(Sys.Date(), "_MeMo_Theca_clock_predictions_Foxl2_scatter_plot.pdf"), width = 7, height = 7)
ggplot(Theca.predictions.Foxl2, 
       aes(x = Actual_Age + X_offset, y = Predicted_Age, color = Color)) +
  geom_point(shape = 17, size = 3, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") + 
  scale_color_identity() + 
  theme_minimal() +
  labs(title = "Predicted vs. Actual Age - Theca",
       x = "Actual Age", y = "Predicted Age") +
  coord_cartesian(xlim = c(0, 100), ylim = c(-20, 120)) + 
  theme(legend.position = "none")  
dev.off()

# Assign groups based on cell ID
assign_group <- function(cell_id) {
  if (grepl("CTL_3m_30d", cell_id)) {
    return("CTL_3m_30d")
  } else if (grepl("VCD_3m_30d", cell_id)) {
    return("VCD_3m_30d")
  } else if (grepl("CTL_3m_90d", cell_id)) {
    return("CTL_3m_90d")
  } else if (grepl("VCD_3m_90d", cell_id)) {
    return("VCD_3m_90d")
  } else if (grepl("CTL_10m_30d", cell_id)) {
    return("CTL_10m_30d")
  } else if (grepl("VCD_10m_30d", cell_id)) {
    return("VCD_10m_30d")
  } else if (grepl("CTL_10m_90d", cell_id)) {
    return("CTL_10m_90d")
  } else if (grepl("VCD_10m_90d", cell_id)) {
    return("VCD_10m_90d")
  } else if (grepl("wt_young", cell_id)) {
    return("wt_young")
  } else if (grepl("het_young", cell_id)) {
    return("het_young")
  } else if (grepl("wt_midage", cell_id)) {
    return("wt_midage")
  } else if (grepl("het_midage", cell_id)) {
    return("het_midage")
  } else {
    return(NA) 
  }
}

# Assign group labels

Granulosa.predictions.combined <- Granulosa.predictions.combined %>%
  mutate(Group = sapply(Cell_ID, assign_group))

Theca.predictions.combined <- Theca.predictions.combined %>%
  mutate(Group = sapply(Cell_ID, assign_group))

Granulosa.predictions.Foxl2 <- Granulosa.predictions.Foxl2 %>%
  mutate(Group = sapply(Cell_ID, assign_group))

Theca.predictions.Foxl2 <- Theca.predictions.Foxl2 %>%
  mutate(Group = sapply(Cell_ID, assign_group))

# Perform Wilcoxon tests for each condition
perform_wilcoxon_test <- function(group1, group2, data) {
  wilcox.test(
    data$Actual_Age[data$Group == group1], 
    data$Predicted_Age[data$Group == group2]
  )$p.value
}

# Compute p-values
Granulosa.wilcoxon_results <- data.frame(
  Comparison = c("VCD_3m_30d", "VCD_3m_90d", "VCD_10m_30d", "VCD_10m_90d", 
                 "Foxl2_young", "Foxl2_midage"),
  P_Value = c(
    perform_wilcoxon_test("CTL_3m_30d", "VCD_3m_30d", Granulosa.predictions.combined),
    perform_wilcoxon_test("CTL_3m_90d", "VCD_3m_90d", Granulosa.predictions.combined),
    perform_wilcoxon_test("CTL_10m_30d", "VCD_10m_30d", Granulosa.predictions.combined),
    perform_wilcoxon_test("CTL_10m_90d", "VCD_10m_90d", Granulosa.predictions.combined),
    perform_wilcoxon_test("wt_young", "het_young", Granulosa.predictions.combined),
    perform_wilcoxon_test("wt_midage", "het_midage", Granulosa.predictions.combined)
  )
)

Theca.wilcoxon_results <- data.frame(
  Comparison = c("VCD_3m_30d", "VCD_3m_90d", "VCD_10m_30d", "VCD_10m_90d", 
                 "Foxl2_young", "Foxl2_midage"),
  P_Value = c(
    perform_wilcoxon_test("CTL_3m_30d", "VCD_3m_30d", Theca.predictions.combined),
    perform_wilcoxon_test("CTL_3m_90d", "VCD_3m_90d", Theca.predictions.combined),
    perform_wilcoxon_test("CTL_10m_30d", "VCD_10m_30d", Theca.predictions.combined),
    perform_wilcoxon_test("CTL_10m_90d", "VCD_10m_90d", Theca.predictions.combined),
    perform_wilcoxon_test("wt_young", "het_young", Theca.predictions.combined),
    perform_wilcoxon_test("wt_midage", "het_midage", Theca.predictions.combined)
  )
)

# Save results
write.table(Granulosa.wilcoxon_results, file = paste0(Sys.Date(), "_Granulosa_Actual_vs_Predicted_Age_Wilcoxon_Test_Results.txt"), 
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(Theca.wilcoxon_results, file = paste0(Sys.Date(), "_Theca_Actual_vs_Predicted_Age_Wilcoxon_Test_Results.txt"), 
            quote = FALSE, sep = "\t", row.names = FALSE)

# Calculate median
median(Granulosa.predictions.Foxl2$Predicted_Age[Granulosa.predictions.Foxl2$Group == "wt_young"])      # 27.31698 
median(Granulosa.predictions.Foxl2$Predicted_Age[Granulosa.predictions.Foxl2$Group == "het_young"])     # 43.00333
median(Granulosa.predictions.Foxl2$Predicted_Age[Granulosa.predictions.Foxl2$Group == "wt_midage"])        # 36.26867
median(Granulosa.predictions.Foxl2$Predicted_Age[Granulosa.predictions.Foxl2$Group == "het_midage"])       # 38.48362

median(Theca.predictions.Foxl2$Predicted_Age[Theca.predictions.Foxl2$Group == "wt_young"])      # 26.7499 
median(Theca.predictions.Foxl2$Predicted_Age[Theca.predictions.Foxl2$Group == "het_young"])     # 34.43742
median(Theca.predictions.Foxl2$Predicted_Age[Theca.predictions.Foxl2$Group == "wt_midage"])        # 34.88917
median(Theca.predictions.Foxl2$Predicted_Age[Theca.predictions.Foxl2$Group == "het_midage"])       # 36.33337

###############################################################
# 7. Assess age acceleration
###############################################################

# Calculate age acceleration: Predicted age - Chronological age
Granulosa.predictions.VCD$Age_accel <- Granulosa.predictions.VCD$Predicted_Age - Granulosa.predictions.VCD$Actual_Age
Theca.predictions.VCD$Age_accel <- Theca.predictions.VCD$Predicted_Age - Theca.predictions.VCD$Actual_Age

Granulosa.predictions.Foxl2$Age_accel <- Granulosa.predictions.Foxl2$Predicted_Age - Granulosa.predictions.Foxl2$Actual_Age
Theca.predictions.Foxl2$Age_accel <- Theca.predictions.Foxl2$Predicted_Age - Theca.predictions.Foxl2$Actual_Age

# Assign group labels

Granulosa.predictions.VCD <- Granulosa.predictions.VCD %>%
  mutate(Group = sapply(Cell_ID, assign_group))
Theca.predictions.VCD <- Theca.predictions.VCD %>%
  mutate(Group = sapply(Cell_ID, assign_group))

Granulosa.predictions.Foxl2 <- Granulosa.predictions.Foxl2 %>%
  mutate(Group = sapply(Cell_ID, assign_group))
Theca.predictions.Foxl2 <- Theca.predictions.Foxl2 %>%
  mutate(Group = sapply(Cell_ID, assign_group))

Granulosa.predictions.VCD$Group <- factor(Granulosa.predictions.VCD$Group, levels=c("CTL_3m_30d", "VCD_3m_30d", "CTL_3m_90d", "VCD_3m_90d", "CTL_10m_30d", "VCD_10m_30d", "CTL_10m_90d", "VCD_10m_90d"))
Theca.predictions.VCD$Group <- factor(Theca.predictions.VCD$Group, levels=c("CTL_3m_30d", "VCD_3m_30d", "CTL_3m_90d", "VCD_3m_90d", "CTL_10m_30d", "VCD_10m_30d", "CTL_10m_90d", "VCD_10m_90d"))

Granulosa.predictions.Foxl2$Group <- factor(Granulosa.predictions.Foxl2$Group, levels=c("wt_young", "het_young", "wt_midage", "het_midage"))
Theca.predictions.Foxl2$Group <- factor(Theca.predictions.Foxl2$Group, levels=c("wt_young", "het_young", "wt_midage", "het_midage"))

# Run stats
VCD.3m.30d.Granulosa.pval <- wilcox.test(Granulosa.predictions.VCD$Age_accel[Granulosa.predictions.VCD$Group == "CTL_3m_30d"], 
                                         Granulosa.predictions.VCD$Age_accel[Granulosa.predictions.VCD$Group == "VCD_3m_30d"])$p.value
VCD.3m.90d.Granulosa.pval <- wilcox.test(Granulosa.predictions.VCD$Age_accel[Granulosa.predictions.VCD$Group == "CTL_3m_90d"], 
                                         Granulosa.predictions.VCD$Age_accel[Granulosa.predictions.VCD$Group == "VCD_3m_90d"])$p.value
VCD.10m.30d.Granulosa.pval <- wilcox.test(Granulosa.predictions.VCD$Age_accel[Granulosa.predictions.VCD$Group == "CTL_10m_30d"], 
                                          Granulosa.predictions.VCD$Age_accel[Granulosa.predictions.VCD$Group == "VCD_10m_30d"])$p.value
VCD.10m.90d.Granulosa.pval <- wilcox.test(Granulosa.predictions.VCD$Age_accel[Granulosa.predictions.VCD$Group == "CTL_10m_90d"], 
                                          Granulosa.predictions.VCD$Age_accel[Granulosa.predictions.VCD$Group == "VCD_10m_90d"])$p.value

VCD.3m.30d.Theca.pval <- wilcox.test(Theca.predictions.VCD$Age_accel[Theca.predictions.VCD$Group == "CTL_3m_30d"], 
                                     Theca.predictions.VCD$Age_accel[Theca.predictions.VCD$Group == "VCD_3m_30d"])$p.value
VCD.3m.90d.Theca.pval <- wilcox.test(Theca.predictions.VCD$Age_accel[Theca.predictions.VCD$Group == "CTL_3m_90d"], 
                                     Theca.predictions.VCD$Age_accel[Theca.predictions.VCD$Group == "VCD_3m_90d"])$p.value
VCD.10m.30d.Theca.pval <- wilcox.test(Theca.predictions.VCD$Age_accel[Theca.predictions.VCD$Group == "CTL_10m_30d"], 
                                      Theca.predictions.VCD$Age_accel[Theca.predictions.VCD$Group == "VCD_10m_30d"])$p.value
VCD.10m.90d.Theca.pval <- wilcox.test(Theca.predictions.VCD$Age_accel[Theca.predictions.VCD$Group == "CTL_10m_90d"], 
                                      Theca.predictions.VCD$Age_accel[Theca.predictions.VCD$Group == "VCD_10m_90d"])$p.value

Foxl2.young.Granulosa.pval <- wilcox.test(Granulosa.predictions.Foxl2$Age_accel[Granulosa.predictions.Foxl2$Group == "wt_young"], Granulosa.predictions.Foxl2$Age_accel[Granulosa.predictions.Foxl2$Group == "het_young"])$p.value
Foxl2.old.Granulosa.pval <- wilcox.test(Granulosa.predictions.Foxl2$Age_accel[Granulosa.predictions.Foxl2$Group == "wt_midage"], Granulosa.predictions.Foxl2$Age_accel[Granulosa.predictions.Foxl2$Group == "het_midage"])$p.value

Foxl2.young.Theca.pval <- wilcox.test(Theca.predictions.Foxl2$Age_accel[Theca.predictions.Foxl2$Group == "wt_young"], Theca.predictions.Foxl2$Age_accel[Theca.predictions.Foxl2$Group == "het_young"])$p.value
Foxl2.old.Theca.pval <- wilcox.test(Theca.predictions.Foxl2$Age_accel[Theca.predictions.Foxl2$Group == "wt_midage"], Theca.predictions.Foxl2$Age_accel[Theca.predictions.Foxl2$Group == "het_midage"])$p.value

# Save values
pvals_summary <- data.frame(
  Comparison = c("Granulosa_3m_30d", "Granulosa_3m_90d", "Granulosa_10m_30d", "Granulosa_10m_90d",
                 "Theca_3m_30d", "Theca_3m_90d", "Theca_10m_30d", "Theca_10m_90d",
                 "Granulosa_wt_young", "Granulosa_wt_midage",
                 "Theca_wt_young", "Theca_wt_midage"),
  P_Value = c(VCD.3m.30d.Granulosa.pval, VCD.3m.90d.Granulosa.pval, VCD.10m.30d.Granulosa.pval, VCD.10m.90d.Granulosa.pval,
              VCD.3m.30d.Theca.pval, VCD.3m.90d.Theca.pval, VCD.10m.30d.Theca.pval, VCD.10m.90d.Theca.pval,
              Foxl2.young.Granulosa.pval, Foxl2.old.Granulosa.pval,
              Foxl2.young.Theca.pval, Foxl2.old.Theca.pval)
)

# Save the p-values as a table
write.table(pvals_summary, file = paste0(Sys.Date(), "_MeMo_transcriptional_clock_age_acceleration_Wilcoxon_pvalues.txt"), 
            quote = FALSE, sep = "\t", row.names = FALSE)

# Generate plots

pdf(file = paste0(Sys.Date(), "_MeMo_Granulosa_clock_age_acceleration_VCD_bean_plot.pdf"), width = 7, height = 7)
ggplot(Granulosa.predictions.VCD, aes(x = Group, y = Age_accel, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +  
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, color = "black") + 
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") + 
  scale_fill_manual(values = c("CTL_3m_30d" = "deeppink", "VCD_3m_30d" = "gold",
                               "CTL_3m_90d" = "deeppink", "VCD_3m_90d" = "gold",
                               "CTL_10m_30d" = "deeppink", "VCD_10m_30d" = "gold",
                               "CTL_10m_90d" = "deeppink", "VCD_10m_90d" = "gold")) + 
  coord_cartesian(ylim = c(-50, 100)) + 
  theme_minimal() +
  labs(title = "Age Acceleration by Group - Granulosa (VCD)",
       x = "Group", y = "Age Acceleration") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()


# Theca Bean Plot for VCD
pdf(file = paste0(Sys.Date(), "_MeMo_Theca_clock_age_acceleration_VCD_bean_plot.pdf"), width = 7, height = 7)
ggplot(Theca.predictions.VCD, aes(x = Group, y = Age_accel, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +  
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, color = "black") + 
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") + 
  scale_fill_manual(values = c("CTL_3m_30d" = "deeppink", "VCD_3m_30d" = "gold",
                               "CTL_3m_90d" = "deeppink", "VCD_3m_90d" = "gold",
                               "CTL_10m_30d" = "deeppink", "VCD_10m_30d" = "gold",
                               "CTL_10m_90d" = "deeppink", "VCD_10m_90d" = "gold")) + 
  coord_cartesian(ylim = c(-50, 100)) + 
  theme_minimal() +
  labs(title = "Age Acceleration by Group - Theca (VCD)",
       x = "Group", y = "Age Acceleration") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()


pdf(file = paste0(Sys.Date(), "_MeMo_Granulosa_clock_age_acceleration_Foxl2_bean_plot.pdf"), width = 7, height = 7)
ggplot(Granulosa.predictions.Foxl2, aes(x = Group, y = Age_accel, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +  
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, color = "black") + 
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") + 
  scale_fill_manual(values = c("wt_young" = "deeppink", "het_young" = "springgreen4", 
                               "wt_midage" = "deeppink", "het_midage" = "springgreen4")) + 
  coord_cartesian(ylim = c(-50, 100)) + 
  theme_minimal() +
  labs(title = "Age Acceleration by Group",
       x = "Group", y = "Age Acceleration") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

pdf(file = paste0(Sys.Date(), "_MeMo_Theca_clock_age_acceleration_Foxl2_bean_plot.pdf"), width = 7, height = 7)
ggplot(Theca.predictions.Foxl2, aes(x = Group, y = Age_accel, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +  
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, color = "black") + 
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") + 
  scale_fill_manual(values = c("wt_young" = "deeppink", "het_young" = "springgreen4", 
                               "wt_midage" = "deeppink", "het_midage" = "springgreen4")) + 
  coord_cartesian(ylim = c(-50, 100)) + 
  theme_minimal() +
  labs(title = "Age Acceleration by Group",
       x = "Group", y = "Age Acceleration") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

################################################################################
sink(file = paste(Sys.Date(), "_Transcriptional_clock_session_Info.txt"))
sessionInfo()
sink()
