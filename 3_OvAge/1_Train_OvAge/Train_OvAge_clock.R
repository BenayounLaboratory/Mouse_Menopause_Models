setwd("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/3_OvAge/")
options(stringsAsFactors = F)

# rm(list = ls())

library(dplyr)
library('beeswarm')
library('ggplot2')
library(pheatmap)
library(RColorBrewer)
library(FSA)

library(caret)
library(e1071)
library(ranger)
library(randomForest)

#######################################
# OvAge clock by weeks
# Include AMH, FSH, and INHBA to train model
# Use 75% training + 25% testing
# Use OOB error correction
#######################################

#######################################
# 1. Import data
#######################################

my.UVA.control.data <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/3_OvAge/Input/UVA_Hormone_data_for_OvAge.txt", header= TRUE, sep = "\t")

my.UVA.control.data <- my.UVA.control.data %>%
  mutate(Model = case_when(
    grepl("AC|Aging", Description, ignore.case = TRUE) ~ "Aging",
    grepl("Foxl2", Description, ignore.case = TRUE) ~ "Foxl2",
    TRUE ~ "VCD"
  ))

my.UVA.control.data.cl <- my.UVA.control.data[,c("Sample_ID", "Age_at_blood_collection_weeks", "AMH", "FSH", "INHBA")]

colnames(my.UVA.control.data.cl) <- c("Sample_ID", "Age_weeks", "AMH", "FSH", "INHBA")

# Import Fshr data - PMID: 36714222

my.Fshr.AMH.data <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/3_OvAge/Input/MeMo_Fshr_AMH_Rawdata.txt", header = TRUE)
my.Fshr.FSH.data <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/3_OvAge/Input/MeMo_Fshr_FSH_Rawdata.txt", header = TRUE)
my.Fshr.INHBA.data <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/3_OvAge/Input//MeMo_Fshr_INHA_Rawdata.txt", header = TRUE)

my.Fshr.data <- my.Fshr.FSH.data[,c("Sample_ID", "Age", "Genotype")]

my.Fshr.data$AMH <- my.Fshr.AMH.data$Concentration
my.Fshr.data$FSH <- my.Fshr.FSH.data$Concentration
my.Fshr.data$INHBA <- my.Fshr.INHBA.data$Concentration
my.Fshr.data$Model <- "Fshr"

# Calculate age by weeks

my.Fshr.data$Age_weeks <- round(as.numeric(gsub("[^0-9]", "", my.Fshr.data$Age)) * 30 / 7)
my.Fshr.data.cl <- my.Fshr.data[my.Fshr.data$Genotype == "wt",]

# head(my.Fshr.data.cl)

# Combine data
combined.data <- rbind(my.UVA.control.data.cl[,c("Sample_ID", "Age_weeks", "AMH", "FSH", "INHBA")], my.Fshr.data.cl[,c("Sample_ID", "Age_weeks", "AMH", "FSH", "INHBA")])
rownames(combined.data) <- combined.data$Sample_ID

combined.controls.with.metadata.2 <- my.UVA.control.data[,c("Sample_ID", "Age_at_blood_collection_weeks", "AMH", "FSH", "INHBA", "Model")]
colnames(combined.controls.with.metadata.2) <- c("Sample_ID", "Age_weeks", "AMH", "FSH", "INHBA", "Model")

combined.controls.with.metadata <- rbind(combined.controls.with.metadata.2, my.Fshr.data.cl[,c("Sample_ID", "Age_weeks", "AMH", "FSH", "INHBA", "Model")])
write.table(combined.controls.with.metadata, file = paste0(Sys.Date(), "_MeMo_combined_hormone_data_for_OvAge_clock_training.txt"), sep = "\t", quote = FALSE)

#######################################
# 2. Get data partition - 3/4 for training
#######################################

set.seed(1234567890)

hormone.train.idx <- createDataPartition(combined.data$Age_weeks, p=0.75, list=FALSE)

hormone.train <- combined.data[hormone.train.idx,]
dim(hormone.train) # [1] 190  11

hormone.test <- combined.data[-hormone.train.idx,]
dim(hormone.test) # [1] 62 11

hormone.test.data.with.metadata <- combined.controls.with.metadata[combined.controls.with.metadata$Sample_ID %in% hormone.test$Sample_ID,]

write.table(hormone.test.data.with.metadata, file = paste0(Sys.Date(), "_MeMo_hormone_data_for_OvAge_clock_testing_with_metadata.txt"), sep = "\t", quote = FALSE)

#######################################
# 3. Model aging clock - RF
#######################################

# Select relevant columns
data <- hormone.train[, c("Age_weeks", "AMH", "INHBA", "FSH")]

# Define predictors (X) and target variable (y)
X <- data[, c("AMH", "INHBA", "FSH")]
y <- data$Age_weeks

# OOB tuning for mtry
oob_control <- trainControl(method = "oob")

oob_rf <- train(
  x = X, y = y,
  method = "rf",
  trControl = oob_control,
  tuneGrid = expand.grid(mtry = c(1, 2, 3)), 
  ntree = 500,
  importance = TRUE
)

# Extract OOB results
best_mtry_oob <- oob_rf$bestTune$mtry             # 2
oob_rmse <- min(oob_rf$results$RMSE)              # 12.30601
oob_r2 <- max(oob_rf$results$Rsquared)            # 0.6745342

# Plot OOB results
plot(oob_rf)

#######################################
# 4. Test model
#######################################

pred.RF.hormone   <- predict(oob_rf, hormone.test[,3:5])
hormone.test$pred.RF   <- pred.RF.hormone

cor.test(hormone.test$Age_weeks, hormone.test$pred.RF, method = 'spearman')      # p-value = 7.352e-08, rho = 0.6207666

# Save clock

OvAge <- oob_rf
save(OvAge, file = paste0(Sys.Date(),"_MeMo_UVA_hormone_OvAge_clock_randomforest.RData"))

# Save OOB prediction

oob_predictions <- OvAge$finalModel$predicted

oob_df <- data.frame(OOB_Prediction = oob_predictions)

rownames(oob_df) <- gsub("^X", "", rownames(oob_df))

oob_df$Sample_ID <- rownames(oob_df)

combined.data$OOB_prediction <- NA

combined.data$OOB_prediction[match(rownames(oob_df), rownames(combined.data))] <- oob_df$OOB_Prediction

write.table(combined.data, file = paste0(Sys.Date(), "_MeMo_OvAge_OOB_prediction_table.txt"), sep = "\t", quote = FALSE)

#######################################
# 5. Generate performance plot
#######################################

pdf(paste(Sys.Date(), "RF_OvAge_clock_performance_actual_vs_pred_age_without_OOB_predictions.pdf", sep = "_"), width = 8, height = 7)
plot(hormone.test$Age_weeks, hormone.test$pred.RF, 
     xlab = "Actual age (weeks)", ylab = "Predicted age (weeks)",
     main = "OvAge clock - test data (n=62)",
     col = "black", cex = 1, pch = 19, xlim = c(0, 100), ylim = c(0, 100),
     cex.lab = 1.5, cex.axis = 1.5)
abline(0, 1, col = "red", lty = "dashed")
text(20, 90, "Rho ~0.6207666", cex = 1.5, col = "blue")
text(20, 85, "pval ~7.352e-08", cex = 1.5, col = "blue")
dev.off()

#######################
sink(file = paste0(Sys.Date(), "_OvAge_clock_session_info.txt"))
sessionInfo()
sink()
