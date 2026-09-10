options(stringsAsFactors = FALSE)

library(caret)
library(randomForest)
library(cocor)
library(ggplot2)

rm(list = ls())

#######################################
# Menopause-model project / revision
# Benchmark: log-linear AMH-only model vs OvAge clock (AMH+FSH+INHBA RF) for
# predicting chronological age
#######################################

#######################################
# 1. Import data
#######################################

# Pre-trained OvAge model + its saved train/test partition 

ovage.dir <- "~/Dropbox/Benayoun_lab/Menopause_model_project/OvAge/2026-02-18"

load(file.path(ovage.dir, "2026-02-18_MeMo_UVA_hormone_OvAge_clock_randomforest.RData"))   

hormone.train <- read.table(file.path(ovage.dir, "2026-02-18_MeMo_OvAge_train_set_hormone_list_with_prediction.txt"),
                             header = TRUE, sep = "\t")
hormone.test  <- read.table(file.path(ovage.dir, "2026-02-18_MeMo_OvAge_test_set_hormone_list_with_prediction.txt"),
                             header = TRUE, sep = "\t")

#######################################
# 2. Fit AMH-only log-linear model
#######################################

amh.log <- lm(Age_weeks ~ log(AMH), data = hormone.train)

#######################################
# 3. Chronological age prediction (test set, n=68)
#######################################

hormone.test$pred.AMH.log <- predict(amh.log, newdata = hormone.test)

cor.AMHlog.age <- cor.test(hormone.test$Age_weeks, hormone.test$pred.AMH.log, method = "spearman", exact = FALSE)
cor.RF.age     <- cor.test(hormone.test$Age_weeks, hormone.test$pred.RF,      method = "spearman", exact = FALSE)

# Calculate RMSE
rmse <- function(actual, pred) sqrt(mean((actual - pred)^2))

rmse.AMHlog.age <- rmse(hormone.test$Age_weeks, hormone.test$pred.AMH.log)
rmse.RF.age     <- rmse(hormone.test$Age_weeks, hormone.test$pred.RF)

# Assess correlation
r.kh.age.log <- cor(hormone.test$pred.AMH.log, hormone.test$pred.RF, method = "spearman")
cocor.age.log <- cocor.dep.groups.overlap(r.jk = as.numeric(cor.AMHlog.age$estimate),
                                           r.jh = as.numeric(cor.RF.age$estimate),
                                           r.kh = r.kh.age.log,
                                           n = nrow(hormone.test), test = "steiger1980")

steiger.table <- data.frame(
  Comparison = "AMH-only (log-linear) vs OvAge",
  Outcome    = "Chronological age",
  r_jk       = cocor.age.log@r.jk,
  r_jh       = cocor.age.log@r.jh,
  r_kh       = cocor.age.log@r.kh,
  diff       = cocor.age.log@diff,
  n          = cocor.age.log@n,
  z          = cocor.age.log@steiger1980$statistic,
  p_value    = cocor.age.log@steiger1980$p.value
)
steiger.table[, c("r_jk", "r_jh", "r_kh", "diff", "z")] <- round(steiger.table[, c("r_jk", "r_jh", "r_kh", "diff", "z")], 4)

write.table(steiger.table, paste0(Sys.Date(), "_MeMo_OvAge_vs_AMH_steiger_test_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

#######################################
# 4. Summary table
#######################################

summary.table <- data.frame(
  Predictor  = c("AMH-only (log-linear)", "OvAge"),
  Outcome    = c("Chronological age", "Chronological age"),
  Rho        = c(as.numeric(cor.AMHlog.age$estimate), as.numeric(cor.RF.age$estimate)),
  p_value    = c(cor.AMHlog.age$p.value, cor.RF.age$p.value),
  RMSE_weeks = c(rmse.AMHlog.age, rmse.RF.age),
  n          = c(nrow(hormone.test), nrow(hormone.test))
)

write.table(summary.table, paste0(Sys.Date(), "_MeMo_OvAge_vs_AMH_summary_table.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

save(hormone.test, hormone.train, summary.table, amh.log, cocor.age.log, boot.rmse.diff.age.log,
     file = paste0(Sys.Date(), "_MeMo_OvAge_vs_AMH_benchmark_results.RData"))

#######################################
# 5. Generate plot
#######################################

pdf(paste0(Sys.Date(), "_MeMo_AMH_only_LOG_LINEAR_clock_performance_actual_vs_pred_age.pdf"), width = 8, height = 7)
plot(hormone.test$Age_weeks, hormone.test$pred.AMH.log,
     xlab = "Actual age (weeks)", ylab = "Predicted age (weeks)",
     main = "AMH-only clock [log-linear] - test data (n=68)",
     col = "black", cex = 1, pch = 19, xlim = c(0, 100), ylim = c(0, 100),
     cex.lab = 1.5, cex.axis = 1.5)
abline(0, 1, col = "red", lty = "dashed")
text(20, 90, paste0("Rho ~", round(cor.AMHlog.age$estimate, 4)), cex = 1.5, col = "blue")
text(20, 85, paste0("pval ~", format(cor.AMHlog.age$p.value, digits = 4)), cex = 1.5, col = "blue")
dev.off()

#######################################

sink(file = paste0(Sys.Date(), "_3_AMH_only_model_session_info.txt"))
sessionInfo()
sink()
