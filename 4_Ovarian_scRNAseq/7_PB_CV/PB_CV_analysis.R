setwd("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/For_github/PB_CV")
options(stringsAsFactors = FALSE)

# Load libraries
library(ggplot2)
library(dplyr)
library(matrixStats)
library(data.table)

rm(list = ls())

theme_set(theme_bw())
set.seed(123)

################################################################################
# PB CV analysis
# Aging, VCD and Foxl2 models
################################################################################

# Define global variables
cell_type_order <- c("Granulosa", "Theca", "Stroma", "Mesenchyme", "Pericyte",
                     "BEC", "LEC", "Epithelial", "Neutrophil", "Monocyte",  
                     "Macrophage", "DC", "ILC", "NK", "NKT", "CD8 NKT",  
                     "CD8 T", "CD4 T", "DNT", "DPT", "B")

bubble_data_all <- list()

################################################################################
# 1. Aging Model: YF vs OF
################################################################################

load("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/For_github/PB_CV/Input/2025-01-14_10x_ovary_Benayoun_lab_AC_PB_VST_counts.RData")

for (cell_type in names(vst.counts)) {
  df <- vst.counts[[cell_type]]
  YF <- df[, grep("YF", colnames(df))]
  OF <- df[, grep("OF", colnames(df))]
  
  cv_YF <- apply(YF, 1, function(x) sd(x) / mean(x))
  cv_OF <- apply(OF, 1, function(x) sd(x) / mean(x))
  
  log2FC <- log2(median(cv_OF, na.rm = TRUE) / median(cv_YF, na.rm = TRUE))
  p <- wilcox.test(cv_YF, cv_OF)$p.value
  
  bubble_data_all[[paste0("Aging_", cell_type)]] <- data.frame(
    Model = "Aging",
    CellType = cell_type,
    Comparison = "YF vs OF",
    log2FC_CV = log2FC,
    PValue = p
  )
}

################################################################################
# 2. VCD Model: VCD vs CTL
################################################################################

load("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/For_github/PB_CV/Input/2025-01-15_10x_ovary_Benayoun_lab_VCD_PB_VST_counts.RData")

vcd_groups <- list(
  "CTL_3m_30d" = "CTL_3m_30d", "VCD_3m_30d" = "VCD_3m_30d",
  "CTL_3m_90d" = "CTL_3m_90d", "VCD_3m_90d" = "VCD_3m_90d",
  "CTL_10m_30d" = "CTL_10m_30d", "VCD_10m_30d" = "VCD_10m_30d",
  "CTL_10m_90d" = "CTL_10m_90d", "VCD_10m_90d" = "VCD_10m_90d"
)

for (cell_type in names(vst.counts)) {
  df <- vst.counts[[cell_type]]
  
  for (grp in c("3m_30d", "3m_90d", "10m_30d", "10m_90d")) {
    ctl <- df[, grep(paste0("CTL_", grp), colnames(df))]
    vcd <- df[, grep(paste0("VCD_", grp), colnames(df))]
    
    if (ncol(ctl) > 0 && ncol(vcd) > 0) {
      cv_ctl <- apply(ctl, 1, function(x) sd(x) / mean(x))
      cv_vcd <- apply(vcd, 1, function(x) sd(x) / mean(x))
      
      log2FC <- log2(median(cv_vcd, na.rm = TRUE) / median(cv_ctl, na.rm = TRUE))
      p <- wilcox.test(cv_ctl, cv_vcd)$p.value
      
      bubble_data_all[[paste0("VCD_", cell_type, "_", grp)]] <- data.frame(
        Model = "VCD",
        CellType = cell_type,
        Comparison = paste0("VCD vs CTL (", grp, ")"),
        log2FC_CV = log2FC,
        PValue = p
      )
    }
  }
}

################################################################################
# 3. Foxl2 Model: het vs wt (young & old)
################################################################################

load("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/For_github/PB_CV/Input/2025-01-15_10x_ovary_Benayoun_lab_Foxl2_NULL_PB_VST_counts.RData")

for (cell_type in names(vst.counts)) {
  df <- vst.counts[[cell_type]]
  
  y_wt <- df[, grep("Foxl2_wt_young", colnames(df))]
  y_het <- df[, grep("Foxl2_het_young", colnames(df))]
  o_wt <- df[, grep("Foxl2_wt_old", colnames(df))]
  o_het <- df[, grep("Foxl2_het_old", colnames(df))]
  
  # Young
  if (ncol(y_wt) > 0 && ncol(y_het) > 0) {
    cv_wt <- apply(y_wt, 1, function(x) sd(x) / mean(x))
    cv_het <- apply(y_het, 1, function(x) sd(x) / mean(x))
    log2FC <- log2(median(cv_het, na.rm = TRUE) / median(cv_wt, na.rm = TRUE))
    p <- wilcox.test(cv_wt, cv_het)$p.value
    
    bubble_data_all[[paste0("Foxl2_young_", cell_type)]] <- data.frame(
      Model = "Foxl2",
      CellType = cell_type,
      Comparison = "young_het vs young_wt",
      log2FC_CV = log2FC,
      PValue = p
    )
  }
  
  # Old
  if (ncol(o_wt) > 0 && ncol(o_het) > 0) {
    cv_wt <- apply(o_wt, 1, function(x) sd(x) / mean(x))
    cv_het <- apply(o_het, 1, function(x) sd(x) / mean(x))
    log2FC <- log2(median(cv_het, na.rm = TRUE) / median(cv_wt, na.rm = TRUE))
    p <- wilcox.test(cv_wt, cv_het)$p.value
    
    bubble_data_all[[paste0("Foxl2_old_", cell_type)]] <- data.frame(
      Model = "Foxl2",
      CellType = cell_type,
      Comparison = "old_het vs old_wt",
      log2FC_CV = log2FC,
      PValue = p
    )
  }
}

################################################################################
# 4. Combine and Plot Bubble Plot
################################################################################

bubble_df <- bind_rows(bubble_data_all)

# Adjust p-values and add significance info
bubble_df$FDR <- p.adjust(bubble_df$PValue, method = "fdr")
bubble_df$log10P <- -log10(pmax(bubble_df$FDR, 1e-300))
bubble_df$CellType <- factor(bubble_df$CellType, levels = cell_type_order)
bubble_df$Model <- factor(bubble_df$Model, levels = c("Aging", "VCD", "Foxl2"))

# Bubble plot
pdf(paste0(Sys.Date(), "_CV_Combined_BubblePlot.pdf"), width = 16, height = 6)
ggplot(bubble_df, aes(x = CellType, y = Comparison, size = log10P, color = log2FC_CV)) +
  geom_point(alpha = 0.85) +
  facet_wrap(~Model, nrow = 1, scales = "free_y") +
  scale_size_continuous(range = c(2, 10)) +
  scale_color_gradientn(
    colors = c("darkblue", "dodgerblue4", "dodgerblue3", "dodgerblue1",
               "white", "lightcoral", "brown1", "firebrick2", "firebrick4"),
    limits = c(-1.5, 1.5),
    oob = scales::squish,
    na.value = "gray90"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Coefficient of Variation (CV) Differences Across Models",
    x = "Cell Type", y = "Comparison",
    size = "-log10(FDR)", color = "log2FC (CV)"
  )
dev.off()

################################################################################
sink(file = paste0(Sys.Date(), "_PB_CV_combined_sessionInfo.txt"))
sessionInfo()
sink()
