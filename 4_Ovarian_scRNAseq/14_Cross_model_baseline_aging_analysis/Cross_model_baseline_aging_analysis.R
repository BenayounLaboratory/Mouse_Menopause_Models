library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(DESeq2)

rm(list = ls())
set.seed(123123)

################################################################################
# Menopause-model project / revision
# Within-model baseline aging-rate correlation across Aging, VCD-vehicle, and Foxl2-WT arms
# Each model's own control population is regressed on chronological age, and then vectors are correlated pairwise (Spearman).
################################################################################

BASE <- "/Volumes/jinho01/Benayoun_lab/Projects/Menopause_model_project/Data/Ovarian_scRNAseq/PB_DESeq2/Without_scTE"

SHARED_CELLTYPES <- c("Granulosa", "Theca", "Stroma", "BEC", "Epithelial", "DNT", "LEC")

MIN_N_GENES <- 500
FDR_CUTOFF  <- 0.05

################################################################################
# 1. Load data
################################################################################

# Aging model: young (YF) vs old (OF) DESeq2 results, rescaled to log2FC/month below
AC.deseq.res.list <- get(load(file.path(BASE, "AC/2026-02-08/DESeq2_results/2026-02-08_10x_ovary_Benayoun_lab_AC_PB_DESeq2_results.RData")))

# Raw pseudobulk counts (pre-VST, pre-SVA)
VCD.counts.pb <- get(load(file.path(BASE,
  "VCD/2026-02-08/2026-02-08_MeMo_VCD_cohort_PB_counts_post_QC.RData")))
Foxl2.counts.pb <- get(load(file.path(BASE,
  "Foxl2/2026-02-08/2026-02-08_MeMo_Foxl2_NULL_cohort_PB_counts_post_QC.RData")))

################################################################################
# 2. Helper functions
################################################################################

# Named log2FC vector
get_l2fc <- function(deseq_list, celltype) {
  res <- as.data.frame(deseq_list[[celltype]])
  v <- res$log2FoldChange
  names(v) <- rownames(res)
  v
}

# Significant genes (FDR < cutoff) for one cell type from a DESeq2 results list
get_sig_genes <- function(deseq_list, celltype, cutoff = FDR_CUTOFF) {
  res <- as.data.frame(deseq_list[[celltype]])
  rownames(res)[!is.na(res$padj) & res$padj < cutoff]
}

# Spearman correlation between two named vectors, restricted to shared names
spearman_pair <- function(v1, v2, restrict_to = NULL) {
  common <- intersect(names(v1), names(v2))
  if (!is.null(restrict_to)) common <- intersect(common, restrict_to)
  if (length(common) < 3) {
    return(list(rho = NA_real_, p = NA_real_, n = length(common)))
  }
  ct <- suppressWarnings(cor.test(v1[common], v2[common], method = "spearman"))
  list(rho = unname(ct$estimate), p = ct$p.value, n = length(common))
}

# VCD library ID (e.g. "CTL_10m_30d_1") -> age at collection (months)
vcd_age_months <- function(sample_id) {
  m <- regmatches(sample_id, regexec("^(CTL|VCD)_([0-9]+)m_([0-9]+)d_([0-9]+)$", sample_id))
  age_inj_m <- as.numeric(vapply(m, `[`, character(1), 3))
  dur_d     <- as.numeric(vapply(m, `[`, character(1), 4))
  age_inj_m + dur_d / 30.44
}
vcd_is_ctl <- function(sample_id) grepl("^CTL_", sample_id)

# Foxl2 library ID (e.g. "Foxl2_wt_young_1") -> age (months); matches the
# 'dataDesign' logic in the Foxl2 DESeq2 pipeline script
foxl2_age_lookup <- c(young = 4, old = 9, supold = 17)
foxl2_age_months <- function(sample_id) {
  m <- regmatches(sample_id, regexec("^Foxl2_(het|wt)_(young|old|supold)_([0-9]+)$", sample_id))
  grp <- vapply(m, `[`, character(1), 3)
  unname(foxl2_age_lookup[grp])
}
foxl2_is_wt <- function(sample_id) grepl("^Foxl2_wt_", sample_id)

# Fit expression ~ age_months per gene 
fit_age_slope_deseq2 <- function(counts_mat, age_months, min_frac_expressed) {
  good_genes <- rowSums(counts_mat > 0) >= ceiling(min_frac_expressed * ncol(counts_mat))
  counts_mat <- counts_mat[good_genes, , drop = FALSE]

  col_data <- data.frame(row.names = colnames(counts_mat), age_months = age_months)

  dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = col_data, design = ~ age_months)
  dds <- DESeq(dds, quiet = TRUE)

  res <- as.data.frame(results(dds, name = "age_months"))
  res <- res[!is.na(res$padj), ]

  data.frame(
    gene   = rownames(res),
    slope  = res$log2FoldChange,
    pvalue = res$pvalue,
    FDR    = res$padj,
    row.names = rownames(res)
  )
}

################################################################################
# 3. Baseline aging-rate (log2FC/month) correlation
################################################################################

# Aging model
AGING_DELTA_MONTHS <- 16

aging_slope_list <- setNames(
  lapply(SHARED_CELLTYPES, function(ct) get_l2fc(AC.deseq.res.list, ct) / AGING_DELTA_MONTHS),
  SHARED_CELLTYPES
)

# VCD-vehicle (CTL-only) model
vcd_slope_list <- list()
vcd_ctl_metadata <- list()

for (ct in SHARED_CELLTYPES) {
  counts_mat <- VCD.counts.pb[[ct]]
  ctl_cols   <- colnames(counts_mat)[vcd_is_ctl(colnames(counts_mat))]
  ctl_counts <- counts_mat[, ctl_cols, drop = FALSE]
  age_m      <- vcd_age_months(ctl_cols)

  fit <- fit_age_slope_deseq2(ctl_counts, age_m, min_frac_expressed = 14/16)

  v <- fit$slope; names(v) <- fit$gene
  vcd_slope_list[[ct]] <- v

  vcd_ctl_metadata[[ct]] <- tibble(
    CellType = ct, n_libraries = length(ctl_cols),
    n_ages = length(unique(round(age_m, 2))),
    age_months_min = round(min(age_m), 2), age_months_max = round(max(age_m), 2),
    age_months_values = paste(sort(unique(round(age_m, 2))), collapse = ", "),
    n_genes_fit = nrow(fit)
  )
  assign(paste0("vcd_sig_", ct), fit$gene[!is.na(fit$FDR) & fit$FDR < FDR_CUTOFF])
}
vcd_ctl_metadata <- bind_rows(vcd_ctl_metadata)

# Foxl2 haploinsufficiency model
foxl2_slope_list <- list()
foxl2_wt_metadata <- list()

for (ct in SHARED_CELLTYPES) {
  counts_mat <- Foxl2.counts.pb[[ct]]
  wt_cols    <- colnames(counts_mat)[foxl2_is_wt(colnames(counts_mat))]
  wt_counts  <- counts_mat[, wt_cols, drop = FALSE]
  age_m      <- foxl2_age_months(wt_cols)

  fit <- fit_age_slope_deseq2(wt_counts, age_m, min_frac_expressed = 12/16)

  v <- fit$slope; names(v) <- fit$gene
  foxl2_slope_list[[ct]] <- v

  age_tab <- table(age_m)
  foxl2_wt_metadata[[ct]] <- tibble(
    CellType = ct, n_libraries = length(wt_cols),
    n_ages = length(age_tab),
    age_months_values = paste(sort(unique(age_m)), collapse = ", "),
    n_per_age_group = paste(names(age_tab), "mo: n=", as.integer(age_tab), collapse = "; "),
    n_in_oldest_group = as.integer(age_tab[as.character(max(age_m))]),
    leverage_flag = as.integer(age_tab[as.character(max(age_m))]) == 1,
    n_genes_fit = nrow(fit)
  )
  assign(paste0("foxl2_sig_", ct), fit$gene[!is.na(fit$FDR) & fit$FDR < FDR_CUTOFF])
}
foxl2_wt_metadata <- bind_rows(foxl2_wt_metadata)

write.csv(vcd_ctl_metadata,   file = paste0(Sys.Date(), "_VCD_CTL_control_library_metadata.csv"), row.names = FALSE)
write.csv(foxl2_wt_metadata,  file = paste0(Sys.Date(), "_Foxl2_WT_control_library_metadata.csv"), row.names = FALSE)

# Perform pairwise Spearman correlation of log2FC/month vectors

slope_lists <- list(Aging = aging_slope_list, VCD = vcd_slope_list, Foxl2 = foxl2_slope_list)

sig_genes_by_model <- function(model, ct) {
  if (model == "Aging") get_sig_genes(AC.deseq.res.list, ct)
  else if (model == "VCD") get(paste0("vcd_sig_", ct))
  else if (model == "Foxl2") get(paste0("foxl2_sig_", ct))
}

model_pairs <- list(
  c("Aging", "VCD"),
  c("Aging", "Foxl2"),
  c("VCD", "Foxl2")
)

results <- list()
i <- 1
for (celltype in SHARED_CELLTYPES) {
  for (pr in model_pairs) {
    m1 <- pr[1]; m2 <- pr[2]
    v1 <- slope_lists[[m1]][[celltype]]
    v2 <- slope_lists[[m2]][[celltype]]

    primary <- spearman_pair(v1, v2)

    sig_union <- union(sig_genes_by_model(m1, celltype), sig_genes_by_model(m2, celltype))
    secondary <- spearman_pair(v1, v2, restrict_to = sig_union)

    label_m1 <- ifelse(m1 == "VCD", "VCD-vehicle", ifelse(m1 == "Foxl2", "Foxl2-WT", m1))
    label_m2 <- ifelse(m2 == "VCD", "VCD-vehicle", ifelse(m2 == "Foxl2", "Foxl2-WT", m2))

    results[[i]] <- tibble(
      CellType   = celltype,
      ModelPair  = paste0(label_m1, "_vs_", label_m2),
      Rho        = primary$rho,
      Pvalue     = primary$p,
      n          = primary$n,
      Rho_sigOnly    = secondary$rho,
      Pvalue_sigOnly = secondary$p,
      n_sigOnly      = secondary$n
    )
    i <- i + 1
  }
}
results <- bind_rows(results) %>%
  mutate(
    FDR      = p.adjust(Pvalue, method = "BH"),
    LowN_flag = n < MIN_N_GENES
  ) %>%
  relocate(FDR, .after = Pvalue) %>%
  relocate(LowN_flag, .after = n)

write.csv(results, file = paste0(Sys.Date(), "_baseline_aging_rate_correlation.csv"), row.names = FALSE)

################################################################################
sink(file = paste0(Sys.Date(), "_Cross_model_baseline_aging_analysis_session_info.txt"))
sessionInfo()
sink()
