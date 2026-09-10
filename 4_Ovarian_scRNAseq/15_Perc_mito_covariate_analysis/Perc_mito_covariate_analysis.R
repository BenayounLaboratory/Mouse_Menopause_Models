# Load libraries
library(dplyr)
library(Seurat)
library(DESeq2)
library(sva)
library(limma)

rm(list = ls())

################################################################################
# Menopause-model project / revision
# %mito covariate-adjustment sensitivity
################################################################################

ALPHA <- 0.05
SHARED_CELLTYPES <- c("Granulosa", "Theca", "Stroma", "BEC", "Epithelial", "DNT")

################################################################################
# 1. Helper functions
################################################################################

# %mito computed directly from pseudobulk counts
compute_pb_mito_pct <- function(counts_mat, pattern = "^mt-") {
  mt_genes <- grep(pattern, rownames(counts_mat), value = TRUE)
  100 * colSums(counts_mat[mt_genes, , drop = FALSE]) / colSums(counts_mat)
}

# SVA + DESeq2 with mito covariate
run_mito_only_deseq2 <- function(counts_full, dataDesign_base, batch_vec,
                                  sva_formula, sva_keep_cols,
                                  final_formula_base, contrast,
                                  min_gene_samples, covar_df) {

  good_genes <- rowSums(counts_full > 0) >= min_gene_samples
  counts <- counts_full[good_genes, , drop = FALSE]

  dataDesign <- dataDesign_base
  dataDesign$batch <- batch_vec
  m <- match(rownames(dataDesign), covar_df$sample_id)
  dataDesign$pct_mito_pb   <- covar_df$pct_mito_pb[m]
  dataDesign$pct_mito_pb_c <- as.numeric(scale(dataDesign$pct_mito_pb, scale = FALSE))

  n_samples <- ncol(counts)

  mod1 <- model.matrix(sva_formula, data = dataDesign)
  set.seed(123123)
  n.sv <- tryCatch(num.sv(counts, mod1, method = "be"), error = function(e) 0)
  my.svseq <- tryCatch(svaseq(as.matrix(counts), mod1, n.sv = n.sv, constant = 0.1),
                        error = function(e) NULL)

  if (!is.null(my.svseq) && !is.null(my.svseq$n.sv) && my.svseq$n.sv > 0) {
    my.clean <- removeBatchEffect(log2(counts + 0.1), batch = dataDesign$batch,
                                   covariates = my.svseq$sv, design = mod1[, sva_keep_cols, drop = FALSE])
  } else {
    my.clean <- removeBatchEffect(log2(counts + 0.1), batch = dataDesign$batch,
                                   design = mod1[, sva_keep_cols, drop = FALSE])
  }
  sva.cleaned <- round(2 ^ my.clean - 0.1)
  sva.cleaned[sva.cleaned < 0] <- 0
  mode(sva.cleaned) <- "integer"

  final_formula <- as.formula(paste(paste(deparse(final_formula_base), collapse = ""), "+ pct_mito_pb_c"))
  design_mat  <- model.matrix(final_formula, data = dataDesign)
  design_rank <- qr(design_mat)$rank
  resid_df    <- n_samples - design_rank

  status <- data.frame(
    n_samples = n_samples, n_genes_tested = nrow(counts),
    design_terms = design_rank, resid_df = resid_df,
    flag = if (resid_df < 1) "SKIPPED_INSUFFICIENT_DF"
           else if (resid_df <= 3) "THIN_MARGIN_DF"
           else "OK"
  )
  if (resid_df < 1) return(list(status = status, res = NULL))

  dds <- DESeqDataSetFromMatrix(countData = sva.cleaned, colData = dataDesign, design = final_formula)
  dds <- tryCatch(DESeq(dds), error = function(e) { status$flag <<- paste0("DESEQ_ERROR: ", conditionMessage(e)); NULL })
  if (is.null(dds)) return(list(status = status, res = NULL))

  res <- results(dds, contrast = contrast)
  res <- res[!is.na(res$padj), ]
  list(status = status, res = res)
}

all_covariates <- list()
deseq.res.list.mito.only <- list()
all_status <- list()

################################################################################
# 2. Aging model
################################################################################

load("~/Benayoun_lab/Projects/Menopause_model_project/Data/Ovarian_scRNAseq/PB_DESeq2/Without_scTE/AC/2026-02-08/2026-02-08_MeMo_Aging_cohort_PB_counts_post_QC.RData")
counts.pb.AC <- counts.pb
rm(counts.pb); gc()

covariates_AC <- do.call(rbind, lapply(SHARED_CELLTYPES, function(ct) {
  mat <- counts.pb.AC[[ct]]
  data.frame(dataset = "AC", celltype = ct, sample_id = colnames(mat),
             pct_mito_pb = compute_pb_mito_pct(mat), row.names = NULL)
}))
all_covariates[["AC"]] <- covariates_AC

for (ct in SHARED_CELLTYPES) {
  cat("  [AC]", ct, "\n")
  counts <- counts.pb.AC[[ct]]
  dataDesign_base <- data.frame(row.names = colnames(counts),
                                 age = ifelse(grepl("YF", colnames(counts)), "YF", "OF"))
  batch_vec <- ifelse(grepl("3", colnames(counts)), "AC_2", "AC_1")
  covar_ct <- covariates_AC[covariates_AC$celltype == ct, ]

  out <- run_mito_only_deseq2(counts, dataDesign_base, batch_vec,
                               sva_formula = ~ age + batch, sva_keep_cols = 1:2,
                               final_formula_base = ~ age, contrast = c("age", "OF", "YF"),
                               min_gene_samples = 4, covar_df = covar_ct)
  out$status$dataset <- "AC"; out$status$celltype <- ct
  all_status[[paste("AC", ct)]] <- out$status
  if (!is.null(out$res)) deseq.res.list.mito.only[["AC"]][[ct]] <- out$res
}
rm(counts.pb.AC); gc()

################################################################################
# 3. VCD model
################################################################################

load("~/Benayoun_lab/Projects/Menopause_model_project/Data/Ovarian_scRNAseq/PB_DESeq2/Without_scTE/VCD/2026-02-08/2026-02-08_MeMo_VCD_cohort_PB_counts_post_QC.RData")
counts.pb.VCD <- counts.pb
rm(counts.pb); gc()

load("~/Benayoun_lab/Projects/Menopause_model_project/Data/Ovarian_scRNAseq/0_Annotated_Seurat_objects/Without_scTE/2024-10-24_10x_ovary_Benayoun_lab_VCD_Seurat_object_with_final_annotation.RData")
md_VCD <- ovary.VCD@meta.data
rm(ovary.VCD); gc()

vcd_metadata <- md_VCD %>% dplyr::select(Library, Age, Treatment, Duration, Batch) %>% distinct()
rownames(vcd_metadata) <- vcd_metadata$Library
rm(md_VCD); gc()

covariates_VCD <- do.call(rbind, lapply(SHARED_CELLTYPES, function(ct) {
  mat <- counts.pb.VCD[[ct]]
  data.frame(dataset = "VCD", celltype = ct, sample_id = colnames(mat),
             pct_mito_pb = compute_pb_mito_pct(mat), row.names = NULL)
}))
all_covariates[["VCD"]] <- covariates_VCD

for (ct in SHARED_CELLTYPES) {
  cat("  [VCD]", ct, "\n")
  counts <- counts.pb.VCD[[ct]]
  dataDesign_base <- data.frame(
    row.names = colnames(counts),
    age       = ifelse(grepl("3m", colnames(counts)), "3m", "10m"),
    treatment = ifelse(grepl("CTL", colnames(counts)), "CTL", "VCD"),
    duration  = ifelse(grepl("30d", colnames(counts)), "30d", "90d")
  )
  batch_vec <- vcd_metadata[colnames(counts), "Batch"]
  covar_ct <- covariates_VCD[covariates_VCD$celltype == ct, ]

  out <- run_mito_only_deseq2(counts, dataDesign_base, batch_vec,
                               sva_formula = ~ age + treatment + duration + batch, sva_keep_cols = 1:4,
                               final_formula_base = ~ age + treatment + duration, contrast = c("treatment", "VCD", "CTL"),
                               min_gene_samples = 14, covar_df = covar_ct)
  out$status$dataset <- "VCD"; out$status$celltype <- ct
  all_status[[paste("VCD", ct)]] <- out$status
  if (!is.null(out$res)) deseq.res.list.mito.only[["VCD"]][[ct]] <- out$res
}
rm(counts.pb.VCD); gc()

################################################################################
# 4. Foxl2 model
################################################################################

load("~/Benayoun_lab/Projects/Menopause_model_project/Data/Ovarian_scRNAseq/PB_DESeq2/Without_scTE/Foxl2/2026-02-08/2026-02-08_MeMo_Foxl2_NULL_cohort_PB_counts_post_QC.RData")
counts.pb.Foxl2 <- counts.pb
rm(counts.pb); gc()

load("~/Benayoun_lab/Projects/Menopause_model_project/Data/Ovarian_scRNAseq/0_Annotated_Seurat_objects/Without_scTE/2025-11-05_10x_ovary_Benayoun_lab_Foxl2_Seurat_object_with_final_annotation.RData")
md_Foxl2 <- ovary.Foxl2@meta.data
rm(ovary.Foxl2); gc()

foxl2_metadata <- md_Foxl2 %>% dplyr::select(Library, Age, Genotype, Batch) %>% distinct()
rownames(foxl2_metadata) <- foxl2_metadata$Library
rm(md_Foxl2); gc()

covariates_Foxl2 <- do.call(rbind, lapply(SHARED_CELLTYPES, function(ct) {
  mat <- counts.pb.Foxl2[[ct]]
  data.frame(dataset = "Foxl2", celltype = ct, sample_id = colnames(mat),
             pct_mito_pb = compute_pb_mito_pct(mat), row.names = NULL)
}))
all_covariates[["Foxl2"]] <- covariates_Foxl2

for (ct in SHARED_CELLTYPES) {
  cat("  [Foxl2]", ct, "\n")
  counts <- counts.pb.Foxl2[[ct]]
  dataDesign_base <- data.frame(
    row.names = colnames(counts),
    age      = ifelse(grepl("young", colnames(counts)), "4m",
                       ifelse(grepl("supold", colnames(counts)), "17m", "9m")),
    genotype = ifelse(grepl("_wt_", colnames(counts)), "wt", "het")
  )
  batch_vec <- foxl2_metadata[colnames(counts), "Batch"]
  covar_ct <- covariates_Foxl2[covariates_Foxl2$celltype == ct, ]

  out <- run_mito_only_deseq2(counts, dataDesign_base, batch_vec,
                               sva_formula = ~ age + genotype + batch, sva_keep_cols = 1:4,
                               final_formula_base = ~ genotype + age, contrast = c("genotype", "het", "wt"),
                               min_gene_samples = 12, covar_df = covar_ct)
  out$status$dataset <- "Foxl2"; out$status$celltype <- ct
  all_status[[paste("Foxl2", ct)]] <- out$status
  if (!is.null(out$res)) deseq.res.list.mito.only[["Foxl2"]][[ct]] <- out$res
}
rm(counts.pb.Foxl2); gc()

################################################################################
# 5. Assemble %mito covariate table; save mito DESeq2 results
################################################################################

covariate_table <- do.call(rbind, all_covariates)
write.table(covariate_table,
            file = paste0(Sys.Date(), "_MeMo_AC_VCD_Foxl2_pseudobulk_mito_covariates.csv"),
            sep = ",", row.names = FALSE, quote = FALSE)

save(covariate_table,
     file = paste0(Sys.Date(), "_MeMo_AC_VCD_Foxl2_pseudobulk_mito_covariates.RData"))

save(deseq.res.list.mito.only,
     file = paste0(Sys.Date(), "_MeMo_AC_VCD_Foxl2_PB_DESeq2_results_MITO_ADJUSTED.RData"))

status_table <- do.call(rbind, all_status); rownames(status_table) <- NULL
status_table <- status_table[, c("dataset","celltype","n_samples","n_genes_tested","design_terms","resid_df","flag")]
write.table(status_table, file = paste0(Sys.Date(), "_MeMo_AC_VCD_Foxl2_mito_only_sample_sizes_and_DF.csv"),
            sep = ",", row.names = FALSE, quote = FALSE)
cat("\n=== Sample sizes / residual df, mito design ===\n")
print(status_table)

################################################################################
# 6. DEG-level comparison against baseline
#    (AC, Foxl2: original 2026-02-08 baseline; VCD: SVA-duration-fixed baseline)
################################################################################

load("~/Benayoun_lab/Projects/Menopause_model_project/Data/Ovarian_scRNAseq/PB_DESeq2/Without_scTE/AC/2026-02-08/DESeq2_results/2026-02-08_10x_ovary_Benayoun_lab_AC_PB_DESeq2_results.RData")
baseline.AC <- deseq.res.list
load("~/Benayoun_lab/Projects/Menopause_model_project/Data/Ovarian_scRNAseq/PB_DESeq2/Without_scTE/Foxl2/2026-02-08/DESeq2_results/2026-02-08_10x_ovary_Benayoun_lab_Foxl2_NULL_PB_DESeq2_object.RData")
baseline.Foxl2 <- deseq.res.list
load("~/Benayoun_lab/Projects/Menopause_model_project/Data/Ovarian_scRNAseq/PB_DESeq2/Without_scTE/VCD/2026-08-17_SVA_duration_protection_fix/2026-08-17_MeMo_VCD_PB_DESeq2_results_SVA_DURATION_FIXED.RData")
baseline.VCD <- deseq.res.list.FIXED
rm(deseq.res.list, deseq.res.list.FIXED); gc()

baselines <- list(AC = baseline.AC, VCD = baseline.VCD, Foxl2 = baseline.Foxl2)

comparison <- list()
for (ds in c("AC", "VCD", "Foxl2")) {
  for (ct in SHARED_CELLTYPES) {
    base_res  <- baselines[[ds]][[ct]]
    mito_res  <- deseq.res.list.mito.only[[ds]][[ct]]

    if (is.null(base_res) || is.null(mito_res)) {
      comparison[[paste(ds, ct)]] <- data.frame(Dataset = ds, CellType = ct, flag = "MISSING",
        n_DEG_baseline = NA, n_DEG_mito_only = NA, pct_retained_mito_only = NA,
        pct_direction_concordant = NA)
      next
    }
    base <- as.data.frame(base_res); base <- base[!is.na(base$padj), ]
    mito <- as.data.frame(mito_res); mito <- mito[!is.na(mito$padj), ]

    deg_base <- rownames(base)[base$padj < ALPHA]
    deg_mito <- rownames(mito)[mito$padj < ALPHA]
    overlap  <- intersect(deg_base, deg_mito)
    pct_dir  <- if (length(overlap) > 0) {
      100 * sum(sign(base[overlap, "log2FoldChange"]) == sign(mito[overlap, "log2FoldChange"])) / length(overlap)
    } else NA

    comparison[[paste(ds, ct)]] <- data.frame(
      Dataset = ds, CellType = ct, flag = "OK",
      n_DEG_baseline = length(deg_base), n_DEG_mito_only = length(deg_mito),
      pct_retained_mito_only = ifelse(length(deg_base) > 0, 100 * length(overlap) / length(deg_base), NA),
      pct_direction_concordant = pct_dir
    )
  }
}

comparison_table <- do.call(rbind, comparison); rownames(comparison_table) <- NULL
write.table(comparison_table,
            file = paste0(Sys.Date(), "_MeMo_AC_VCD_Foxl2_DEG_comparison_BASELINE_vs_MITO.csv"),
            sep = ",", row.names = FALSE, quote = FALSE)
cat("\n=== DEG comparison: baseline vs mito ===\n")
print(comparison_table)

################################################################################
# 7. Strip plots (log2FC jitter per cell type, colored by significance/direction)
################################################################################

make_stripplot <- function(deseq.res.list, cell_type_order, outfile, ylab, ylim,
                            up_label, up_color, down_label, down_color, main_title) {

  available_cell_types <- intersect(cell_type_order, names(deseq.res.list))
  deseq.res.list <- deseq.res.list[match(available_cell_types, names(deseq.res.list))]

  xlab <- character(length(deseq.res.list))
  for (i in seq_along(deseq.res.list)) {
    sig_genes <- sum(deseq.res.list[[i]]$padj < 0.05, na.rm = TRUE)
    xlab[i] <- paste(available_cell_types[i], "\n(", sig_genes, " sig.)", sep = "")
  }

  pdf(outfile, width = 6, height = 5)
  par(mar = c(3.1, 4.1, 2, 1))
  par(oma = c(6, 2, 1, 1))

  plot(x = 1, y = 1, type = "n",
       xlim = c(0.5, length(deseq.res.list) + 0.5), ylim = ylim,
       axes = FALSE, xlab = "", ylab = ylab, main = main_title, cex.main = 0.9)

  abline(h = 0)
  abline(h = seq(ylim[1], ylim[2], by = 5)[seq(ylim[1], ylim[2], by = 5) != 0], lty = "dotted", col = "grey")

  for (i in seq_along(deseq.res.list)) {
    current_result <- deseq.res.list[[i]]
    if (is.null(current_result) || nrow(current_result) == 0) next

    sig_genes <- current_result$padj < 0.05
    colors <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), nrow(current_result))
    colors[sig_genes & current_result$log2FoldChange > 0] <- up_color
    colors[sig_genes & current_result$log2FoldChange < 0] <- down_color

    points(x = jitter(rep(i, nrow(current_result)), amount = 0.2),
           y = rev(current_result$log2FoldChange),
           pch = 16, col = rev(colors), cex = 0.5, bg = rev(colors))
  }

  axis(1, at = 1:length(deseq.res.list), tick = FALSE, las = 2, lwd = 0, labels = xlab, cex.axis = 0.7)
  axis(2, las = 1, at = seq(ylim[1], ylim[2], 5))
  legend("topright", legend = c(up_label, down_label), col = c(up_color, down_color),
         pch = 16, pt.cex = 1, cex = 0.6, bty = "n")
  box()
  dev.off()
}

set.seed(123123123)

make_stripplot(
  deseq.res.list.mito.only[["AC"]], SHARED_CELLTYPES,
  outfile = paste0(Sys.Date(), "_MeMo_AC_stripplot_OF_vs_YF_MITO_ADJUSTED.pdf"),
  ylab = "Log2 FC", ylim = c(-20, 20),
  up_label = "Up in OF", up_color = "deeppink4",
  down_label = "Up in YF", down_color = "deeppink1",
  main_title = "AC (mito%-adjusted): OF vs YF"
)

make_stripplot(
  deseq.res.list.mito.only[["VCD"]], SHARED_CELLTYPES,
  outfile = paste0(Sys.Date(), "_MeMo_VCD_stripplot_VCD_vs_CTL_MITO_ADJUSTED.pdf"),
  ylab = "Log2 fold change (VCD / CTL)", ylim = c(-15, 20),
  up_label = "Up in VCD", up_color = "yellow",
  down_label = "Up in CTL", down_color = "deeppink1",
  main_title = "VCD (SVA-fixed, mito%-adjusted): VCD vs CTL"
)

make_stripplot(
  deseq.res.list.mito.only[["Foxl2"]], SHARED_CELLTYPES,
  outfile = paste0(Sys.Date(), "_MeMo_Foxl2_stripplot_het_vs_wt_MITO_ADJUSTED.pdf"),
  ylab = "Log2 fold change (het / wt)", ylim = c(-15, 20),
  up_label = "Up in het", up_color = "springgreen",
  down_label = "Up in wt", down_color = "deeppink1",
  main_title = "Foxl2 (mito%-adjusted): het vs wt"
)

################################################################################
sink(file = paste0(Sys.Date(), "_Perc_mito_covariate_analysis_session_info.txt"))
sessionInfo()
sink()
