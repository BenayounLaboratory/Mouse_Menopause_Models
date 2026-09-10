library(readxl)
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(ARTool)
library(vegan)
library(effsize)
library(beeswarm)

rm(list = ls())

#######################################
# Menopause-model project / revision
# PCA per model (VCD, Foxl2), using follicle counts AND serum hormone
# levels (AMH, FSH, INHBA) jointly, fit on each model's OWN internal
# controls (CTL for VCD, wt for Foxl2), with treated animals (VCD, het)
# projected onto that control-defined space. Features are MFA
# block-weighted (5 follicle vars vs 3 hormone vars). 
#######################################

follicle_vars <- c("Primordial", "Primary", "Secondary", "Antral", "CL")
hormone_vars  <- c("AMH", "FSH", "INHBA")
all_vars <- c(follicle_vars, hormone_vars)
blocks   <- list(follicle = follicle_vars, hormone = hormone_vars)

VCD.age.levels   <- c("3m", "6m", "8m", "10m")
Foxl2.age.levels <- c("Y", "O", "SO")

# Palettes 
CTL.palette <- setNames(c("#ed2a91", "#d01779", "#b81f6c", "#a21f60"), VCD.age.levels)
VCD.palette <- setNames(c("#f6eb16", "#e0d120", "#dacb25", "#c6b92f"), VCD.age.levels)
wt.palette  <- setNames(c("#ed2a91", "#b01f68", "#600a38"), Foxl2.age.levels)
het.palette <- setNames(c("#72c16a", "#098c45", "#166635"), Foxl2.age.levels)

# Shape by age
VCD.shape.map   <- setNames(c(21, 22, 23, 24), VCD.age.levels)
Foxl2.shape.map <- setNames(c(21, 22, 24), Foxl2.age.levels)

# Boxplot colors
VCD.colors   <- c("deeppink", "yellow")
Foxl2.colors <- c("deeppink", "springgreen")

fmt <- function(p) formatC(p, format = "e", digits = 2)

#######################################
# 1. Import data
#######################################

VCD.follicle.raw   <- read.table("~/Dropbox/Benayoun_lab/Menopause_model_project/Hormone_histology_OHI/Input_data/2025-02-27_MeMo_VCD_follicle_data.txt", header = TRUE)
VCD.hormone.raw    <- read.table("~/Dropbox/Benayoun_lab/Menopause_model_project/Hormone_histology_OHI/Input_data/2025-02-27_MeMo_VCD_hormone_data.txt", header = TRUE)
Foxl2.follicle.raw <- read.table("~/Dropbox/Benayoun_lab/Menopause_model_project/Hormone_histology_OHI/Input_data/2026-02-11_Foxl2_follicle_medians_updated.txt", header = TRUE)
Foxl2.hormone.raw  <- read.table("~/Dropbox/Benayoun_lab/Menopause_model_project/Hormone_histology_OHI/Input_data/2026-02-11_MeMo_Foxl2_hormone_data_updated.txt", header = TRUE)

state_to_duration <- c(Post_I_1M = "30d", Post_I_3M = "90d", Post_I_5M = "150d")

make_VCD_timepoint <- function(state_label, duration_label) {
  h.tp <- VCD.hormone.raw[VCD.hormone.raw$State == state_label, ]
  f.tp <- VCD.follicle.raw[VCD.follicle.raw$Mouse_ID %in% h.tp$Mouse_ID, ]
  rownames(h.tp) <- h.tp$Mouse_ID
  rownames(f.tp) <- f.tp$Mouse_ID
  common <- intersect(rownames(h.tp), rownames(f.tp))
  dat <- cbind(f.tp[common, ], h.tp[common, c("Age_at_injection", "Treatment", hormone_vars)])
  dat$Treatment <- ifelse(dat$Treatment == "Safflower_oil", "CTL", dat$Treatment)
  dat %>% rename(GroupVar = Treatment) %>%
    mutate(AgeVar = tolower(Age_at_injection), Timepoint = duration_label)
}

VCD.combined <- bind_rows(
  make_VCD_timepoint("Post_I_1M", "30d"),
  make_VCD_timepoint("Post_I_3M", "90d"),
  make_VCD_timepoint("Post_I_5M", "150d")
) %>% mutate(GroupVar = factor(GroupVar, levels = c("CTL", "VCD")), AgeVar = factor(AgeVar, levels = VCD.age.levels))

Foxl2.hormone.N <- Foxl2.hormone.raw[Foxl2.hormone.raw$Issues == "N", ]
Foxl2.matched.samples <- intersect(Foxl2.follicle.raw$Mouse_ID, Foxl2.hormone.N$Mouse_ID)
Foxl2.h.cl <- Foxl2.hormone.N[Foxl2.hormone.N$Mouse_ID %in% Foxl2.matched.samples, c("Mouse_ID", "Genotype", "Age_group", hormone_vars)]
Foxl2.f.cl <- Foxl2.follicle.raw[Foxl2.follicle.raw$Mouse_ID %in% Foxl2.matched.samples, ]

Foxl2.combined <- Foxl2.f.cl %>% inner_join(Foxl2.h.cl, by = "Mouse_ID") %>%
  filter(Age_group %in% Foxl2.age.levels) %>%
  rename(GroupVar = Genotype, AgeVar = Age_group) %>%
  mutate(GroupVar = factor(GroupVar, levels = c("wt", "het")), AgeVar = factor(AgeVar, levels = Foxl2.age.levels), Timepoint = NA_character_)

VCD.ctrl   <- VCD.combined %>% filter(GroupVar == "CTL")
VCD.trt    <- VCD.combined %>% filter(GroupVar == "VCD")
Foxl2.ctrl <- Foxl2.combined %>% filter(GroupVar == "wt")
Foxl2.trt  <- Foxl2.combined %>% filter(GroupVar == "het")

#######################################
# 2. Helpers: MFA block-weighted PCA fit on internal controls (projecting
#    treated animals onto it), PC1 span metric, plotting, factorial tests
#######################################

# After z-scoring each variable (using CONTROL mean/sd), further divide each
# block's variables by sqrt(lambda1_block). 

fit_and_project_internal_MFA <- function(block_list, ctrl_data, trt_data) {
  feature_cols <- unlist(block_list, use.names = FALSE)
  ctrl_cc <- ctrl_data[complete.cases(ctrl_data[, feature_cols]), ]
  trt_cc  <- trt_data[complete.cases(trt_data[, feature_cols]), ]

  ctrl_mean <- sapply(ctrl_cc[, feature_cols], mean)
  ctrl_sd   <- sapply(ctrl_cc[, feature_cols], sd)

  block_weight <- setNames(numeric(length(feature_cols)), feature_cols)
  for (blk in names(block_list)) {
    vars_b <- block_list[[blk]]
    z_b <- scale(ctrl_cc[, vars_b], center = ctrl_mean[vars_b], scale = ctrl_sd[vars_b])
    lambda1 <- prcomp(z_b, center = FALSE, scale. = FALSE)$sdev[1]^2
    block_weight[vars_b] <- 1 / sqrt(lambda1)
  }

  prep <- function(d) {
    z <- scale(d[, feature_cols], center = ctrl_mean[feature_cols], scale = ctrl_sd[feature_cols])
    sweep(z, 2, block_weight[feature_cols], `*`)
  }

  ctrl_prepped <- prep(ctrl_cc)
  pca <- prcomp(ctrl_prepped, center = FALSE, scale. = FALSE)
  varexp <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)
  ctrl_scores <- cbind(ctrl_cc[, setdiff(colnames(ctrl_cc), feature_cols)], as.data.frame(pca$x))

  trt_prepped <- prep(trt_cc)
  proj <- predict(pca, newdata = trt_prepped)
  trt_scores <- cbind(trt_cc[, setdiff(colnames(trt_cc), feature_cols)], as.data.frame(proj))

  list(pca = pca, varexp = varexp, ctrl_scores = ctrl_scores, trt_scores = trt_scores, block_weight = block_weight)
}

internal_span <- function(ctrl_scores, age_col, youngest, oldest) {
  y <- median(ctrl_scores$PC1[ctrl_scores[[age_col]] == youngest], na.rm = TRUE)
  o <- median(ctrl_scores$PC1[ctrl_scores[[age_col]] == oldest], na.rm = TRUE)
  list(young = y, old = o, span = o - y)
}

pct_span_internal <- function(ctrl_scores, trt_scores, span, age_col, age_levels_ordered, timepoint_col = NULL) {
  has_tp <- !is.null(timepoint_col) && timepoint_col %in% colnames(ctrl_scores)
  tps <- if (has_tp) unique(na.omit(ctrl_scores[[timepoint_col]])) else NA_character_
  rows <- list()
  for (tp in tps) {
    ctrl_tp <- if (has_tp) ctrl_scores[ctrl_scores[[timepoint_col]] == tp, ] else ctrl_scores
    trt_tp  <- if (has_tp) trt_scores[trt_scores[[timepoint_col]] == tp, ]  else trt_scores
    young_ref <- median(ctrl_tp$PC1[ctrl_tp[[age_col]] == age_levels_ordered[1]], na.rm = TRUE)
    for (age in age_levels_ordered) {
      ctrl_age <- ctrl_tp$PC1[ctrl_tp[[age_col]] == age]
      trt_age  <- trt_tp$PC1[trt_tp[[age_col]] == age]
      if (length(ctrl_age) == 0) next
      cd <- if (length(trt_age) >= 2 && length(ctrl_age) >= 2) cliff.delta(trt_age, ctrl_age) else NULL
      rows[[length(rows) + 1]] <- data.frame(
        Timepoint = tp, Age = age,
        pct_CTL = 100 * (median(ctrl_age, na.rm = TRUE) - young_ref) / span,
        pct_Treated = if (length(trt_age) > 0) 100 * (median(trt_age, na.rm = TRUE) - young_ref) / span else NA_real_,
        gap_pct_of_span = if (length(trt_age) > 0) 100 * (median(trt_age, na.rm = TRUE) - median(ctrl_age, na.rm = TRUE)) / span else NA_real_,
        cliffs_delta = if (!is.null(cd)) as.numeric(cd$estimate) else NA_real_,
        magnitude = if (!is.null(cd)) as.character(cd$magnitude) else NA_character_,
        n_CTL = length(ctrl_age), n_Treated = length(trt_age)
      )
    }
  }
  do.call(rbind, rows)
}

plot_internal_projection <- function(ctrl_scores, trt_scores, age_col, age_levels, ctrl_label, trt_label,
                                      ctrl_palette, trt_palette, shape_map, title, subtitle, out_file, varexp, facet_col = NULL) {
  ctrl_scores[[age_col]] <- factor(ctrl_scores[[age_col]], levels = age_levels)
  trt_scores[[age_col]]  <- factor(trt_scores[[age_col]], levels = age_levels)

  # Dummy invisible layer for legend
  group_legend_df <- data.frame(PC1 = NA_real_, PC2 = NA_real_,
                                 Group = factor(c(ctrl_label, trt_label), levels = c(ctrl_label, trt_label)))
  group_colors <- setNames(c(ctrl_palette[[length(ctrl_palette)]], trt_palette[[length(trt_palette)]]), c(ctrl_label, trt_label))

  p <- ggplot() +
    geom_point(data = ctrl_scores, aes(x = PC1, y = PC2, fill = .data[[age_col]], shape = .data[[age_col]]), size = 2.6, alpha = 0.85, color = "grey20") +
    scale_fill_manual(values = ctrl_palette, guide = "none") +
    scale_shape_manual(values = shape_map, name = "Age") +
    ggnewscale::new_scale_fill() +
    geom_point(data = trt_scores, aes(x = PC1, y = PC2, fill = .data[[age_col]], shape = .data[[age_col]]), size = 2.6, alpha = 0.85, color = "grey20") +
    scale_fill_manual(values = trt_palette, guide = "none") +
    geom_point(data = group_legend_df, aes(x = PC1, y = PC2, color = Group), size = 3, na.rm = TRUE) +
    scale_color_manual(values = group_colors, name = "Group") +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 3))) +
    labs(x = paste0("PC1 (", varexp[1], "%, control-defined aging axis)"), y = paste0("PC2 (", varexp[2], "%)"),
         title = title, subtitle = subtitle) +
    theme_bw(base_size = 11)
  if (!is.null(facet_col)) { p <- p + facet_wrap(vars(.data[[facet_col]])); ggsave(out_file, p, width = 12, height = 5) }
  else ggsave(out_file, p, width = 6.5, height = 5.5)
}

run_factorial_tests <- function(ctrl_scores, trt_scores, ctrl_level, trt_level, timepoint_col = NULL) {
  combined <- bind_rows(ctrl_scores, trt_scores) %>% droplevels()
  pc_cols <- grep("^PC", colnames(combined), value = TRUE)
  has_tp <- !is.null(timepoint_col) && timepoint_col %in% colnames(combined) && !all(is.na(combined[[timepoint_col]]))
  tps <- if (has_tp) unique(na.omit(combined[[timepoint_col]])) else "ALL"

  out <- lapply(tps, function(tp) {
    d <- if (has_tp) combined[combined[[timepoint_col]] == tp, ] else combined
    d <- droplevels(d)
    art_out <- tryCatch({
      fit <- art(PC1 ~ GroupVar * AgeVar, data = d)
      aov_tab <- anova(fit)
      data.frame(Timepoint = tp, Test = "ART_ANOVA_PC1", Term = c("Group", "Age", "Group:Age"),
                 p_value = c(aov_tab["GroupVar", "Pr(>F)"], aov_tab["AgeVar", "Pr(>F)"], aov_tab["GroupVar:AgeVar", "Pr(>F)"]), R2 = NA_real_)
    }, error = function(e) data.frame(Timepoint = tp, Test = "ART_ANOVA_PC1", Term = c("Group","Age","Group:Age"), p_value = NA_real_, R2 = NA_real_))

    perm_out <- tryCatch({
      dmat <- dist(as.matrix(d[, pc_cols]), method = "euclidean")
      ad <- adonis2(dmat ~ GroupVar * AgeVar, data = d, permutations = 999, by = "terms")
      data.frame(Timepoint = tp, Test = "PERMANOVA_allPCs", Term = c("Group", "Age", "Group:Age"),
                 p_value = ad[c("GroupVar", "AgeVar", "GroupVar:AgeVar"), "Pr(>F)"],
                 R2 = ad[c("GroupVar", "AgeVar", "GroupVar:AgeVar"), "R2"])
    }, error = function(e) data.frame(Timepoint = tp, Test = "PERMANOVA_allPCs", Term = c("Group","Age","Group:Age"), p_value = NA_real_, R2 = NA_real_))

    bind_rows(art_out, perm_out)
  })
  bind_rows(out)
}

#######################################
# 3. VCD - MFA-weighted PCA fit on CTL, project VCD, per Duration
#######################################

res_VCD_MFA <- fit_and_project_internal_MFA(blocks, VCD.ctrl, VCD.trt)

plot_internal_projection(res_VCD_MFA$ctrl_scores, res_VCD_MFA$trt_scores, "AgeVar", VCD.age.levels, "CTL", "VCD",
                          CTL.palette, VCD.palette, VCD.shape.map,
                          "VCD combined follicle+hormone (MFA block-weighted) projected onto internal-CTL PCA",
                          "Shape = age (shared across groups), color = group (pink=CTL, yellow=VCD). Block-balanced 5-follicle/3-hormone. Matched samples only.",
                          paste0(Sys.Date(), "_MeMo_VCD_combined_MFA_projected_onto_internalCTL_PCA.pdf"),
                          varexp = res_VCD_MFA$varexp, facet_col = "Timepoint")

span_VCD_MFA <- internal_span(res_VCD_MFA$ctrl_scores, "AgeVar", "3m", "10m")
VCD.span.table.MFA <- pct_span_internal(res_VCD_MFA$ctrl_scores, res_VCD_MFA$trt_scores, span_VCD_MFA$span, "AgeVar", VCD.age.levels, "Timepoint")
VCD.factorial.MFA  <- run_factorial_tests(res_VCD_MFA$ctrl_scores, res_VCD_MFA$trt_scores, "CTL", "VCD", "Timepoint")

write.table(VCD.span.table.MFA, paste0(Sys.Date(), "_MeMo_VCD_combined_internalCTL_MFA_pct_span.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(VCD.factorial.MFA,  paste0(Sys.Date(), "_MeMo_VCD_combined_internalCTL_MFA_factorial_tests.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

#######################################
# 4. Foxl2 - MFA-weighted PCA fit on wt, project het
#######################################

res_Foxl2_MFA <- fit_and_project_internal_MFA(blocks, Foxl2.ctrl, Foxl2.trt)

plot_internal_projection(res_Foxl2_MFA$ctrl_scores, res_Foxl2_MFA$trt_scores, "AgeVar", Foxl2.age.levels, "wt", "het",
                          wt.palette, het.palette, Foxl2.shape.map,
                          "Foxl2 combined follicle+hormone (MFA block-weighted) projected onto internal-wt PCA",
                          "Shape = age group (shared across groups), color = genotype (pink=wt, green=het). Block-balanced 5-follicle/3-hormone. Matched samples only.",
                          paste0(Sys.Date(), "_MeMo_Foxl2_combined_MFA_projected_onto_internalWT_PCA.pdf"),
                          varexp = res_Foxl2_MFA$varexp)

span_Foxl2_MFA <- internal_span(res_Foxl2_MFA$ctrl_scores, "AgeVar", "Y", "SO")
Foxl2.span.table.MFA <- pct_span_internal(res_Foxl2_MFA$ctrl_scores, res_Foxl2_MFA$trt_scores, span_Foxl2_MFA$span, "AgeVar", Foxl2.age.levels)
Foxl2.factorial.MFA  <- run_factorial_tests(res_Foxl2_MFA$ctrl_scores, res_Foxl2_MFA$trt_scores, "wt", "het")

write.table(Foxl2.span.table.MFA, paste0(Sys.Date(), "_MeMo_Foxl2_combined_internalWT_MFA_pct_span.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Foxl2.factorial.MFA,  paste0(Sys.Date(), "_MeMo_Foxl2_combined_internalWT_MFA_factorial_tests.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

#######################################
# 5. PC1 boxplots by group + ART ANOVA
#######################################

plot_PC1_panel <- function(dat, group_col, age_col, group_levels, age_levels, group_label,
                            colors, pch, ylab, main) {
  dat <- dat %>%
    mutate(
      !!group_col := factor(.data[[group_col]], levels = group_levels),
      !!age_col   := factor(.data[[age_col]],   levels = age_levels),
      Group = factor(paste(.data[[group_col]], .data[[age_col]], sep = "_"),
                      levels = unlist(lapply(age_levels, function(a) paste(group_levels, a, sep = "_"))))
    )

  fit <- art(as.formula(paste0("PC1 ~ ", group_col, " * ", age_col)), data = dat)
  aov_tab <- anova(fit)
  label <- paste0(
    "ART ANOVA: ", group_label, " = ", fmt(aov_tab[group_col, "Pr(>F)"]), " | ",
    "Age = ", fmt(aov_tab[age_col, "Pr(>F)"]), " | ",
    group_label, "×Age = ", fmt(aov_tab[paste0(group_col, ":", age_col), "Pr(>F)"])
  )

  pvals <- sapply(age_levels, function(a) {
    sub <- dat[dat[[age_col]] == a, ]
    wilcox.test(as.formula(paste0("PC1 ~ ", group_col)), data = sub)$p.value
  })
  p.adj <- p.adjust(pvals, method = "BH")

  rng <- range(dat$PC1, na.rm = TRUE)
  pad <- diff(rng) * 0.18
  ylim <- c(rng[1] - pad * 0.3, rng[2] + pad)
  y.hi <- rng[2] + pad * 0.95
  y.lo <- rng[2] + pad * 0.55

  op <- par(mar = c(8, 4.5, 3, 1) + 0.1)
  boxplot(PC1 ~ Group, dat,
          outline = FALSE, ylim = ylim,
          col = rep(colors, length(age_levels)), las = 2,
          xlab = "", ylab = ylab, main = main)
  beeswarm(PC1 ~ Group, dat, pch = pch, col = "black", add = TRUE, cex = 1)

  x.mid <- seq(1.5, by = 2, length.out = length(age_levels))
  y.txt <- rep(c(y.hi, y.lo), length.out = length(age_levels))
  text(x.mid, y.txt, paste0("p=", fmt(p.adj)))

  mtext(label, side = 1, line = 6.5, cex = 0.8)
  par(op)

  invisible(list(art = aov_tab, wilcox_p = pvals, wilcox_p_adj = p.adj))
}

# VCD
VCD.dat <- bind_rows(res_VCD_MFA$ctrl_scores, res_VCD_MFA$trt_scores)

pdf(paste0(Sys.Date(), "_MeMo_VCD_combined_MFA_PC1_boxplot_internalCTL.pdf"), width = 8, height = 7)
VCD.stats <- list()
for (tp in c("30d", "90d", "150d")) {
  dat_tp <- VCD.dat %>% filter(Timepoint == tp)
  VCD.stats[[tp]] <- plot_PC1_panel(
    dat_tp, "GroupVar", "AgeVar", c("CTL", "VCD"), VCD.age.levels, "Group",
    VCD.colors, pch = 2,
    ylab = paste0("PC1 (", res_VCD_MFA$varexp[1], "%, MFA block-weighted, control-defined aging axis)"),
    main = paste0("VCD combined follicle+hormone (MFA), ", tp, " post-injection")
  )
}
dev.off()

# Foxl2
Foxl2.dat <- bind_rows(res_Foxl2_MFA$ctrl_scores, res_Foxl2_MFA$trt_scores)

pdf(paste0(Sys.Date(), "_MeMo_Foxl2_combined_MFA_PC1_boxplot_internalWT.pdf"), width = 8, height = 7)
Foxl2.stats <- plot_PC1_panel(
  Foxl2.dat, "GroupVar", "AgeVar", c("wt", "het"), Foxl2.age.levels, "Genotype",
  Foxl2.colors, pch = 5,
  ylab = paste0("PC1 (", res_Foxl2_MFA$varexp[1], "%, MFA block-weighted, control-defined aging axis)"),
  main = "Foxl2 combined follicle+hormone (MFA)"
)
dev.off()

#######################################
sink(file = paste0(Sys.Date(), "_PCA_per_model_internal_control_combined_session_info.txt"))
sessionInfo()
sink()
