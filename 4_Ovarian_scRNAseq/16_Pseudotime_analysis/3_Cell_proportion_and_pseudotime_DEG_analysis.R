library(monocle3)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ARTool)
library(patchwork)
library(pheatmap)
library(viridisLite)

set.seed(1234)

###############################################################################
# Menopause-model project / revision
# Two downstream pseudotime analyses (granulosa, theca): 
#    1. Marker-based subtype composition 
#    2. Genes differentially expressed along pseudotime 
###############################################################################

cell_types <- list(
  granulosa = list(
    cds_path = "2026-08-13_granulosa_cds_rooted.rds",
    subtype_markers = list(
      Preantral = c("Kctd14"),
      Antral    = c("Slc18a2", "Rgs2", "Ccnd2", "Serpine2", "Inhbb", "Cyp19a1"),
      Atretic   = c("Ghr", "Pik3r1")
    ),
    subtype_levels = c("Preantral", "Antral", "Atretic", "Other")
  ),
  theca = list(
    cds_path = "2026-08-13_theca_cds_rooted.rds",
    subtype_markers = list(
      TPC = c("Wt1", "Ptch1"),
      TC1 = c("Cyp17a1"),
      TC2 = c("Star", "Lhcgr")
    ),
    subtype_levels = c("TPC", "TC1", "TC2", "Other")
  )
)

age_bin_levels <- c("<9m", "9-15m", ">15m")

condition_label <- tribble(
  ~Model,   ~Group,          ~Condition_label,
  "VCD",    "Control",       "CTL",
  "VCD",    "Perturbation",  "VCD",
  "Foxl2",  "Control",       "wt",
  "Foxl2",  "Perturbation",  "het"
)

###############################################################################
# Helpers
###############################################################################

# Model/Group come from existing metadata columns
add_group_model <- function(df) {
  df %>%
    mutate(
      Model = recode(Dataset, AC = "Aging", VCD = "VCD", Foxl2 = "Foxl2"),
      Group = case_when(
        Model == "Aging" ~ "Control",
        Model == "VCD"   ~ ifelse(Treatment == "VCD", "Perturbation", "Control"),
        Model == "Foxl2" ~ ifelse(Genotype == "het",   "Perturbation", "Control"),
        TRUE ~ NA_character_
      ),
      Group = factor(Group, levels = c("Control", "Perturbation"))
    )
}

# Classify each cell into a subtype by marker positivity 
classify_subtype <- function(norm_counts, cells, marker_list, levels) {
  is_positive <- lapply(marker_list, function(genes) {
    genes <- intersect(genes, rownames(norm_counts))
    if (length(genes) == 0) return(rep(FALSE, length(cells)))
    mat <- norm_counts[genes, cells, drop = FALSE]
    Matrix::colSums(mat > 0) > 0
  })
  subtype <- rep("Other", length(cells))
  for (nm in rev(names(marker_list))) {
    subtype[is_positive[[nm]]] <- nm
  }
  factor(subtype, levels = levels)
}

# Library-level ART ANOVA
run_art_by_model <- function(lib_df, response_col, extra_group_col = NULL,
                              age_col = "Age_numeric", age_term_name = "Age_factor") {
  results <- list()
  split_col <- if (is.null(extra_group_col)) "Model" else c("Model", extra_group_col)
  groups <- lib_df %>% distinct(across(all_of(split_col)))
  for (i in seq_len(nrow(groups))) {
    key <- groups[i, , drop = FALSE]
    key_label <- paste(sapply(key, as.character), collapse = " / ")
    sub_df <- dplyr::inner_join(lib_df, key, by = split_col)
    sub_df[[age_term_name]] <- droplevels(factor(sub_df[[age_col]]))
    sub_df$Group <- droplevels(factor(sub_df$Group, levels = c("Control", "Perturbation")))
    if (nlevels(sub_df[[age_term_name]]) < 2 || length(unique(sub_df[[response_col]])) < 2) {
      cat("  Skipping", key_label, "-- insufficient factor levels (n =", nrow(sub_df), ")\n")
      next
    }
    formula_str <- if (nlevels(sub_df$Group) < 2) {
      paste(response_col, "~", age_term_name)
    } else {
      paste(response_col, "~", age_term_name, "* Group")
    }
    fit <- tryCatch(
      art(as.formula(formula_str), data = sub_df),
      error = function(e) { cat("    skipped (", conditionMessage(e), ")\n"); NULL }
    )
    if (is.null(fit)) next
    a <- as.data.frame(anova(fit))
    for (col in names(key)) a[[col]] <- key[[col]]
    results[[key_label]] <- a
  }
  bind_rows(results)
}

lib_subtype_all <- list()

###############################################################################
# Per cell type: marker-based subtype classification + composition by age bin, faceted by Model x Group
###############################################################################

for (label in names(cell_types)) {

  cfg <- cell_types[[label]]

  cds <- readRDS(cfg$cds_path)

  cd <- as.data.frame(colData(cds)) %>%
    mutate(Cell = colnames(cds), Pseudotime = pseudotime(cds)) %>%
    add_group_model()

  cd <- cd[is.finite(cd$Pseudotime), ]
  cd$Age_factor <- factor(cd$Age_numeric)

  cd <- cd %>%
    mutate(
      Age_bin = case_when(
        Age_numeric < 9  ~ "<9m",
        Age_numeric <= 15 ~ "9-15m",
        TRUE ~ ">15m"
      ),
      Age_bin = factor(Age_bin, levels = age_bin_levels)
    )

  norm_counts <- normalized_counts(cds)
  cd$Subtype <- classify_subtype(norm_counts, cd$Cell, cfg$subtype_markers, cfg$subtype_levels)

  subtype_pal <- setNames(
    c("#440154", "#31688E", "#35B779", "#CCCCCC")[seq_along(cfg$subtype_levels)],
    cfg$subtype_levels
  )

  lib_subtype <- cd %>%
    count(Model, Library, Age_numeric, Age_bin, Group, Subtype, name = "n") %>%
    complete(Subtype, nesting(Model, Library, Age_numeric, Age_bin, Group), fill = list(n = 0)) %>%
    group_by(Model, Library) %>%
    mutate(total = sum(n), pct = 100 * n / total) %>%
    ungroup() %>%
    filter(!is.na(Model))

  write.csv(lib_subtype, paste0(Sys.Date(), "_", label, "_subtype_proportion_by_library.csv"),
            row.names = FALSE)
  lib_subtype_all[[label]] <- lib_subtype

  art_subtype <- run_art_by_model(lib_subtype, "pct", extra_group_col = "Subtype",
                                   age_col = "Age_bin", age_term_name = "Age_bin")
  write.csv(art_subtype, paste0(Sys.Date(), "_", label, "_subtype_ART_ANOVA.csv"), row.names = FALSE)

  # Age_bin x Group interaction 
  subtype_p_label <- art_subtype %>%
    group_by(Model) %>%
    group_modify(~ {
      has_group <- "Group" %in% .x$Term
      term_wanted <- if (has_group) "Age_bin:Group" else "Age_bin"
      term_short  <- if (has_group) "Age x Grp p" else "Age p"
      lines <- vapply(cfg$subtype_levels, function(st) {
        p <- .x[["Pr(>F)"]][.x$Term == term_wanted & .x$Subtype == st]
        if (length(p) == 0) return(NA_character_)
        sprintf("%s: %s = %.3g", st, term_short, p)
      }, character(1))
      tibble(label = paste(lines[!is.na(lines)], collapse = "\n"))
    }) %>%
    ungroup() %>%
    mutate(Group = factor("Control", levels = levels(cd$Group)))

  prop_subtype <- cd %>%
    count(Model, Age_bin, Group, Subtype, .drop = FALSE) %>%
    group_by(Model, Age_bin, Group) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    filter(!is.na(Model), !is.nan(pct))

  fig_subtype_bar <- ggplot(prop_subtype, aes(x = Age_bin, y = pct, fill = Subtype)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(pct >= 4, paste0(round(pct, 1), "%"), "")),
              position = position_stack(vjust = 0.5), size = 2.8, color = "white", fontface = "bold") +
    geom_text(data = subtype_p_label, aes(x = -Inf, y = 101, label = label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1, size = 2.2) +
    scale_fill_manual(values = subtype_pal, name = "Subtype") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 145), breaks = c(0, 25, 50, 75, 100)) +
    facet_grid(Group ~ Model, scales = "free_x", space = "free_x") +
    theme_bw() +
    labs(x = "Age bin", y = "Proportion of cells (%)",
         title = paste0(label, " -- subtype composition by age bin"),
         caption = paste0("Priority: ", paste(cfg$subtype_levels, collapse = " > "),
                           ". P-values: library-level ART ANOVA on subtype %, per Model x Subtype"))

  ggsave(paste0(Sys.Date(), "_", label, "_subtype_composition.pdf"), fig_subtype_bar,
         width = 10, height = 7.5)

  write.csv(cd %>% select(-any_of(c("UMAP1", "UMAP2"))),
            paste0(Sys.Date(), "_", label, "_pseudotime_age_subtype_percell.csv"),
            row.names = FALSE)

  rm(cds, cd, norm_counts, lib_subtype, prop_subtype)
  gc()
}

###############################################################################
# Controls only, pooled across models 
###############################################################################

for (label in names(cell_types)) {

  cfg <- cell_types[[label]]

  lib_subtype <- lib_subtype_all[[label]] %>% filter(Group == "Control")

  prop_subtype <- lib_subtype %>%
    group_by(Age_bin, Subtype, .drop = FALSE) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(Age_bin) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

  subtype_pal <- setNames(
    c("#440154", "#31688E", "#35B779", "#CCCCCC")[seq_along(cfg$subtype_levels)],
    cfg$subtype_levels
  )

  art_subtype <- list()
  for (st in cfg$subtype_levels) {
    sub_df <- lib_subtype %>% filter(Subtype == st) %>% mutate(Age_bin = droplevels(Age_bin))
    if (nlevels(sub_df$Age_bin) < 2) next
    fit <- tryCatch(art(pct ~ Age_bin, data = sub_df), error = function(e) NULL)
    if (is.null(fit)) next
    a <- as.data.frame(anova(fit))
    a$Subtype <- st
    art_subtype[[st]] <- a
  }
  art_subtype <- bind_rows(art_subtype)
  write.csv(art_subtype, paste0(Sys.Date(), "_", label, "_subtype_composition_controls_pooled_ART_ANOVA.csv"),
            row.names = FALSE)

  subtype_p_label <- art_subtype %>%
    mutate(label = sprintf("%s: Age p = %.3g", Subtype, `Pr(>F)`)) %>%
    summarise(label = paste(label, collapse = "\n"))

  fig <- ggplot(prop_subtype, aes(x = Age_bin, y = pct, fill = Subtype)) +
    geom_col(width = 0.55, color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(pct >= 2, paste0(round(pct, 1), "%"), "")),
              position = position_stack(vjust = 0.5), size = 3.2, color = "white", fontface = "bold") +
    annotate("text", x = -Inf, y = 101, label = subtype_p_label$label, hjust = -0.05, vjust = 1, size = 2.6) +
    scale_fill_manual(values = subtype_pal, name = "Subtype") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 115), breaks = c(0, 25, 50, 75, 100)) +
    theme_bw() +
    labs(x = "Age bin", y = "Proportion of cells (%)",
         title = paste0(label, " -- subtype composition, controls only (AC + VCD + Foxl2 pooled)"),
         caption = paste0("Priority: ", paste(cfg$subtype_levels, collapse = " > "),
                           ". Control arms only, pooled across models per age bin.",
                           " P-values: library-level ART ANOVA on subtype %% ~ Age_bin (Model pooled)"))

  ggsave(paste0(Sys.Date(), "_", label, "_subtype_composition_controls_pooled.pdf"), fig,
         width = 6.5, height = 6.5)
}

###############################################################################
# VCD and Foxl2 perturbation vs. their own control
###############################################################################

for (label in names(cell_types)) {

  cfg <- cell_types[[label]]

  lib_subtype <- lib_subtype_all[[label]] %>%
    filter(Model %in% c("VCD", "Foxl2")) %>%
    left_join(condition_label, by = c("Model", "Group")) %>%
    mutate(
      Age_factor      = factor(Age_numeric),
      Condition_label = factor(Condition_label, levels = c("CTL", "VCD", "wt", "het"))
    )

  subtype_pal <- setNames(
    c("#440154", "#31688E", "#35B779", "#CCCCCC")[seq_along(cfg$subtype_levels)],
    cfg$subtype_levels
  )

  # Cells summed across libraries within each Model x Age x Group
  prop_subtype <- lib_subtype %>%
    group_by(Model, Age_factor, Group, Condition_label, Subtype, .drop = FALSE) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(Model, Age_factor, Group) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    filter(!is.nan(pct))

  art_subtype <- list()
  for (mdl in c("VCD", "Foxl2")) {
    for (st in cfg$subtype_levels) {
      sub_df <- lib_subtype %>% filter(Model == mdl, Subtype == st) %>%
        mutate(Age_factor = droplevels(Age_factor), Group = droplevels(Group))
      if (nlevels(sub_df$Age_factor) < 2 || nlevels(sub_df$Group) < 2) next
      fit <- tryCatch(art(pct ~ Age_factor * Group, data = sub_df), error = function(e) NULL)
      if (is.null(fit)) next
      a <- as.data.frame(anova(fit))
      a$Model <- mdl
      a$Subtype <- st
      art_subtype[[paste(mdl, st)]] <- a
    }
  }
  art_subtype <- bind_rows(art_subtype)
  print(art_subtype)
  write.csv(art_subtype, paste0(Sys.Date(), "_", label, "_subtype_composition_perturbation_vs_control_ART_ANOVA.csv"),
            row.names = FALSE)

  subtype_p_label <- art_subtype %>%
    filter(Term == "Age_factor:Group") %>%
    group_by(Model) %>%
    group_modify(~ {
      lines <- vapply(cfg$subtype_levels, function(st) {
        p <- .x[["Pr(>F)"]][.x$Subtype == st]
        if (length(p) == 0) return(NA_character_)
        sprintf("%s: Age x Grp p = %.3g", st, p)
      }, character(1))
      tibble(label = paste(lines[!is.na(lines)], collapse = "\n"))
    }) %>%
    ungroup()

  # VCD and Foxl2 
  build_model_panel <- function(mdl) {
    df <- prop_subtype %>% filter(Model == mdl) %>% mutate(Age_factor = droplevels(Age_factor))
    lbl <- subtype_p_label %>% filter(Model == mdl) %>%
      mutate(Age_factor = factor(levels(df$Age_factor)[1], levels = levels(df$Age_factor)))
    ggplot(df, aes(x = Condition_label, y = pct, fill = Subtype)) +
      geom_col(width = 0.65, color = "white", linewidth = 0.3) +
      geom_text(aes(label = ifelse(pct >= 4, paste0(round(pct, 1), "%"), "")),
                position = position_stack(vjust = 0.5), size = 2.8, color = "white", fontface = "bold") +
      geom_text(data = lbl, aes(x = -Inf, y = 101, label = label),
                inherit.aes = FALSE, hjust = -0.05, vjust = 1, size = 2.2) +
      scale_fill_manual(values = subtype_pal, name = "Subtype") +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 130), breaks = c(0, 25, 50, 75, 100)) +
      facet_grid(Model ~ Age_factor) +
      theme_bw() +
      labs(x = "Condition", y = "Proportion of cells (%)")
  }

  panels <- lapply(c("VCD", "Foxl2"), build_model_panel)

  fig <- wrap_plots(panels, ncol = 1, guides = "collect") +
    plot_annotation(
      title = paste0(label, " -- subtype composition, perturbation vs. its own control, by age"),
      caption = paste0("Priority: ", paste(cfg$subtype_levels, collapse = " > "),
                        ". Every distinct age shown separately. P-values: library-level ART ANOVA,",
                        " Age x Group interaction, per Model x Subtype")
    )

  ggsave(paste0(Sys.Date(), "_", label, "_subtype_composition_perturbation_vs_control.pdf"), fig,
         width = 10, height = 8)
}

###############################################################################
sink(file = paste0(Sys.Date(), "_3_Cell_proportion_and_pseudotime_DEG_analysis_session_info.txt"))
sessionInfo()
sink()
