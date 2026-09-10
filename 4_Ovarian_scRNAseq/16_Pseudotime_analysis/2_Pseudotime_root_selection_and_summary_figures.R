library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)

###############################################################################
# Menopause-model project / revision
# Interactive pseudotime root selection (granulosa, theca), summary figures,
# and marker-gene identity validation
###############################################################################

###############################################################################
# 1. Granulosa root selection
###############################################################################

if (!exists("cds_granulosa")) {
  granulosa_file <- list.files(pattern = "_granulosa_cds_prerooted\\.rds$", full.names = TRUE)
  granulosa_file <- granulosa_file[order(granulosa_file, decreasing = TRUE)]
  cds_granulosa <- readRDS(granulosa_file[1])
}

cds_granulosa <- order_cells(cds_granulosa)

# Sanity check: pseudotime should increase with Age_numeric 
cor.test(pseudotime(cds_granulosa), colData(cds_granulosa)$Age_numeric, method = "spearman")

plot_cells(cds_granulosa,
           color_cells_by          = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves            = FALSE,
           label_branch_points     = FALSE) +
  ggtitle("Granulosa - pseudotime (post-rooting)")

saveRDS(cds_granulosa, paste0(Sys.Date(), "_granulosa_cds_rooted.rds"))

###############################################################################
# 2. Theca root selection
###############################################################################

if (!exists("cds_theca")) {
  theca_file <- list.files(pattern = "_theca_cds_prerooted\\.rds$", full.names = TRUE)
  theca_file <- theca_file[order(theca_file, decreasing = TRUE)]
  cds_theca <- readRDS(theca_file[1])
}

cds_theca <- order_cells(cds_theca)

cor.test(pseudotime(cds_theca), colData(cds_theca)$Age_numeric, method = "spearman")

plot_cells(cds_theca,
           color_cells_by          = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves            = FALSE,
           label_branch_points     = FALSE) +
  ggtitle("Theca - pseudotime (post-rooting)")

saveRDS(cds_theca, paste0(Sys.Date(), "_theca_cds_rooted.rds"))

###############################################################################
# 3. Pseudotime summary figures
###############################################################################

extract_pt_df <- function(cds, cell_type_label) {
  data.frame(
    Cell            = colnames(cds),
    Pseudotime      = pseudotime(cds),
    Library         = colData(cds)$Library,
    Age_numeric     = colData(cds)$Age_numeric,
    Seurat_cluster  = colData(cds)$seurat_clusters,
    Monocle_cluster = as.character(clusters(cds)),
    CellType        = cell_type_label,
    stringsAsFactors = FALSE
  )
}

pt_all <- bind_rows(
  extract_pt_df(cds_granulosa, "Granulosa"),
  extract_pt_df(cds_theca,     "Theca")
)

n_before <- nrow(pt_all)
pt_all <- pt_all %>% filter(is.finite(Pseudotime))

pt_all <- pt_all %>%
  mutate(
    Group = case_when(
      grepl("VCD", Library) ~ "Perturbation",
      grepl("het", Library) ~ "Perturbation",
      TRUE                  ~ "Control"
    ),
    Model = case_when(
      grepl("^YF|^OF",       Library) ~ "Aging",
      grepl("^CTL_|^VCD_",   Library) ~ "VCD",
      grepl("^Foxl2_",       Library) ~ "Foxl2",
      TRUE                             ~ NA_character_
    )
  )

# Plot: controls, ordered by age

pt_control <- pt_all %>% filter(Group == "Control")

lib_order <- pt_control %>%
  distinct(Library, Age_numeric) %>%
  arrange(Age_numeric, Library) %>%
  pull(Library)
pt_control$Library <- factor(pt_control$Library, levels = lib_order)

fig1 <- ggplot(pt_control, aes(x = Library, y = Pseudotime, fill = Model)) +
  geom_violin(trim = TRUE, scale = "width") +
  facet_wrap(~ CellType, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Library (ordered by age)", y = "Pseudotime", fill = "Model",
       title = "Control samples: pseudotime by library, ordered by age")

ggsave(paste0(Sys.Date(), "_control_pseudotime_by_age.pdf"), fig1, width = 12, height = 8)

# Plot: VCD / Foxl2-het vs. age-matched controls

pt_pert <- pt_all %>% filter(Model %in% c("VCD", "Foxl2"))
pt_pert$Age_factor <- factor(pt_pert$Age_numeric)
pt_pert$Group <- factor(pt_pert$Group, levels = c("Control", "Perturbation"))

fig2 <- ggplot(pt_pert, aes(x = Age_factor, y = Pseudotime, fill = Group)) +
  geom_violin(position = position_dodge(width = 0.8), trim = TRUE, scale = "width") +
  facet_grid(CellType ~ Model, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(x = "Age (months)", y = "Pseudotime", fill = "Group",
       title = "Pseudotime: VCD / Foxl2-het vs. age-matched controls")

ggsave(paste0(Sys.Date(), "_perturbation_vs_control_pseudotime.pdf"), fig2, width = 10, height = 8)

###############################################################################
# 4. Marker gene expression
###############################################################################

# Granulosa cells
granulosa_markers <- c("Cyp19a1", "Foxl2", "Amh", "Fshr", "Lhcgr", "Inha", "Ptgs2", "Nr2f2")

if (!exists("granulosa_obj")) {
  f <- list.files(pattern = "_granulosa_reintegrated\\.rds$", full.names = TRUE)
  f <- f[order(f, decreasing = TRUE)]
  granulosa_obj <- readRDS(f[1])
}

DefaultAssay(granulosa_obj) <- "RNA"
granulosa_obj <- NormalizeData(granulosa_obj, verbose = FALSE)

present_markers <- intersect(granulosa_markers, rownames(granulosa_obj))
missing_markers <- setdiff(granulosa_markers, present_markers)

FeaturePlot(granulosa_obj,
            features = present_markers,
            ncol     = 4,
            pt.size  = 0.3,
            order    = TRUE,
            cols     = c("lightgrey", "darkred")) &
  theme(legend.position = "right",
        axis.title      = element_blank(),
        axis.text       = element_blank(),
        axis.ticks      = element_blank())

DimPlot(granulosa_obj, group.by = "seurat_clusters", label = TRUE)
DimPlot(granulosa_obj, group.by = "Age", label = TRUE)

# Theca cells
theca_markers <- c("Cyp17a1", "Star", "Insl3", "Nr5a1")

if (!exists("theca_obj")) {
  f <- list.files(path = c(".", ".."), pattern = "_theca_reintegrated\\.rds$", full.names = TRUE)
  f <- f[order(f, decreasing = TRUE)]
  theca_obj <- readRDS(f[1])
}

DefaultAssay(theca_obj) <- "RNA"
theca_obj <- NormalizeData(theca_obj, verbose = FALSE)

present_markers_theca <- intersect(theca_markers, rownames(theca_obj))
missing_markers_theca <- setdiff(theca_markers, present_markers_theca)

FeaturePlot(theca_obj,
            features = present_markers_theca,
            ncol     = 4,
            pt.size  = 0.3,
            order    = TRUE,
            cols     = c("lightgrey", "darkred")) &
  theme(legend.position = "right",
        axis.title      = element_blank(),
        axis.text       = element_blank(),
        axis.ticks      = element_blank())

DimPlot(theca_obj, group.by = "seurat_clusters", label = TRUE)
DimPlot(theca_obj, group.by = "Age", label = TRUE)

###############################################################################
sink(file = paste0(Sys.Date(), "_2_Pseudotime_root_selection_and_figures_session_info.txt"))
sessionInfo()
sink()
