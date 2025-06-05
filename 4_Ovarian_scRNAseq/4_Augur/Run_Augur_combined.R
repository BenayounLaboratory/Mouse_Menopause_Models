setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/4_Augur")
options(stringsAsFactors = F)

library(Seurat)
library(SeuratWrappers)
library('Augur')
library(viridis)
library(ggplot2)
library('ggrastr')
library(dplyr)
library(tidyr)
library(gridExtra)

# rm(list = ls())

################################################################################
# Perform Augur - cell type prioritization
# Aging, VCD and Foxl2 haploinsufficiency models
################################################################################

################################################################################
# 1. Load & subset data
################################################################################

#######################
# Aging (AC)
#######################

load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/4_Augur/Input/2024-10-15_10x_ovary_Benayoun_lab_AC_Seurat_object_with_final_annotation.RData")

#######################
# VCD
#######################

load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/4_Augur/Input/2024-10-24_10x_ovary_Benayoun_lab_VCD_Seurat_object_with_final_annotation.RData")

# Subset dataset

ovary.VCD.30d <- subset(ovary.VCD, subset = Duration == "30d")
ovary.VCD.90d <- subset(ovary.VCD, subset = Duration == "90d")

ovary.VCD.30d.3m <- subset(ovary.VCD.30d, subset = Age == "3m")
ovary.VCD.30d.10m <- subset(ovary.VCD.30d, subset = Age == "10m")
ovary.VCD.90d.3m <- subset(ovary.VCD.90d, subset = Age == "3m")
ovary.VCD.90d.10m <- subset(ovary.VCD.90d, subset = Age == "10m")

# Clear up memory
rm(ovary.VCD)

#######################
# Foxl2 haploinsufficiency
#######################

load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/4_Augur/Input/2024-10-24_10x_ovary_Benayoun_lab_Foxl2_Seurat_object_with_final_annotation.RData")

# Add Group metadata

Group <- rep("NA", length(colnames(ovary.Foxl2@assays$RNA)))
Group[grep("wt_young", colnames(ovary.Foxl2@assays$RNA))]  <- "wt_young"
Group[grep("wt_midage", colnames(ovary.Foxl2@assays$RNA))]    <- "wt_midage"
Group[grep("het_young", colnames(ovary.Foxl2@assays$RNA))] <- "het_young"
Group[grep("het_midage", colnames(ovary.Foxl2@assays$RNA))]   <- "het_midage"

Group <- data.frame(Group)
rownames(Group) <- colnames(ovary.Foxl2@assays$RNA)

# update Seurat with metadata

ovary.Foxl2 <- AddMetaData(object = ovary.Foxl2, metadata = as.vector(Group), col.name = "Group")

# Subset seurat based on age and genotype at injection

ovary.Foxl2.Young    <- subset(ovary.Foxl2, subset = Age == "Young")
ovary.Foxl2.Midage   <- subset(ovary.Foxl2, subset = Age == "Midage")

# Clear up memory
rm(ovary.Foxl2)

################################################################################
# 2. Run Augur
################################################################################

#######################
# Aging
#######################

augur.ovary.AC <- calculate_auc(as.matrix(ovary.AC@assays$SCT@data),
                                ovary.AC@meta.data, 
                                cell_type_col = "celltype.level2", 
                                label_col = "Age",
                                n_threads = 3)

save(augur.ovary.AC, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_all_cells_Augur_object.RData"))

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_all_cells_Augur_UMAP.pdf"), height = 5, width = 6)
plot_umap(augur.ovary.AC, ovary.AC, cell_type_col = "celltype.level2")
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_AC_all_cells_Augur_Lollipop.pdf"))
plot_lollipop(augur.ovary.AC)
dev.off()

# Clear up memory
rm(ovary.AC)

#######################
# VCD
#######################

# 30d, 3m dataset
augur.ovary.VCD.30d.3m <- calculate_auc(as.matrix(ovary.VCD.30d.3m@assays$SCT@data),
                                        ovary.VCD.30d.3m@meta.data, 
                                        cell_type_col = "celltype.level2", 
                                        label_col = "Treatment",
                                        n_threads = 3)
# augur.ovary.VCD.30d.3m

save(augur.ovary.VCD.30d.3m, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_30d_3m_all_cells_Augur_object.RData"))

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_30d_3m_all_cells_Augur_UMAP.pdf"), height = 5, width = 6)
plot_umap(augur.ovary.VCD.30d.3m, ovary.VCD.30d.3m, cell_type_col = "celltype.level2")
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_30d_3m_all_cells_Augur_Lollipop.pdf"))
plot_lollipop(augur.ovary.VCD.30d.3m)
dev.off()

rm(ovary.VCD.30d.3m)

# augur.ovary.VCD.90d.3m

augur.ovary.VCD.90d.3m <- calculate_auc(as.matrix(ovary.VCD.90d.3m@assays$SCT@data),
                                        ovary.VCD.90d.3m@meta.data, 
                                        cell_type_col = "celltype.level2", 
                                        label_col = "Treatment",
                                        n_threads = 3)

save(augur.ovary.VCD.90d.3m, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_90d_3m_all_cells_Augur_object.RData"))

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_90d_3m_all_cells_Augur_UMAP.pdf"), height = 5, width = 6)
plot_umap(augur.ovary.VCD.90d.3m, ovary.VCD.90d.3m, cell_type_col = "celltype.level2")
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_90d_3m_all_cells_Augur_Lollipop.pdf"))
plot_lollipop(augur.ovary.VCD.90d.3m)
dev.off()

rm(ovary.VCD.90d.3m)

# augur.ovary.VCD.30d.10m

augur.ovary.VCD.30d.10m <- calculate_auc(as.matrix(ovary.VCD.30d.10m@assays$SCT@data),
                                         ovary.VCD.30d.10m@meta.data, 
                                         cell_type_col = "celltype.level2", 
                                         label_col = "Treatment",
                                         n_threads = 3)

save(augur.ovary.VCD.30d.10m, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_30d_10m_all_cells_Augur_object.RData"))

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_30d_10m_all_cells_Augur_UMAP.pdf"), height = 5, width = 6)
plot_umap(augur.ovary.VCD.30d.10m, ovary.VCD.30d.10m, cell_type_col = "celltype.level2")
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_30d_10m_all_cells_Augur_Lollipop.pdf"))
plot_lollipop(augur.ovary.VCD.30d.10m)
dev.off()

rm(ovary.VCD.30d.10m)

# augur.ovary.VCD.90d.10m

augur.ovary.VCD.90d.10m <- calculate_auc(as.matrix(ovary.VCD.90d.10m@assays$SCT@data),
                                         ovary.VCD.90d.10m@meta.data, 
                                         cell_type_col = "celltype.level2", 
                                         label_col = "Treatment",
                                         n_threads = 3)

save(augur.ovary.VCD.90d.10m, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_90d_10m_all_cells_Augur_object.RData"))

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_90d_10m_all_cells_Augur_UMAP.pdf"), height = 5, width = 6)
plot_umap(augur.ovary.VCD.90d.10m, ovary.VCD.90d.10m, cell_type_col = "celltype.level2")
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_VCD_90d_10m_all_cells_Augur_Lollipop.pdf"))
plot_lollipop(augur.ovary.VCD.90d.10m)
dev.off()

rm(ovary.VCD.90d.10m)

#######################
# Foxl2 haploinsufficiency
#######################

# Young

augur.ovary.Foxl2.Young <-  calculate_auc(as.matrix(ovary.Foxl2.Young@assays$SCT@data),
                                          ovary.Foxl2.Young@meta.data, 
                                          cell_type_col = "celltype.level2", 
                                          label_col = "Group",
                                          n_threads = 3)

save(augur.ovary.Foxl2.Young, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_Young_Augur_object.RData"))

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_Young_Augur_UMAP.pdf"), height = 5, width = 6)
plot_umap(augur.ovary.Foxl2.Young, ovary.Foxl2.Young, cell_type_col = "celltype.level2")
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_Young_Augur_Lollipop.pdf"))
plot_lollipop(augur.ovary.Foxl2.Young)
dev.off()

# Mid-age

augur.ovary.Foxl2.Midage <-  calculate_auc(as.matrix(ovary.Foxl2.Midage@assays$SCT@data),
                                           ovary.Foxl2.Midage@meta.data, 
                                           cell_type_col = "celltype.level2", 
                                           label_col = "Group",
                                           n_threads = 3)

save(augur.ovary.Foxl2.Midage, file = paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_Midage_Augur_object.RData"))

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_Midage_Augur_UMAP.pdf"), height = 5, width = 6)
plot_umap(augur.ovary.Foxl2.Midage, ovary.Foxl2.Midage, cell_type_col = "celltype.level2")
dev.off()

pdf(paste0(Sys.Date(),"_10x_ovary_Benayoun_lab_Foxl2_Midage_Augur_Lollipop.pdf"))
plot_lollipop(augur.ovary.Foxl2.Midage)
dev.off()

################################################################################
# 3. Generate UMAPs 
################################################################################

# Aging

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_Aging_all_cells_AUGUR_UMAP_100dpi_.pdf"), height = 5, width = 6)
p <- plot_umap(augur.AC, ovary.AC, cell_type_col = "celltype.level2")
p + rasterize(geom_point(size = 0.05), dpi = 100) + 
  scale_color_gradientn(colors = c("darkorange", "purple"), limits = AUG_range)
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_Aging_all_cells_AUGUR_UMAP.pdf"), height = 5, width = 6)
p <- plot_umap(augur.AC, ovary.AC, cell_type_col = "celltype.level2")
p + scale_color_gradientn(colors = c("darkorange", "purple"), limits = AUG_range)
dev.off()


# VCD

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_VCD_all_cells_3m_30d_AUGUR_UMAP_100dpi_.pdf"), height = 5, width = 6)
p <- plot_umap(augur.ovary.VCD.30d.3m, ovary.VCD, cell_type_col = "celltype.level2")
p + rasterize(geom_point(size = 0.05), dpi = 100) + 
  scale_color_gradientn(colors = c("darkorange", "purple"), limits = AUG_range)
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_VCD_all_cells_3m_90d_AUGUR_UMAP_100dpi_.pdf"), height = 5, width = 6)
p <- plot_umap(augur.ovary.VCD.90d.3m, ovary.VCD, cell_type_col = "celltype.level2")
p + rasterize(geom_point(size = 0.05), dpi = 100) + 
  scale_color_gradientn(colors = c("darkorange", "purple"), limits = AUG_range)
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_VCD_all_cells_10m_30d_AUGUR_UMAP_100dpi_.pdf"), height = 5, width = 6)
p <- plot_umap(augur.ovary.VCD.30d.10m, ovary.VCD, cell_type_col = "celltype.level2")
p + rasterize(geom_point(size = 0.05), dpi = 100) + 
  scale_color_gradientn(colors = c("darkorange", "purple"), limits = AUG_range)
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_VCD_all_cells_10m_90d_AUGUR_UMAP_100dpi_.pdf"), height = 5, width = 6)
p <- plot_umap(augur.ovary.VCD.90d.10m, ovary.VCD, cell_type_col = "celltype.level2")
p + rasterize(geom_point(size = 0.05), dpi = 100) + 
  scale_color_gradientn(colors = c("darkorange", "purple"), limits = AUG_range)
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_VCD_all_cells_10m_30d_AUGUR_UMAP_100dpi_.pdf"), height = 5, width = 6)
p <- plot_umap(augur.ovary.VCD.30d.10m, ovary.VCD, cell_type_col = "celltype.level2")
p + rasterize(geom_point(size = 0.05), dpi = 100) + 
  scale_color_gradientn(colors = c("darkorange", "purple"), limits = AUG_range)
dev.off()

# Foxl2

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_Foxl2_all_cells_young_AUGUR_UMAP_100dpi_.pdf"), height = 5, width = 6)
p <- plot_umap(augur.ovary.Foxl2.Young, ovary.Foxl2, cell_type_col = "celltype.level2")
p + rasterize(geom_point(size = 0.1), dpi = 100) + 
  scale_color_gradientn(colors = c("darkorange", "purple"), limits = AUG_range)
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_Foxl2_all_cells_mid-age_AUGUR_UMAP_100dpi_.pdf"), height = 5, width = 6)
p <- plot_umap(augur.ovary.Foxl2.Midage, ovary.Foxl2, cell_type_col = "celltype.level2")
p + rasterize(geom_point(size = 0.05), dpi = 100) + 
  scale_color_gradientn(colors = c("darkorange", "purple"), limits = AUG_range)
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_Foxl2_all_cells_young_AUGUR_UMAP.pdf"), height = 5, width = 6)
p <- plot_umap(augur.ovary.Foxl2.Young, ovary.Foxl2, cell_type_col = "celltype.level2")
p + scale_color_gradientn(colors = c("darkorange", "purple"), limits = AUG_range)
dev.off()

pdf(paste0(Sys.Date(), "_10x_ovary_Benayoun_lab_Foxl2_all_cells_mid-age_AUGUR_UMAP.pdf"), height = 5, width = 6)
p <- plot_umap(augur.ovary.Foxl2.Midage, ovary.Foxl2, cell_type_col = "celltype.level2")
p + scale_color_gradientn(colors = c("darkorange", "purple"), limits = AUG_range)
dev.off()

################################################################################
# 4. Generate combined scatter plots - VCD & Foxl2 haploinsufficiency
################################################################################

# Define a fixed AUC range for all plots
AUG_range <- c(0.5, 0.9) 

#######################
# VCD
#######################

# Combine AUC results
auc_data.VCD.30d <- full_join(
  augur.ovary.VCD.30d.3m$AUC %>% rename(AUC_VCD_3m_30d = auc),
  augur.ovary.VCD.30d.10m$AUC %>% rename(AUC_VCD_10m_30d = auc),
  by = "cell_type"
)

auc_data.VCD.90d <- full_join(
  augur.ovary.VCD.90d.3m$AUC %>% rename(AUC_VCD_3m_90d = auc),
  augur.ovary.VCD.90d.10m$AUC %>% rename(AUC_VCD_10m_90d = auc),
  by = "cell_type"
)

# Generate plot
p.VCD.30d <- ggplot(auc_data.VCD.30d, aes(x = AUC_VCD_3m_30d, y = AUC_VCD_10m_30d, label = cell_type)) +
  geom_point(color = "gold", size = 6, shape = 17) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") +  
  geom_text(vjust = -0.5, size = 3) + 
  theme_minimal() +
  labs(title = "AUC Comparison: CTL vs VCD - 30d",
       x = "AUC (3m)",
       y = "AUC (10m)") +
  coord_fixed(ratio = 1, xlim = c(0.5, 0.9), ylim = c(0.5, 0.9))

p.VCD.90d <- ggplot(auc_data.VCD.90d, aes(x = AUC_VCD_3m_90d, y = AUC_VCD_10m_90d, label = cell_type)) +
  geom_point(color = "gold", size = 6, shape = 17) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") +  
  geom_text(vjust = -0.5, size = 3) + 
  theme_minimal() +
  labs(title = "AUC Comparison: CTL vs VCD - 90d",
       x = "AUC (3m)",
       y = "AUC (10m)") +
  coord_fixed(ratio = 1, xlim = c(0.5, 0.9), ylim = c(0.5, 0.9))

#######################
# Foxl2 haploinsufficiency
#######################

# Combine AUC results
auc_data.Foxl2 <- full_join(
  augur.ovary.Foxl2.Young$AUC %>% rename(AUC_Young = auc),
  augur.ovary.Foxl2.Midage$AUC %>% rename(AUC_Midage = auc),
  by = "cell_type"
)

# Generate plot
p.Foxl2 <- ggplot(auc_data.Foxl2, aes(x = AUC_Young, y = AUC_Midage, label = cell_type)) +
  geom_point(color = "springgreen", size = 6, shape = 18) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") +  
  geom_text(vjust = -0.5, size = 3) + 
  theme_minimal() +
  labs(title = "AUC Comparison: Young vs Mid-age",
       x = "AUC (Young)",
       y = "AUC (Mid-age)") +
  coord_fixed(ratio = 1, xlim = c(0.5, 0.9), ylim = c(0.5, 0.9))


# Save
pdf(paste0(Sys.Date(), "_MeMo_Augur_scatterplots_combined.pdf"), width = 20, height = 5) 
grid.arrange(p.VCD.30d, p.VCD.90d, p.Foxl2, ncol = 3) 
dev.off()

################################################################################
# 4-1. Generate combined scatter plots - VCD & Foxl2 haploinsufficiency
#      Treat NAs as 0.5
################################################################################

# VCD 
auc_data.VCD.30d <- full_join(
  augur.ovary.VCD.30d.3m$AUC %>% rename(AUC_VCD_3m_30d = auc),
  augur.ovary.VCD.30d.10m$AUC %>% rename(AUC_VCD_10m_30d = auc),
  by = "cell_type"
) %>%
  mutate(across(starts_with("AUC_VCD"), ~ replace_na(.x, 0.5)))  

auc_data.VCD.90d <- full_join(
  augur.ovary.VCD.90d.3m$AUC %>% rename(AUC_VCD_3m_90d = auc),
  augur.ovary.VCD.90d.10m$AUC %>% rename(AUC_VCD_10m_90d = auc),
  by = "cell_type"
) %>%
  mutate(across(starts_with("AUC_VCD"), ~ replace_na(.x, 0.5)))  

# Generate scatter plots
p.VCD.30d <- ggplot(auc_data.VCD.30d, aes(x = AUC_VCD_3m_30d, y = AUC_VCD_10m_30d, label = cell_type)) +
  geom_point(color = "gold", size = 6, shape = 17) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") +  
  geom_text(vjust = -0.5, size = 3) + 
  theme_minimal() +
  labs(title = "AUC Comparison: CTL vs VCD - 30d",
       x = "AUC (3m)",
       y = "AUC (10m)") +
  coord_fixed(ratio = 1, xlim = c(0.5, 0.9), ylim = c(0.5, 0.9))

p.VCD.90d <- ggplot(auc_data.VCD.90d, aes(x = AUC_VCD_3m_90d, y = AUC_VCD_10m_90d, label = cell_type)) +
  geom_point(color = "gold", size = 6, shape = 17) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") +  
  geom_text(vjust = -0.5, size = 3) + 
  theme_minimal() +
  labs(title = "AUC Comparison: CTL vs VCD - 90d",
       x = "AUC (3m)",
       y = "AUC (10m)") +
  coord_fixed(ratio = 1, xlim = c(0.5, 0.9), ylim = c(0.5, 0.9))

# Foxl2
auc_data.Foxl2 <- full_join(
  augur.ovary.Foxl2.Young$AUC %>% rename(AUC_Young = auc),
  augur.ovary.Foxl2.Midage$AUC %>% rename(AUC_Midage = auc),
  by = "cell_type"
) %>%
  mutate(across(starts_with("AUC"), ~ replace_na(.x, 0.5)))

p.Foxl2 <- ggplot(auc_data.Foxl2, aes(x = AUC_Young, y = AUC_Midage, label = cell_type)) +
  geom_point(color = "springgreen", size = 6, shape = 18) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") +  
  geom_text(vjust = -0.5, size = 3) + 
  theme_minimal() +
  labs(title = "AUC Comparison: Young vs Mid-age",
       x = "AUC (Young)",
       y = "AUC (Mid-age)") +
  coord_fixed(ratio = 1, xlim = c(0.5, 0.9), ylim = c(0.5, 0.9))

# Save 
pdf(paste0(Sys.Date(), "_MeMo_Augur_scatterplots_combined_NAs_as_0.5.pdf"), width = 20, height = 5)  
grid.arrange(p.VCD.30d, p.VCD.90d, p.Foxl2, ncol = 3)  
dev.off()

################################################################################################################################################################
sink(file = paste(Sys.Date(),"_Run_Augur_combined_session_Info.txt", sep =""))
sessionInfo()
sink()
