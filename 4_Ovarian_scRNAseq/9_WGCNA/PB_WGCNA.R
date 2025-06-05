setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/6_WGCNA/")

library(WGCNA)
library(pheatmap)
library(ggplot2)
library(dplyr)

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
library(ggplot2)

options(mc.cores = 3)

rm(list = ls())

################################################################################
# 1. Load data
################################################################################

# Load datasets
load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/6_WGCNA/Input/Aging_PB_VST_counts.RData")
VST_Aging <- vst.counts
rm(vst.counts)

load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/6_WGCNA/Input/VCD_PB_VST_counts.RData")
VST_VCD <- vst.counts
rm(vst.counts)

load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/6_WGCNA/Input/Foxl2_PB_VST_counts.RData")
VST_Foxl2 <- vst.counts
rm(vst.counts)

################################################################################
# 2. Define cell types and create directories
################################################################################

# Define cell types
cell_types <- intersect(names(VST_Aging), intersect(names(VST_VCD), names(VST_Foxl2)))
# [1] "BEC"        "DNT"        "Epithelial" "Granulosa"  "Stroma"     "Theca"     

# Create folders for each cell type
for (cell in cell_types) {
  dir.create(file.path(getwd(), cell), showWarnings = FALSE, recursive = TRUE)
}

################################################################################
# 3. Iterate over cell types and perform WGCNA
################################################################################

enriched_GO_list <- list()

for (cell in cell_types) {
  
  # Define output directory
  output_dir <- file.path(getwd(), cell)
  
  # Extract data
  Aging_data <- as.data.frame(VST_Aging[cell])
  VCD_data <- as.data.frame(VST_VCD[cell])
  Foxl2_data <- as.data.frame(VST_Foxl2[cell])
  
  # Find overlapping genes
  overlapping.genes <- intersect(rownames(Aging_data), intersect(rownames(VCD_data), rownames(Foxl2_data)))
  
  # Prepare expression matrix
  datExpr <- as.data.frame(t(Aging_data[overlapping.genes,]))
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  
  # Save cleaned data
  save(datExpr, file = file.path(output_dir, paste0(Sys.Date(), "_", cell, "_cleaned_WGCNA_input.RData")))
  
  ################################################################################
  # 3-1. WGCNA: Soft thresholding and network construction
  ################################################################################
  
  powers <- c(1:50)
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  softPower <- sft$powerEstimate
  
  adjacency_mat <- adjacency(datExpr, power = softPower)
  TOM <- TOMsimilarity(adjacency_mat)
  dissTOM <- 1 - TOM
  
  save(TOM, file = file.path(output_dir, paste0(Sys.Date(), "_", cell, "_TOM_Matrix.RData")))
  
  ################################################################################
  # 3-2. Identify modules
  ################################################################################
  
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  minModuleSize <- 30
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, minClusterSize = minModuleSize)
  moduleColors <- labels2colors(dynamicMods)
  
  save(moduleColors, file = file.path(output_dir, paste0(Sys.Date(), "_", cell, "_Module_Colors.RData")))
  
  # Generate Gene Dendrogram
  pdf(file = file.path(output_dir, paste0(Sys.Date(), "_", cell, "_Gene_Dendrogram.pdf")), width = 10, height = 6)
  plot(geneTree,
       main = paste("Gene Clustering Dendrogram -", cell),
       xlab = "", sub = "", labels = FALSE)
  dev.off()
  
  pdf(file = file.path(output_dir, paste0(Sys.Date(), "_", cell, "_Gene_Dendrogram_with_module_colors.pdf")), width = 10, height = 6)
  plotDendroAndColors(geneTree, 
                      colors = moduleColors, 
                      groupLabels = "Modules", 
                      main = paste("Gene Clustering Dendrogram -", cell),
                      cex.axis = 0.5, 
                      cex.lab = 0.7, 
                      dendroLabels = FALSE) 
  dev.off()
  
  ################################################################################
  # 3-3. Extract eigengenes and correlate with traits
  ################################################################################
  
  MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  sampleNames <- rownames(MEs)
  
  traitMatrix <- data.frame(
    Group1 = as.numeric(grepl("OF", sampleNames)), 
    Group2 = as.numeric(grepl("YF", sampleNames))
  )
  rownames(traitMatrix) <- sampleNames
  
  moduleTraitCor <- cor(MEs, traitMatrix, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
  
  save(moduleTraitCor, moduleTraitPvalue, 
       file = file.path(output_dir, paste0(Sys.Date(), "_", cell, "_Module_Trait_Correlations.RData")))
  
  ################################################################################
  # 3-4. Plot module-trait heatmap
  ################################################################################
  
  cluster_rows <- ifelse(nrow(moduleTraitCor) > 1, TRUE, FALSE)
  cluster_cols <- ifelse(ncol(moduleTraitCor) > 1, TRUE, FALSE)
  
  pdf(file = file.path(output_dir, paste0(Sys.Date(), "_", cell, "_Module_Trait_Correlation_Heatmap.pdf")))
  
  textMatrix <- paste0(round(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  pheatmap(moduleTraitCor, 
           color = colorRampPalette(c("blue", "white", "red"))(100),
           display_numbers = textMatrix, 
           cluster_rows = cluster_rows, 
           cluster_cols = cluster_cols, 
           main = paste(cell, "Module-Trait Correlations"))
  
  dev.off()
  
  ################################################################################
  # 3-5. Extract and analyze significant modules (pval < 0.1)
  ################################################################################
  
  significant_modules <- unique(rownames(which(moduleTraitPvalue < 0.1, arr.ind = TRUE)))
  module_gene_list <- list()
  
  for (module in significant_modules) {
    
    module_clean <- sub("^ME", "", module) 
    module_genes <- names(datExpr)[moduleColors == module_clean]
    
    enriched_GO <- enrichGO(gene = module_genes, 
                            OrgDb = org.Mm.eg.db, 
                            keyType = "SYMBOL", 
                            universe = overlapping.genes,
                            ont = "ALL", 
                            pvalueCutoff = 0.1)
    
    module_gene_list[[module]] <- module_genes
    
    enriched_GO_list[[cell]][[module]] <- enriched_GO
    
  }
  
  ################################################################################
  # 3-6. Save results
  ################################################################################
  
  save(module_gene_list, 
       file = file.path(output_dir, paste0(Sys.Date(), "_", cell, "_Module_Gene_Lists.RData")))
  
  write.csv(significant_modules, 
            file = file.path(output_dir, paste0(Sys.Date(), "_", cell, "_Significant_Modules.csv")))

}

################################################################################
# 4. Generate combined dotplot - ORA
################################################################################

# Initialize list
all_go_terms_list <- list()

# Loop through each cell type
for (cell in names(enriched_GO_list)) {
  
  # Loop through each module
  for (module in names(enriched_GO_list[[cell]])) {
    
    enriched_result <- enriched_GO_list[[cell]][[module]]
    
    if (!is.null(enriched_result) && nrow(enriched_result@result) > 0) {
        go_data <- enriched_result@result %>%
        dplyr::select(ID, Description, p.adjust, Count) %>%
        mutate(CellType = cell, Module = module) 
      
      all_go_terms_list[[paste0(cell, "_", module)]] <- go_data
    }
  }
}

# Combine data
all_go_terms_df <- bind_rows(all_go_terms_list)

# Calculate -log10(padj)
all_go_terms_df$logP <- -log10(all_go_terms_df$p.adjust)

# Generate Bubble Plot
pdf(paste0(Sys.Date(), "_MeMo_WGCNA_ORA_bubbleplots.pdf"), width = 20, height = 30)
ggplot(all_go_terms_df, aes(x = Module, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) + 
  scale_size_continuous(range = c(3, 10)) + 
  scale_color_gradientn(colors = c("darkorange2", "purple"), limits = c(0, 0.1)) +
  facet_wrap(~ CellType, scales = "free_x") +  
  theme_minimal(base_size = 9) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(size = "Gene Count", color = "p.adjust", x = "Module", y = "GO Term", title = "GO Enrichment Bubble Plot")
dev.off()

################################################################################
# Save Session Info
sink(file = paste0(Sys.Date(), "_PB_WGCNA_session_info.txt"))
sessionInfo()
sink()
