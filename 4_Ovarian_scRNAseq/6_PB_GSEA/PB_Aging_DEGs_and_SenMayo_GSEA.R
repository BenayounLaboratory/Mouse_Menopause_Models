setwd("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/5_PB_GSEA")
options(stringsAsFactors = FALSE)

library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(stringr)

theme_set(theme_bw())
set.seed(123123)

rm(list = ls())

################################################################################
# 1. Load DESeq2 data and filter to shared cell types
################################################################################

Aging.DESeq2 <- get(load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/5_PB_GSEA/Input/2025-01-14_10x_ovary_Benayoun_lab_AC_PB_DESeq2_object.RData"))
VCD.DESeq2 <- get(load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/5_PB_GSEA/Input/2025-01-15_10x_ovary_Benayoun_lab_VCD_PB_DESeq2_object.RData"))
Foxl2.DESeq2 <- get(load("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/5_PB_GSEA/Input/2025-01-15_10x_ovary_Benayoun_lab_Foxl2_NULL_PB_DESeq2_object.RData"))

common_celltypes <- Reduce(unique, list(names(Aging.DESeq2), names(VCD.DESeq2), names(Foxl2.DESeq2)))

Aging <- Aging.DESeq2[common_celltypes]
VCD <- VCD.DESeq2[common_celltypes]
Foxl2 <- Foxl2.DESeq2[common_celltypes]

################################################################################
# 2. Aging DEG GSEA (AGE_UP and AGE_DOWN based on Aging)
################################################################################

positive_lfc <- lapply(Aging, function(x) rownames(x[x$padj < 0.05 & x$log2FoldChange > 0, ]))
negative_lfc <- lapply(Aging, function(x) rownames(x[x$padj < 0.05 & x$log2FoldChange < 0, ]))

deg_plot_data <- data.frame()

for (model in c("Foxl2", "VCD")) {
  for (ct in names(datasets[[model]])) {
    
    df <- data.frame(datasets[[model]][[ct]])
    df$Gene <- rownames(df)
    
    if (!"stat" %in% colnames(df) || all(is.na(df$stat))) {
      next
    }
    gene_list <- sort(setNames(df$stat, df$Gene), decreasing = TRUE)
    
    TERM2GENE_UP <- data.frame(TERM = "AGE_UP", GENE = positive_lfc[[ct]])
    TERM2GENE_DOWN <- data.frame(TERM = "AGE_DOWN", GENE = negative_lfc[[ct]])
    
    for (term in c("AGE_UP", "AGE_DOWN")) {
      term_genes <- if (term == "AGE_UP") TERM2GENE_UP else TERM2GENE_DOWN
      res <- tryCatch({
        GSEA(gene_list, TERM2GENE = term_genes, minGSSize = 5, maxGSSize = 600,
             eps = 1e-100, pvalueCutoff = 1, pAdjustMethod = 'BH', by = 'fgsea')
      }, error = function(e) NULL)
      
      if (!is.null(res) && term %in% res@result$ID) {
        r <- res@result[res@result$ID == term, ]
        deg_plot_data <- rbind(deg_plot_data, data.frame(
          Dataset = model, CellType = ct, Term = term,
          NES = r$NES, logFDR = -log10(r$p.adjust)
        ))
      }
    }
  }
}

# Bubble plot - Aging DEGs
deg_plot_data$Dataset <- factor(deg_plot_data$Dataset, levels = c("VCD", "Foxl2"))
deg_plot_data$Term <- factor(deg_plot_data$Term, levels = c("AGE_UP", "AGE_DOWN"))

pdf(paste0(Sys.Date(), "_Aging_DEGs_GSEA_BubblePlot.pdf"))
ggplot(deg_plot_data, aes(x = Dataset, y = CellType, size = logFDR, color = NES)) +
  geom_point(na.rm = TRUE) +
  scale_color_gradientn(colors = c("darkblue", "dodgerblue4", "dodgerblue3", "dodgerblue1",
                                   "white", "lightcoral", "brown1", "firebrick2", "firebrick4"),
                        limits = c(-1, 1), oob = scales::squish, na.value = "gray90") +
  facet_wrap(~Term, nrow = 1) +
  theme_minimal() +
  labs(title = "GSEA Bubble Plot: Aging DEGs in VCD & Foxl2", x = "Dataset", y = "Cell Type",
       size = "-log10(FDR)", color = "NES") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

################################################################################
# 3. SenMayo GSEA
################################################################################

Sym.SEN.MAYO <- read.gmt("~/Dropbox/Benayoun_lab/Projects/Menopause_model_project/Data/For_Github/4_Ovarian_scRNAseq/5_PB_GSEA/Input/SAUL_SEN_MAYO.v2023.2.Mm.gmt")
TERM2GENE_SEN <- data.frame(TERM = "SenMayo", GENE = Sym.SEN.MAYO$gene)

datasets <- list(Aging = Aging, VCD = VCD, Foxl2 = Foxl2)
senmayo_plot_data <- data.frame()

for (model in names(datasets)) {
  for (ct in names(datasets[[model]])) {
    df <- data.frame(datasets[[model]][[ct]])
    df$Gene <- rownames(df)
    
    if (!"stat" %in% colnames(df) || all(is.na(df$stat))) {
      next
    }
    
    gene_list <- sort(setNames(df$stat, df$Gene), decreasing = TRUE)
    
    res <- tryCatch({
      GSEA(gene_list, TERM2GENE = TERM2GENE_SEN, minGSSize = 15, maxGSSize = 500,
           eps = 1e-100, pvalueCutoff = 1, pAdjustMethod = 'BH', by = 'fgsea')
    }, error = function(e) NULL)
    
    if (!is.null(res) && "SenMayo" %in% res@result$ID) {
      sen <- res@result[res@result$ID == "SenMayo", ]
      senmayo_plot_data <- rbind(senmayo_plot_data, data.frame(
        Dataset = model, CellType = ct, NES = sen$NES, logFDR = -log10(sen$p.adjust)
      ))
    }
  }
}

# Bubble plot - SenMayo
senmayo_plot_data$Dataset <- factor(senmayo_plot_data$Dataset, levels = c("Aging", "VCD", "Foxl2"))

pdf(paste0(Sys.Date(), "_SenMayo_GSEA_BubblePlot.pdf"))
ggplot(senmayo_plot_data, aes(x = Dataset, y = CellType, size = logFDR, color = NES)) +
  geom_point(na.rm = TRUE) +
  scale_color_gradientn(
    colors = c("darkblue", "dodgerblue4", "dodgerblue3", "dodgerblue1",
               "white", "lightcoral", "brown1", "firebrick2", "firebrick4"),
    limits = c(-1, 1), oob = scales::squish, na.value = "gray90") + 
  theme_minimal() +
  labs(title = "SenMayo GSEA Bubble Plot", x = "Dataset", y = "Cell Type",
       size = "-log10(FDR)", color = "NES") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

################################################################################
sink(file = paste0(Sys.Date(), "_PB_Aging_DEGs_and_SenMayo_GSEA_sessionInfo.txt"))
sessionInfo()
sink()    
