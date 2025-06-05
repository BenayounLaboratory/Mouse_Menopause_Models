setwd("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/Shiny_app/")

library(Seurat)
library(ShinyCell)

###############################################################################
# Menopause-model project
# Shiny app - Aging, VCD and Foxl2 haploinsufficiency model scRNAseq 
###############################################################################

###############################################################################
# 1. Import data
###############################################################################

# Aging model

load("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/Shiny_app/2024-10-15_10x_ovary_Benayoun_lab_AC_Seurat_object_with_final_annotation.RData")

ovary.AC@meta.data$Model <- "Aging"

# VCD model

load("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/Shiny_app/2024-10-24_10x_ovary_Benayoun_lab_VCD_Seurat_object_with_final_annotation.RData")

ovary.VCD@meta.data$Model <- "VCD"

# Foxl2 haploinsufficiency model

load("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/Shiny_app/2024-10-24_10x_ovary_Benayoun_lab_Foxl2_Seurat_object_with_final_annotation.RData")

ovary.Foxl2@meta.data$Model <- "Foxl2"

###############################################################################
# 2. Create ShinyCell configuration and configure settings
###############################################################################

# Aging model

scConf1 <- createConfig(ovary.AC)
scConf1 <- delMeta(scConf1, c("orig.ident", "SCT_snn_res.2", "DoubletFinder", "scds_hybrid", "DoubletCall", "PtprcPos_UCell", "is.pure", "is.pure.level1"))
scConf1 <- modColours(scConf1, meta.to.mod <- "Library", 
                      new.colours<- c("deeppink", "deeppink", "deeppink4", "deeppink4"))
scConf1 <- modDefault(scConf1, default1 = "celltype.level2", default2 = "celltype.level1")
makeShinyFiles(ovary.AC, scConf1, gex.assay <- "RNA", gex.slot <- "data",
               gene.mapping <- TRUE, shiny.prefix <- "sc1",
               shiny.dir <- "shinyAppMulti/",
               default.gene1 <- "Cyp19a1", default.gene2 <- "Cyp11a1",
               default.multigene <- c("Amh", "Cyp19a1", "Cyp17a1", "Cyp11a1", "Pdgfra", "Col1a1",       
                                      "Acta2", "Notch3", "Rgs5", "Flt1",  "Prox1", "Upk1b"),
               default.dimred <- c("UMAP_1", "UMAP_2")) 

# VCD model

DefaultAssay(ovary.VCD) <- "SCT"

scConf2 <- createConfig(ovary.VCD)
scConf2 <- delMeta(scConf2, c("orig.ident", "DoubletFinder", "scds_hybrid", "DoubletCall", "immune_UCell", "is.pure", "is.pure.level1"))
scConf2 <- modColours(scConf2, meta.to.mod <- "Library", 
                      new.colours<- c("#A40C5E","#A40C5E",
                                      "#970B57","#970B57",
                                      "#FF1493","#FF1493",
                                      "#CB0F75","#CB0F75",
                                      "#B18B05","#B18B05",
                                      "#8B6508","#8B6508",
                                      "#FFD700","#FFD700",
                                      "#D8B102","#D8B102"))
scConf2 <- modDefault(scConf2, default1 = "celltype.level2", default2 = "celltype.level1")
makeShinyFiles(ovary.VCD, scConf2, gex.assay <- "SCT", gex.slot <- "data",
               gene.mapping <- TRUE, shiny.prefix <- "sc2",
               shiny.dir <- "shinyAppMulti/",
               default.gene1 <- "Cyp19a1", default.gene2 <- "Cyp11a1",
               default.multigene <- c("Amh", "Cyp19a1", "Cyp17a1", "Cyp11a1", "Pdgfra", "Col1a1",       
                                      "Acta2", "Notch3", "Rgs5", "Flt1",  "Prox1", "Upk1b"),
               default.dimred <- c("UMAP_1", "UMAP_2")) 

# Foxl2 haploinsufficiency model

DefaultAssay(ovary.Foxl2) <- "SCT"

scConf3 <- createConfig(ovary.Foxl2)
scConf3 <- delMeta(scConf3, c("orig.ident", "DoubletFinder", "scds_hybrid", "DoubletCall", "immune_UCell", "is.pure", "is.pure.level1"))
scConf3 <- modColours(scConf3, meta.to.mod <- "Library", 
                      new.colours<- c("#508B8E","#508B8E","#508B8E",
                                      "#ADDFE4","#ADDFE4",
                                      "#B10D66","#B10D66","#B10D66",
                                      "#FF1493","#FF1493"))
scConf3 <- modDefault(scConf3, default1 = "celltype.level2", default2 = "celltype.level1")
makeShinyFiles(ovary.Foxl2, scConf3, gex.assay <- "RNA", gex.slot <- "data",
               gene.mapping <- TRUE, shiny.prefix <- "sc3",
               shiny.dir <- "shinyAppMulti/",
               default.gene1 <- "Cyp19a1", default.gene2 <- "Cyp11a1",
               default.multigene <- c("Amh", "Cyp19a1", "Cyp17a1", "Cyp11a1", "Pdgfra", "Col1a1",       
                                      "Acta2", "Notch3", "Rgs5", "Flt1",  "Prox1", "Upk1b"),
               default.dimred <- c("UMAP_1", "UMAP_2")) 

###############################################################################
# 3. Generate code for Shiny app
###############################################################################

citation = list(
  author  = "Kim M., Bhala R., Wang J. et al.",
  title   = "Systematic characterization of mouse menopause models reveals shared and model-specific signatures.")
makeShinyCodesMulti(
  shiny.title = "Mouse menopause models", shiny.footnotes = citation,
  shiny.prefix = c("sc1", "sc2", "sc3"),
  shiny.headers = c("Aging model", "VCD model", "Foxl2 haploinsufficiency model"), 
  shiny.dir = "shinyAppMulti/")

###############################################################################
sink(file = paste(Sys.Date(), "_MeMo_Shiny_app_Session_Info.txt"))
sessionInfo()
sink()
