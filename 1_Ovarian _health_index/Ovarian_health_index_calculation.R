library(readxl)
library(dplyr)
library(beeswarm)

rm(list = ls())

###############################
# Menopause-model project
# Calculate ovarian health index for AC, VCD & Foxl2 haploinsufficiency models
###############################

###############################
# 1. Import data
###############################

# Data

AC.follicle    <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/Hormone_histology_OHI/Input_data/2025-02-27_MeMo_AC_follicle_data.txt", header = TRUE)
VCD.follicle   <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/Hormone_histology_OHI/Input_data/2025-02-27_MeMo_VCD_follicle_data.txt", header = TRUE)
Foxl2.follicle <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/Hormone_histology_OHI/Input_data/2025-02-27_MeMo_Foxl2_follicle_data.txt", header = TRUE)

AC.hormone    <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/Hormone_histology_OHI/Input_data/2025-02-27_MeMo_AC_hormone_data.txt", header = TRUE)
VCD.hormone   <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/Hormone_histology_OHI/Input_data/2025-02-27_MeMo_VCD_hormone_data.txt", header = TRUE)
Foxl2.hormone <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/Hormone_histology_OHI/Input_data/2025-02-27_MeMo_Foxl2_hormone_data.txt", header = TRUE)

# Metadata

VCD.metadata   <- read_excel("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/Animal_metadata/VCD_animals_metadata.xlsx")
Foxl2.metadata <- read_excel("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/Animal_metadata/Foxl2_haplo_animals_metadata.xlsx")

###############################
# 2. Process data
###############################

# Aging model

AC.follicle.cl <- AC.follicle[AC.follicle$OHI == "Y",]
AC.hormone.cl <- AC.hormone[AC.hormone$OHI == "Y",]

# VCD model

VCD.hormone.5M  <- VCD.hormone[VCD.hormone$State == "Post_I_5M",]
VCD.follicle.5M <- VCD.follicle[VCD.follicle$Mouse_ID %in% VCD.hormone.5M$Mouse_ID,]

rownames(VCD.hormone.5M)  <- VCD.hormone.5M$Mouse_ID
rownames(VCD.follicle.5M) <- VCD.follicle.5M$Mouse_ID

common.VCD <- intersect(rownames(VCD.hormone.5M), rownames(VCD.follicle.5M))     # 1 less sample in follicle counts due to processing error

VCD.hormone.5M.cl  <- VCD.hormone.5M[rownames(VCD.hormone.5M) %in% common.VCD,]
VCD.follicle.5M.cl <- VCD.follicle.5M[rownames(VCD.follicle.5M) %in% common.VCD,]

identical(rownames(VCD.hormone.5M.cl), rownames(VCD.follicle.5M.cl))             # TRUE

head(VCD.hormone.5M.cl)
head(VCD.follicle.5M.cl)

VCD.data <- cbind(VCD.follicle.5M.cl, VCD.hormone.5M.cl[,c("Age_at_injection", "Treatment", "AMH", "FSH", "LH", "INHBA")])

VCD.data$follicles_combined <- rowSums(VCD.data[, c("Primordial", "Primary", "Secondary", "Antral", "CL")], na.rm=TRUE)

VCD.data <- VCD.data %>%
  mutate(
    Treatment = ifelse(Treatment == "Safflower_oil", "CTL", Treatment),
    Group = paste(Treatment, Age_at_injection, sep = "_")
  )

# Foxl2 haploinsufficiency model

Foxl2.hormone.N <- Foxl2.hormone[Foxl2.hormone$Issues == "N",]

Foxl2.hormone.cl  <- Foxl2.hormone.N[Foxl2.hormone.N$OHI == "Y",]
Foxl2.follicle.cl <- Foxl2.follicle[Foxl2.follicle$OHI == "Y",]

identical(Foxl2.hormone.cl$Mouse_ID, Foxl2.follicle.cl$Mouse_ID)             # TRUE

Foxl2.data <- Foxl2.follicle.cl %>%
  inner_join(Foxl2.hormone.cl, by = "Mouse_ID") %>%
  select(Mouse_ID, 
         Primordial, Primary, Secondary, Antral, CL, 
         Genotype, Age_group, AMH, FSH, LH, INHBA)

Foxl2.data$follicles_combined <- rowSums(Foxl2.data[, c("Primordial", "Primary", "Secondary", "Antral", "CL")], na.rm=TRUE)

Foxl2.data <- Foxl2.data %>%
  mutate(Group = paste(Genotype, Age_group, sep = "_"))

write.table(VCD.data, file = paste0(Sys.Date(), "_MeMo_VCD_data_clean_for_ovarian_health_index_calculation.txt"), quote = FALSE, sep = "\t")
write.table(Foxl2.data, file = paste0(Sys.Date(), "_MeMo_Foxl2_Haplo_data_clean_for_ovarian_health_index_calculation.txt"), quote = FALSE, sep = "\t")

###############################
# 3. Calculate ovarian health index
###############################

# Assign scores based on AC data

AC.follicle.cl$follicles_combined <- rowSums(AC.follicle.cl[, c("Primordial", "Primary", "Secondary", "Antral", "CL")], na.rm=TRUE)

AC.data <- cbind(AC.hormone.cl, AC.follicle.cl[, c("follicles_combined")])
colnames(AC.data)[ncol(AC.data)] <- "follicles_combined"

# Calculate medians for YF and OF groups
medians <- AC.data %>%
  group_by(Group) %>%
  summarise(across(c(AMH, FSH, INHBA, follicles_combined), median, .names = "{.col}_median"))

YF_medians <- filter(medians, Group == "YF")
OF_medians <- filter(medians, Group == "OF")

# Function to assign scores
calculate_scores <- function(value, YF_median, OF_median, increasing) {
  if (increasing) {
    return(case_when(
      value <= YF_median ~ 3,
      value > YF_median & value < OF_median ~ 2,
      value >= OF_median ~ 1,
      TRUE ~ NA_real_
    ))
  } else {
    return(case_when(
      value <= OF_median ~ 1,
      value > OF_median & value < YF_median ~ 2,
      value >= YF_median ~ 3,
      TRUE ~ NA_real_
    ))
  }
}

# Function to process datasets and compute scores
process_data <- function(data) {
  data <- data %>%
    rowwise() %>%
    mutate(
      AMH_score = calculate_scores(AMH, YF_medians$AMH_median, OF_medians$AMH_median, FALSE),
      FSH_score = calculate_scores(FSH, YF_medians$FSH_median, OF_medians$FSH_median, TRUE),
      INHBA_score = calculate_scores(INHBA, YF_medians$INHBA_median, OF_medians$INHBA_median, FALSE),
      hormone_avg_score = (AMH_score + FSH_score + INHBA_score) / 3,
      follicle_score = calculate_scores(follicles_combined, YF_medians$follicles_combined_median, OF_medians$follicles_combined_median, FALSE),
      total_score = (hormone_avg_score + follicle_score) / 6 * 100
    )
  return(data)
}

# Process all datasets
AC.data <- process_data(AC.data)
VCD.data <- process_data(VCD.data)
Foxl2.data <- process_data(Foxl2.data)

# Save processed data
write.table(AC.data, "MeMo_Ovarian_health_index_AC_cohort.txt", quote = FALSE, sep = "\t")
write.table(VCD.data, "MeMo_Ovarian_health_index_VCD_cohort.txt", quote = FALSE, sep = "\t")
write.table(Foxl2.data, "MeMo_Ovarian_health_index_Foxl2_haplo_cohort.txt", quote = FALSE, sep = "\t")

###############################
# 4. Generate plots
###############################

AC.data$Group    <- factor(AC.data$Group, levels = c("YF", "OF"))
VCD.data$Group   <- factor(VCD.data$Group, levels = c("CTL_3M", "VCD_3M", "CTL_6M", "VCD_6M", "CTL_8M", "VCD_8M", "CTL_10M", "VCD_10M"))
Foxl2.data$Group <- factor(Foxl2.data$Group, levels = c("wt_Y", "het_Y", "wt_O", "het_O"))

# AC data

p_value.AC <- wilcox.test(AC.data$total_score[AC.data$Group == "YF"], AC.data$total_score[AC.data$Group == "OF"])$p.value

pdf(paste0(Sys.Date(), "_MeMo_AC_data_ovarian_health_index.pdf"), width = 8, height = 7)
boxplot(total_score ~ Group, AC.data, 
        outline = FALSE, ylim = c(0, 120), 
        col = c("deeppink", "deeppink4"), las = 1, 
        ylab = "Ovarian health index")
beeswarm(total_score ~ Group, AC.data, pch = 1, col = "black", add = TRUE, cex = 1)
text(1.5, 120, paste("p-value =", format(p_value.AC, digits = 5)))
dev.off()

# VCD data

p_value.VCD.3m <- wilcox.test(VCD.data$total_score[VCD.data$Group == "CTL_3M"], VCD.data$total_score[VCD.data$Group == "VCD_3M"])$p.value   
p_value.VCD.6m <- wilcox.test(VCD.data$total_score[VCD.data$Group == "CTL_6M"], VCD.data$total_score[VCD.data$Group == "VCD_6M"])$p.value      
p_value.VCD.8m <- wilcox.test(VCD.data$total_score[VCD.data$Group == "CTL_8M"], VCD.data$total_score[VCD.data$Group == "VCD_8M"])$p.value    
p_value.VCD.10m <- wilcox.test(VCD.data$total_score[VCD.data$Group == "CTL_10M"], VCD.data$total_score[VCD.data$Group == "VCD_10M"])$p.value 

pdf(paste0(Sys.Date(), "_MeMo_VCD_data_ovarian_health_index.pdf"), width = 8, height = 7)
boxplot(total_score ~ Group, VCD.data, 
        outline = FALSE, ylim = c(0, 120), 
        col = c("deeppink", "yellow"), las = 1, 
        ylab = "Ovarian health index")
beeswarm(total_score ~ Group, VCD.data, pch = 2, col = "black", add = TRUE, cex = 1)
text(1.5, 110, paste("p-value =", format(p_value.VCD.3m, digits = 5)))
text(3.5, 120, paste("p-value =", format(p_value.VCD.6m, digits = 5)))
text(5.5, 110, paste("p-value =", format(p_value.VCD.8m, digits = 5)))
text(7.5, 120, paste("p-value =", format(p_value.VCD.10m, digits = 5)))
dev.off()

# Foxl2 data

p_value.Foxl2.Y <- wilcox.test(Foxl2.data$total_score[Foxl2.data$Group == "wt_Y"], Foxl2.data$total_score[Foxl2.data$Group == "het_Y"])$p.value   
p_value.Foxl2.O <- wilcox.test(Foxl2.data$total_score[Foxl2.data$Group == "wt_O"], Foxl2.data$total_score[Foxl2.data$Group == "het_O"])$p.value      

pdf(paste0(Sys.Date(), "_MeMo_Foxl2_haplo_data_ovarian_health_index.pdf"), width = 8, height = 7)
boxplot(total_score ~ Group, Foxl2.data, 
        outline = FALSE, ylim = c(0, 120), 
        col = c("deeppink", "springgreen"), las = 1, 
        ylab = "Ovarian health index")
beeswarm(total_score ~ Group, Foxl2.data, pch = 5, col = "black", add = TRUE, cex = 1)
text(1.5, 120, paste("p-value =", format(p_value.Foxl2.Y, digits = 5)))
text(3.5, 120, paste("p-value =", format(p_value.Foxl2.O, digits = 5)))
dev.off()

###############################
# Save session info
sink(file = paste(Sys.Date(),"Ovarian_health_index_calculation_session_Info.txt", sep =""))
sessionInfo()
sink()
