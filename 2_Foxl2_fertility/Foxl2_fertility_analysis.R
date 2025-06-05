options(stringsAsFactors = FALSE)

library(dplyr)
library(survival)
library(survminer)
library(beeswarm)
library(tidyverse)
library(outliers)

###############################################
# Foxl2 haplo model fertility data - Pup count & Latency
# Kaplan-Meier Survival Curves by Genotype
###############################################

###############################################
# 1. Import data
###############################################

Foxl2.haplo.fertility.data <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_model_project/Foxl2_haplo_fertility_data/Input_data/Foxl2_haplo_fertility_data_filtered.txt",
                                         header = TRUE, sep = "\t")

# head(Foxl2.haplo.fertility.data)

Foxl2.haplo.fertility.data.cl <- Foxl2.haplo.fertility.data[Foxl2.haplo.fertility.data$Comments == "N",]

###############################################
# 2. Generate boxplots - pup counts
###############################################

# Ensure Genotype is a factor
Foxl2.haplo.fertility.data.cl$Genotype <- factor(Foxl2.haplo.fertility.data.cl$Genotype, levels = c("Wt", "Het"))

# Perform Wilcoxon test
wilcox_test <- wilcox.test(Foxl2.haplo.fertility.data.cl$Num_pup[Foxl2.haplo.fertility.data.cl$Genotype == "Wt"], Foxl2.haplo.fertility.data.cl$Num_pup[Foxl2.haplo.fertility.data$Genotype == "Het"])     # p-value = 0.224

# Generate a boxplot for all litter sizes combined
pdf(paste(Sys.Date(), "MeMo_Foxl2_haplo_fertility_1st_litter_filtered_data.pdf", sep = "_"), width = 7, height = 7)
boxplot(Num_pup ~ Genotype, data = Foxl2.haplo.fertility.data.cl, 
        outline = TRUE, 
        col = c("deeppink", "springgreen"),
        las = 1, ylim = c(0,15),
        main = "Foxl2 Haplo Fertility - 1st litter (filtered data)", 
        ylab = "1st litter size (count)", 
        xlab = "Genotype")
beeswarm(Num_pup ~ Genotype, data = Foxl2.haplo.fertility.data.cl, pch = 18, col = "black", add = TRUE, cex = 1.0)
text(1.5, 15, paste("p =", format(wilcox_test$p.value, digits = 4, nsmall = 4)))
dev.off()

# Outlier test
grubbs.test(Foxl2.haplo.fertility.data.cl$Num_pup[Foxl2.haplo.fertility.data.cl$Genotype == "Wt"])     # p-value = 0.0002496
grubbs.test(Foxl2.haplo.fertility.data.cl$Num_pup[Foxl2.haplo.fertility.data.cl$Genotype == "Het"])    # p-value = 0.2019

# Remove outlier

Foxl2.haplo.fertility.data.cl.cl <- Foxl2.haplo.fertility.data.cl[Foxl2.haplo.fertility.data.cl$Breeder != "2110",]

# Perform Wilcoxon test
wilcox_test <- wilcox.test(Foxl2.haplo.fertility.data.cl.cl$Num_pup[Foxl2.haplo.fertility.data.cl.cl$Genotype == "Wt"], Foxl2.haplo.fertility.data.cl.cl$Num_pup[Foxl2.haplo.fertility.data$Genotype == "Het"])     # p-value = 0.07489

# Generate a boxplot for all litter sizes combined
pdf(paste(Sys.Date(), "MeMo_Foxl2_haplo_fertility_1st_litter_filtered_data_outlier_removed.pdf", sep = "_"), width = 7, height = 7)
boxplot(Num_pup ~ Genotype, data = Foxl2.haplo.fertility.data.cl.cl, 
        outline = TRUE, 
        col = c("deeppink", "springgreen"),
        las = 1, ylim = c(0,15),
        main = "Foxl2 Haplo Fertility - 1st litter (filtered data & outlier removed)", 
        ylab = "1st litter size (count)", 
        xlab = "Genotype")
beeswarm(Num_pup ~ Genotype, data = Foxl2.haplo.fertility.data.cl.cl, pch = 18, col = "black", add = TRUE, cex = 1.0)
text(1.5, 15, paste("p =", format(wilcox_test$p.value, digits = 4, nsmall = 4)))
dev.off()

###############################################
# 3. Generate latency curve
###############################################

# Add Event column
Foxl2.haplo.fertility.data.cl.cl <- Foxl2.haplo.fertility.data.cl.cl %>%
  mutate(Event = ifelse(is.na(Num_pup) | Num_pup == 0, 0, 1))

# Create a survival object for first litter
surv_object_1st <- Surv(time = Foxl2.haplo.fertility.data.cl.cl$Latency, event = Foxl2.haplo.fertility.data.cl.cl$Event)

# Fit a Kaplan-Meier survival curve by Genotype
fit_1st <- survfit(surv_object_1st ~ Genotype, data = Foxl2.haplo.fertility.data.cl.cl)

# Plot survival curve for first litter
pdf(paste(Sys.Date(), "MeMo_Foxl2_haplo_fertility_latency_1st_litter_plot.pdf", sep = "_"))
ggsurvplot(fit_1st, 
           data = Foxl2.haplo.fertility.data.cl.cl, 
           pval = TRUE,  
           conf.int = FALSE, 
           xlab = "Latency (days)", 
           ylab = "Nulliparous (%)", 
           ggtheme = theme_minimal(), 
           title = "Foxl2 Haplo Fertility Latency - First Litter Only", 
           legend.title = "Genotype", 
           legend.labs = c("wt", "het"), 
           palette = c("deeppink", "springgreen"))
dev.off()

################################################################################
sink(file = paste(Sys.Date(),"Foxl2_fertility_analysis_session_Info.txt",sep="_"))
sessionInfo()
sink()
