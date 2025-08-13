# Author: Sunni Patton
# Date: 10/21/24 (last edited)
# Title: Differential Abundance Analysis 

## Set seed ====
set.seed(123)

## Load packages ====
library(dplyr)
library(ggplot2)
library(tidyr)
library(ANCOMBC)
library(phyloseq)
library(mia)
library(lme4)
library(lmerTest)

## Load data ====
ps.pruned <- readRDS(here::here("Output Files/04 - Phyloseq - Output/ps.pruned.RDS"))

## Compare blank over time ====
ps.pruned.blank <- subset_samples(ps.pruned, Treatment == "0")
## No treatment T0 - T4 
ps.pruned.blank.T0T4 <- subset_samples(ps.pruned.blank, Time_TotalDays == "0" | Time_TotalDays == "2" | Time_TotalDays == "4")
ps.pruned.blank.T0T4@sam_data$Time_TotalDays <- factor(ps.pruned.blank.T0T4@sam_data$Time_TotalDays, levels = c("0", "2", "4"))
## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.pruned.blank.T0T4)
## ANCOMBC2
### Doing pairwise within these 
output_blank_T0T4_pw <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Time_TotalDays", rand_formula = "(1|TankID)",  p_adj_method = "fdr", group = "Time_TotalDays",
                         alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = TRUE) ## mdfdr not needed for pairwise
                      # got the is.singular warning
                      # got >50 warnings

View(output_blank_T0T4_pw$res_pair)
  # nothing sig diff


### Doing dunnett test to compare against T0
output_blank_T0T4 <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Time_TotalDays", rand_formula = "(1|TankID)",  p_adj_method = "fdr", group = "Time_TotalDays",
                                 alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, dunnet = TRUE, 
                              mdfdr_control = list(fwer_ctrl_method = "fdr", B =100))
                    # got the is.singular warning

View(output_blank_T0T4$res_dunn)
  # nothing sig diff

# No treatment T9-T19
ps.pruned.blank.T9T14T19 <- subset_samples(ps.pruned.blank, Time_TotalDays == "9" | Time_TotalDays == "14" | Time_TotalDays == "19")
ps.pruned.blank.T9T14T19@sam_data$Time_TotalDays <- factor(ps.pruned.blank.T9T14T19@sam_data$Time_TotalDays, levels = c("9", "14", "19"))
## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.pruned.blank.T9T14T19)
## ANCOMBC2
### Doing pairwise within these 
output_blank_T9T14T19_pw <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Time_TotalDays", rand_formula = "(1|TankID)",  p_adj_method = "fdr", group = "Time_TotalDays",
                                 alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = TRUE) ## mdfdr not needed for pairwise
                            # got is.singular warning
view(output_blank_T9T14T19_pw$res_pair)
# nothing diff abundant 

### Doing dunnett test to compare against T9
output_blank_T9T14T19 <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Time_TotalDays", rand_formula = "(1|TankID)",  p_adj_method = "fdr", group = "Time_TotalDays",
                              alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, dunnet = TRUE, 
                              mdfdr_control = list(fwer_ctrl_method = "fdr", B =100))
                            # got the is.singular warning
view(output_blank_T9T14T19$res_dunn)
# nothing diff abundant 

# no treatment T24-T34
ps.pruned.blank.T24T29T34 <- subset_samples(ps.pruned.blank, Time_TotalDays == "24" | Time_TotalDays == "29" | Time_TotalDays == "34")
ps.pruned.blank.T24T29T34@sam_data$Time_TotalDays <- factor(ps.pruned.blank.T24T29T34@sam_data$Time_TotalDays, levels = c("24", "29", "34"))
## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.pruned.blank.T24T29T34)
## ANCOMBC2
### Doing pairwise within these 
output_blank_T24T29T34_pw <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Time_TotalDays", rand_formula = "(1|TankID)",  p_adj_method = "fdr", group = "Time_TotalDays",
                                     alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = TRUE) ## mdfdr not needed for pairwise
                                    # got is.singular warning
view(output_blank_T24T29T34_pw$res_pair)

## ANCOMBC formula
output_blank <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Time", rand_formula = "(1|TankID)",  p_adj_method = "fdr", group = "Time",
                         alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = TRUE)) ## mdfdr not needed for pairwise
View(output_blank$res)
## no sig taxa in T2 compared to T0
## no sig taxa in T4 compared to T0
## when i did pairwise, nothing significant for T0-T2, T0-T4, or T2-T4

## no sig taxa in T14 compared to T9
## no sig taxa in T19 compared to T9
## when i did pairwise, nothing significant for T14-T9, T19-T9, or T14-T19

## no sig taxa in T29 compared to T24
## 3 sig taxa in T34 compared T24, but they didn't pass sensitivity threshold
## when i did pairwise, nothing significant for T24-T29, T24-T34, or T29-T34

## but i am getting boundary singular on all of these
## but this could probably be enough justification for grouping these into three groups:
## antibiotics, temp ramp, temp


saveRDS(output_blank, here::here("Output Files/08 - Differential Abundance - Output/ANCOMBC_blank.rds"))

## several ASVs at different time points are differentially abundant compared to T0, further indicating that there is some temporal driver of microbiome diversity (which we've seen in the other metrics too)

### Subset only T0, T2, T4
ps.pruned@sam_data$groups <- paste0(ps.pruned@sam_data$Time_TotalDays)

ps.pruned@sam_data$groups[ps.pruned@sam_data$groups == "0"] <- "Antibiotic challenge"  
ps.pruned@sam_data$groups[ps.pruned@sam_data$groups == "2"] <- "Antibiotic challenge"  
ps.pruned@sam_data$groups[ps.pruned@sam_data$groups == "4"] <- "Antibiotic challenge"  
ps.pruned@sam_data$groups[ps.pruned@sam_data$groups == "9"] <- "Temperature ramp"  
ps.pruned@sam_data$groups[ps.pruned@sam_data$groups == "14"] <- "Temperature ramp"  
ps.pruned@sam_data$groups[ps.pruned@sam_data$groups == "19"] <- "Temperature challenge"  
ps.pruned@sam_data$groups[ps.pruned@sam_data$groups == "24"] <- "Temperature challenge"  
ps.pruned@sam_data$groups[ps.pruned@sam_data$groups == "29"] <- "Temperature challenge"  
ps.pruned@sam_data$groups[ps.pruned@sam_data$groups == "34"] <- "Temperature challenge"  

ps.pruned.T0T2T4 <- subset_samples(ps.pruned, Time_TotalDays == "0" |
                                         Time_TotalDays == "2" |
                                         Time_TotalDays == "4")
ps.pruned.T9T14T19 <- subset_samples(ps.pruned, Time_TotalDays == "9" | 
                                             Time_TotalDays == "14" | 
                                             Time_TotalDays == "19")
ps.pruned.T24T29T34 <- subset_samples(ps.pruned, Time_TotalDays == "24" | 
                                              Time_TotalDays == "29" | 
                                              Time_TotalDays == "34")

ps.pruned.T0T2T4@sam_data$Pretreatment <- factor(ps.pruned.T0T2T4@sam_data$Pretreatment, levels = c("No Treatment", "Antibiotics"))

tse <- mia::convertFromPhyloseq(ps.pruned.T0T2T4)

## ANCOMBC formula
output_T0T2T4 <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment", 
                          rand_formula = "(1|TankID)",  p_adj_method = "fdr",
                          group = "Pretreatment", alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE,
                          dunnet = TRUE,
                          mdfdr_control = list(fwer_ctrl_method = "fdr", B =100))

### still boundary is single 
### T0T2T4 has 16 abx 32 no treat
view(output_T0T2T4$res)
## ASV1, ASV2, ASV91 increased in antibiotic group compared to control (and each passed)
## go back make sure dunnett was done 

ps.pruned.T9T14T19@sam_data$Pretreatment <- factor(ps.pruned.T9T14T19@sam_data$Pretreatment, levels = c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"))
tse <- mia::convertFromPhyloseq(ps.pruned.T9T14T19)
output.T9T14T19 <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment", 
                          rand_formula = "(1|TankID)",  p_adj_method = "fdr",
                          group = "Pretreatment", alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE,
                          dunnet = TRUE,
                          mdfdr_control = list(fwer_ctrl_method = "fdr", B =100))

## all groups n = 12
view(output.T9T14T19$res)
view(output.T9T14T19$res_dunn)


## Compare Treatments at T0 ====
ps.pruned.T0 <- subset_samples(ps.pruned, Time == "0")

## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.pruned.T0)

## Set order for Treatment
tse$Treatment_Long <- factor(tse$Treatment_Long, levels = c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"))


## ANCOMBC formula
output_T0 <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Treatment_Long",  p_adj_method = "fdr", group = "Treatment_Long",
                         alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 
View(output_T0$res)
saveRDS(output_T0, here::here("Output Files/08 - Differential Abundance - Output/output_T0.rds"))

## no significant difference in relative abundance at T0 between 'treatments' which is good

## Compare treatments at T2 ====
ps.pruned.T2 <- subset_samples(ps.pruned, Time_TotalDays== "2")

## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.pruned.T2)

## Set order for Treatment
tse$Pretreatment <- factor(tse$Pretreatment, levels = c("No Treatment", "Antibiotics"))


## ANCOMBC formula
output_T2 <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment", p_adj_method = "fdr", group = "Pretreatment",
                      alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 
View(output_T2$res)
saveRDS(output_T2, here::here("Output Files/08 - Differential Abundance - Output/output_T2.rds"))


## Compare treatments at T4 ====
ps.pruned.T4 <- subset_samples(ps.pruned, Time_TotalDays == "4")

## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.pruned.T4)

## Set order for Treatment
tse$Pretreatment <- factor(tse$Pretreatment, levels = c("No Treatment", "Antibiotics"))


## ANCOMBC formula
output_T4 <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment", p_adj_method = "fdr", group = "Pretreatment",
                      alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 
View(output_T4$res)
saveRDS(output_T4, here::here("Output Files/08 - Differential Abundance - Output/output_T4.rds"))


## Compare treatments at T9 ====
ps.pruned.T9 <- subset_samples(ps.pruned, Time_TotalDays == "9")

## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.pruned.T9)

## Set order for Treatment
tse$Pretreatment <- factor(tse$Pretreatment, levels = c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"))

## ANCOMBC formula
output_T9 <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment", p_adj_method = "fdr", group = "Pretreatment",
                      alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 
View(output_T9$res)
saveRDS(output_T9, here::here("Output Files/08 - Differential Abundance - Output/output_T9.rds"))

## Compare treatments at T14 ====
ps.pruned.T14 <- subset_samples(ps.pruned, Time_TotalDays == "14")

## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.pruned.T14)

## Set order for Treatment
tse$Pretreatment <- factor(tse$Pretreatment, levels = c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"))


## ANCOMBC formula
output_T14 <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment", p_adj_method = "fdr", group = "Pretreatment",
                      alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 
View(output_T14$res)

saveRDS(output_T14, here::here("Output Files/08 - Differential Abundance - Output/output_T14.rds"))


## Compare treatments at T19 ====
ps.pruned.T19 <- subset_samples(ps.pruned, Time_TotalDays == "19")

## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.pruned.T19)

## Set order for Treatment
tse$Pretreatment <- factor(tse$Pretreatment, levels = c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"))


## ANCOMBC formula
output_T19 <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment", p_adj_method = "fdr", group = "Pretreatment",
                       alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 
View(output_T19$res)


## Compare treatments at T24 ====
ps.pruned.T24 <- subset_samples(ps.pruned, Time_TotalDays == "24")

## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.pruned.T24)

## Set order for Treatment
tse$Pretreatment <- factor(tse$Pretreatment, levels = c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"))


## ANCOMBC formula
output_T24 <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment", p_adj_method = "fdr", group = "Pretreatment",
                       alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 
View(output_T24$res)

## Compare treatments at T29 ====
ps.pruned.T29 <- subset_samples(ps.pruned, Time_TotalDays == "29")

## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.pruned.T29)

## Set order for Treatment
tse$Pretreatment <- factor(tse$Pretreatment, levels = c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"))


## ANCOMBC formula
output_T29 <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment", p_adj_method = "fdr", group = "Pretreatment",
                       alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 
View(output_T29$res)

## Compare treatments at T34 ====
ps.pruned.T34 <- subset_samples(ps.pruned, Time_TotalDays == "34")

## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.pruned.T34)

## Set order for Treatment
tse$Pretreatment <- factor(tse$Pretreatment, levels = c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"))


## ANCOMBC formula
output_T34 <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment", p_adj_method = "fdr", group = "Pretreatment",
                       alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 
View(output_T34$res)


## Attempting differential abundance -- globally ====
## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.pruned)

## Set order for Treatment
tse$Pretreatment <- factor(tse$Pretreatment, levels = c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"))

## ANCOMBC formula
output <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment + Time + Pretreatment*Time", p_adj_method = "fdr", group = "Pretreatment",
                       alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 

## Using this formula ideally takes both treatment, time and their interaction into consideration (along with tank as random effect) to allow us to make overall comparisons about treatment (while accounting for time)
## Get the warning boundary (singular) fit
View(output$res)

## Doesn't really give us anything; seems like time is the biggest driver of 


## try aldex for diff abundance 

## Attempting differential abundance by group ===
## All T0 vs active abx
ps.abx <- subset_samples(ps.pruned, Time_TotalDays == "0" | (Time_TotalDays == "2" & Pretreatment == "Antibiotics") |
                           (Time_TotalDays == "4" & Pretreatment == "Antibiotics"))

## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.abx)

## Set order for Treatment
tse$Pretreatment <- factor(tse$Pretreatment, levels = c("No Treatment", "Antibiotics"))


## ANCOMBC formula
output <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment + Time", rand_formula = "(1|TankID)", p_adj_method = "fdr", group = "Pretreatment",
                   alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 

View(output$res)

## blank T2-T4 vs antibiotic T2-T4
ps.abx <- subset_samples(ps.pruned, (Time_TotalDays == "2" & Pretreatment == "No Treatment") | 
                                  (Time_TotalDays == "4" & Pretreatment == "No Treatment") |
                                  (Time_TotalDays == "2" & Pretreatment == "Antibiotics") |
                                  (Time_TotalDays == "4" & Pretreatment == "Antibiotics"))

## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.abx)

## Set order for Treatment
tse$Pretreatment <- factor(tse$Pretreatment, levels = c("No Treatment", "Antibiotics"))


## ANCOMBC formula
output <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment + Time", rand_formula = "(1|TankID)", p_adj_method = "fdr", group = "Pretreatment",
                   alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 

View(output$res)


### Compare active antibiotics (T2-T4) against all T0 and blank T2-T4
ps.abx.active <- subset_samples(ps.pruned, Time_TotalDays == "0" | (Time_TotalDays == "2" & Pretreatment == "No Treatment") | 
                                  (Time_TotalDays == "4" & Pretreatment == "No Treatment") |
                                  (Time_TotalDays == "2" & Pretreatment == "Antibiotics") |
                                  (Time_TotalDays == "4" & Pretreatment == "Antibiotics"))

## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.abx.active)

## Set order for Treatment
tse$Pretreatment <- factor(tse$Pretreatment, levels = c("No Treatment", "Antibiotics"))


## ANCOMBC formula
output <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment + Time", rand_formula = "(1|TankID)", p_adj_method = "fdr", group = "Pretreatment",
                   alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 

View(output$res)

### Compare antibiotic recovery (T9-T34) against all T0 and blank T9-T34
ps.abx.recover <- subset_samples(ps.pruned, (Time_TotalDays != "0" & Time_TotalDays != "2" & Time_TotalDays != "4") & (Pretreatment == "No Treatment" | Pretreatment == "Antibiotics"))
  
## Convert phyloseq object to tree summarized experiment
tse <- mia::convertFromPhyloseq(ps.abx.recover)

## Set order for Treatment
tse$Pretreatment <- factor(tse$Pretreatment, levels = c("No Treatment", "Antibiotics"))


## ANCOMBC formula
output <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment + Time", rand_formula = "(1|TankID)", p_adj_method = "fdr", group = "Pretreatment",
                   alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 

View(output$res)
