# Author: Sunni Patton
# Date: 10/20/24 (last edited)
# Title: Beta Diversity

## Set seed ====
set.seed(123)

## Load libraries ====
library(microViz)
library(remotes)
library(pairwiseAdonis)
library(phyloseq)
library(ggpubr)
library(cowplot)
library(speedyseq)
library(dplyr)
library(ggplot2)

## Load data ====
ps.pruned <- readRDS(here::here("Output Files/04 - Phyloseq - Output/ps.pruned.RDS")) ## pruned data

## Edit data and subset ====
ps.pruned@sam_data$Phase <- paste0(ps.pruned@sam_data$Time_TotalDays)

ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "0"] <- "Antibiotic challenge"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "2"] <- "Antibiotic challenge"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "4"] <- "Antibiotic challenge"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "9"] <- "Temperature ramp"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "14"] <- "Temperature ramp"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "19"] <- "Temperature challenge"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "24"] <- "Temperature challenge"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "29"] <- "Temperature challenge"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "34"] <- "Temperature challenge"  


tax_fix(
  ps.pruned,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
) -> ps.pruned

# subset ps object by phase
subset_samples(ps.pruned, Phase == "Antibiotic challenge") -> ps.pruned.abx
subset_samples(ps.pruned, Phase == "Temperature ramp") -> ps.pruned.ramp
subset_samples(ps.pruned, Phase == "Temperature challenge") -> ps.pruned.temp

## CLR transform ====
ps.pruned.abx.clr <- ps.pruned.abx %>%
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA")

ps.pruned.ramp.clr <- ps.pruned.ramp %>%
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA")

ps.pruned.temp.clr <- ps.pruned.temp %>%
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA")

## Plot ====
ps.pruned.abx.clr@sam_data$Pretreatment <- factor(ps.pruned.abx.clr@sam_data$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
ps.pruned.abx.clr@sam_data$Phase <- factor(ps.pruned.abx.clr@sam_data$Phase, c("Antibiotic challenge", "Temperature ramp", "Temperature challenge"), ordered = TRUE)

ps.pruned.ramp.clr@sam_data$Pretreatment <- factor(ps.pruned.ramp.clr@sam_data$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
ps.pruned.ramp.clr@sam_data$Phase <- factor(ps.pruned.ramp.clr@sam_data$Phase, c("Antibiotic challenge", "Temperature ramp", "Temperature challenge"), ordered = TRUE)

ps.pruned.temp.clr@sam_data$Pretreatment <- factor(ps.pruned.temp.clr@sam_data$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
ps.pruned.temp.clr@sam_data$Phase <- factor(ps.pruned.temp.clr@sam_data$Phase, c("Antibiotic challenge", "Temperature ramp", "Temperature challenge"), ordered = TRUE)

## Antibiotic challenge
plot_abx <- ps.pruned.abx.clr %>%
  ord_plot( 
    colour = "Pretreatment",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + stat_ellipse(aes(group = Pretreatment, color = Pretreatment), linewidth = 1.5) + facet_wrap(~Phase) + theme_bw() + 
  scale_color_manual(values = c("black", "#672998", "#EF8737", "darkred")) + theme_bw() + 
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1)) +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold"),
        legend.position = "bottom") + 
  theme(strip.text = element_text(face = "bold", size = 12)) + labs(color = "Treatment")

## Temperature ramp
plot_ramp <- ps.pruned.ramp.clr %>%
  ord_plot( 
    colour = "Pretreatment",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + stat_ellipse(aes(group = Pretreatment, color = Pretreatment), linewidth = 1.5) + facet_wrap(~Phase) + theme_bw() + 
  scale_color_manual(values = c("black", "#672998", "#EF8737", "darkred")) + theme_bw() + 
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1)) +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold"),
        legend.position = "bottom") + 
  theme(strip.text = element_text(face = "bold", size = 12)) + labs(color = "Treatment")

## Temperature challenge
plot_temp <- ps.pruned.temp.clr %>%
  ord_plot( 
    colour = "Pretreatment",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + stat_ellipse(aes(group = Pretreatment, color = Pretreatment), linewidth = 1.5) + facet_wrap(~Phase) + theme_bw() + 
  scale_color_manual(values = c("black", "#672998", "#EF8737", "darkred")) + theme_bw() + 
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1)) +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold"),
        legend.position = "bottom") + 
  theme(strip.text = element_text(face = "bold", size = 12)) + labs(color = "Treatment")

x <- ggarrange(plot_abx, plot_ramp, plot_temp, ncol = 3, legend = "none")

ggplot2::ggsave(here::here("Output Files/07 - Beta Diversity - Output/plot_legend.png"), plot_temp,
                height = 300, width = 400, units = "mm",
                scale = 0.5, dpi = 1000)
ggplot2::ggsave(here::here("Output Files/07 - Beta Diversity - Output/plot_beta_phases_individ.png"), x,
                height = 250, width = 400, units = "mm",
                scale = 0.5, dpi = 1000)

## PERMANOVA and beta dispersion - Antibiotic challenge ====
# Calculate distance 
clr_euc_abx <- phyloseq::distance(ps.pruned.abx.clr, method = "euclidean")
# Extract sample data
sampledf_abx <- data.frame(sample_data(ps.pruned.abx.clr))
# adonis2
adonis_output_abx <- adonis2(clr_euc_abx ~ Pretreatment, data = sampledf_abx, perm = 9999, strata = sampledf_abx$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Antibiotic challenge", Comp = "No Treatment vs Antibiotics") # add columns
# dispersion stats
disp_abx <- betadisper(clr_euc_abx, sampledf_abx$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_abx) # significant 
stats_disp_abx <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_abx, sampledf_abx$Pretreatment, type = "centroid", 
                                                            bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
# dispersion dataframe
disp_df_abx <- data.frame(disp_abx$distances)
colnames(disp_df_abx)[colnames(disp_df_abx) == 'disp_abx.distances'] <- "Distances"
disp_df_abx <- cbind(disp_df_abx, sampledf_abx)
# plot
disp_df_abx$Pretreatment <- factor(disp_df_abx$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
plot_disp_abx <-  ggplot(disp_df_abx, aes(x = Pretreatment, y = Distances, color = Pretreatment, linewidth = 1.5))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + 
  scale_color_manual(values = c("black", "#672998","#EF8737","darkred"))

plot_disp_abx <- plot_disp_abx + geom_point(aes(color = Pretreatment), alpha = 0.5, 
                                            position = position_jitterdodge(jitter.width = 0.1)) +theme_bw()
plot_disp_abx <- plot_disp_abx + xlab("") + ylab("Distance to Centroid") #+ ylim(10, 35)
plot_disp_abx <- plot_disp_abx + theme(axis.text.x = element_blank(),
                                       axis.text.y = element_text(face = "bold", size = 12),
                                       axis.title = element_text(face = "bold", size = 12), 
                                       title = element_text(face = "bold"),
                                       panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
                                       strip.background = element_rect(
                                         color = "black", # Border color
                                         fill = "white", # Background fill color of the strip
                                         size = 1)) + 
  guides(color = guide_legend(title = "Treatment")) + theme(strip.text = element_text(face = "bold", size = 12)) +
  facet_wrap(~Phase)

## PERMANOVA and beta dispersion - Temperature ramp ====
# calculate distance
clr_euc_ramp <- phyloseq::distance(ps.pruned.ramp.clr, method = "euclidean")
# Extract sample data
sampledf_ramp <- data.frame(sample_data(ps.pruned.ramp.clr))

### No treatment vs antibiotics
selected_treatments <- c("No Treatment", "Antibiotics")
meta_pair <- sampledf_ramp[sampledf_ramp$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_ramp)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_ramp_NoTreatAbx <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature ramp", Comp = "No Treatment vs Antibiotics") # add columns
#### dispersion stats
disp_ramp_NoTreatAbx <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_ramp_NoTreatAbx) # not significant 
#### dispersion dataframe
disp_ramp_NoTreatAbx <- data.frame(disp_ramp_NoTreatAbx$distances)
colnames(disp_ramp_NoTreatAbx)[colnames(disp_ramp_NoTreatAbx) == 'disp_ramp_NoTreatAbx.distances'] <- "Distances"
disp_df_ramp_NoTreatAbx <- cbind(disp_ramp_NoTreatAbx, subset(sampledf_ramp, Pretreatment == "No Treatment" | Pretreatment == "Antibiotics"))

### No treatment vs temperature 
selected_treatments <- c("No Treatment", "Temperature")
meta_pair <- sampledf_ramp[sampledf_ramp$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_ramp)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_ramp_NoTreatTemp <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature ramp", Comp = "No Treatment vs Temperature")
#### dispersion stats
disp_ramp_NoTreatTemp <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_ramp_NoTreatTemp) # not sig
#### dispersion dataframe
disp_ramp_NoTreatTemp <- data.frame(disp_ramp_NoTreatTemp$distances)
colnames(disp_ramp_NoTreatTemp)[colnames(disp_ramp_NoTreatTemp) == 'disp_ramp_NoTreatTemp.distances'] <- "Distances"
disp_df_ramp_NoTreatTemp <- cbind(disp_ramp_NoTreatTemp, subset(sampledf_ramp, Pretreatment == "No Treatment" | Pretreatment == "Temperature"))

### No treatment vs abx + temp
selected_treatments <- c("No Treatment", "Antibiotics + Temperature")
meta_pair <- sampledf_ramp[sampledf_ramp$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_ramp)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_ramp_NoTreatAbxTemp <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature ramp", Comp = "No Treatment vs Antibiotics + Temperature")
#### dispersion stats
disp_ramp_NoTreatAbxTemp <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_ramp_NoTreatAbxTemp) 
#### dispersion dataframe
disp_ramp_NoTreatAbxTemp <- data.frame(disp_ramp_NoTreatAbxTemp$distances)
colnames(disp_ramp_NoTreatAbxTemp)[colnames(disp_ramp_NoTreatAbxTemp) == 'disp_ramp_NoTreatAbxTemp.distances'] <- "Distances"
disp_df_ramp_NoTreatAbxTemp <- cbind(disp_ramp_NoTreatAbxTemp, subset(sampledf_ramp, Pretreatment == "No Treatment" | Pretreatment == "Antibiotics + Temperature"))

### Antibiotics vs temp
selected_treatments <- c("Antibiotics", "Temperature")
meta_pair <- sampledf_ramp[sampledf_ramp$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_ramp)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_ramp_AbxTemp <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature ramp", Comp = "Antibiotics vs Temperature")
#### dispersion stats
disp_ramp_AbxTemp <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_ramp_AbxTemp)
#### dispersion dataframe
disp_ramp_AbxTemp <- data.frame(disp_ramp_AbxTemp$distances)
colnames(disp_ramp_AbxTemp)[colnames(disp_ramp_AbxTemp) == 'disp_ramp_AbxTemp.distances'] <- "Distances"
disp_df_ramp_AbxTemp <- cbind(disp_ramp_AbxTemp, subset(sampledf_ramp, Pretreatment == "Antibiotics" | Pretreatment == "Temperature"))

### Antibiotics vs temp + abx 
selected_treatments <- c("Antibiotics", "Antibiotics + Temperature")
meta_pair <- sampledf_ramp[sampledf_ramp$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_ramp)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_ramp_AbxAbxTemp <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature ramp", Comp = "Antibiotics vs Antibiotics + Temperature")
#### dispersion stats
disp_ramp_AbxAbxTemp <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_ramp_AbxAbxTemp) # not significant 
#### dispersion dataframe
disp_ramp_AbxAbxTemp <- data.frame(disp_ramp_AbxAbxTemp$distances)
colnames(disp_ramp_AbxAbxTemp)[colnames(disp_ramp_AbxAbxTemp) == 'disp_ramp_AbxAbxTemp.distances'] <- "Distances"
disp_df_ramp_AbxAbxTemp <- cbind(disp_ramp_AbxAbxTemp, subset(sampledf_ramp, Pretreatment == "Antibiotics" | Pretreatment == "Antibiotics + Temperature"))

### Antibiotics + Temp vs temp 
selected_treatments <- c("Temperature", "Antibiotics + Temperature")
meta_pair <- sampledf_ramp[sampledf_ramp$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_ramp)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_ramp_TempAbxTemp <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature ramp", Comp = "Temperature vs Antibiotics + Temperature")
#### dispersion stats
disp_ramp_TempAbxTemp <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_ramp_TempAbxTemp)
#### dispersion dataframe
disp_ramp_TempAbxTemp <- data.frame(disp_ramp_TempAbxTemp$distances)
colnames(disp_ramp_TempAbxTemp)[colnames(disp_ramp_TempAbxTemp) == 'disp_ramp_TempAbxTemp.distances'] <- "Distances"
disp_df_ramp_TempAbxTemp <- cbind(disp_ramp_TempAbxTemp, subset(sampledf_ramp, Pretreatment == "Temperature" | Pretreatment == "Antibiotics + Temperature"))

## combine files 
adonis_output_ramp <- rbind(adonis_output_abx, adonis_output_ramp_NoTreatAbx, adonis_output_ramp_NoTreatTemp,
                            adonis_output_ramp_NoTreatAbxTemp, adonis_output_ramp_AbxTemp, adonis_output_ramp_AbxAbxTemp,
                            adonis_output_ramp_TempAbxTemp)

disp_output_ramp <- rbind(disp_df_ramp_AbxTemp, disp_df_ramp_NoTreatAbx, disp_df_ramp_NoTreatTemp,
                          disp_df_ramp_NoTreatAbxTemp, disp_df_ramp_AbxAbxTemp, disp_df_ramp_TempAbxTemp)

# plot
disp_output_ramp$Pretreatment <- factor(disp_output_ramp$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
plot_disp_ramp <-  ggplot(disp_output_ramp, aes(x = Pretreatment, y = Distances, color = Pretreatment, linewidth = 1.5))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + 
  scale_color_manual(values = c("black", "#672998","#EF8737","darkred"))

plot_disp_ramp <- plot_disp_ramp + geom_point(aes(color = Pretreatment), alpha = 0.5, 
                                              position = position_jitterdodge(jitter.width = 0.1)) +theme_bw()
plot_disp_ramp <- plot_disp_ramp + xlab("") + ylab("Distance to Centroid") + ylim(10, 35)
plot_disp_ramp <- plot_disp_ramp + theme(axis.text.x = element_blank(),
                                         axis.text.y = element_text(face = "bold", size = 12),
                                         axis.title = element_text(face = "bold", size = 12), 
                                         panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
                                         strip.background = element_rect(
                                           color = "black", # Border color
                                           fill = "white", # Background fill color of the strip
                                           size = 1),
                                         title = element_text(face = "bold")) + 
  guides(color = guide_legend(title = "Treatment")) + theme(strip.text = element_text(face = "bold", size = 12)) +
  facet_wrap(~Phase)

## PERMANOVA and beta dispersion - Temperature challenge ====
# calculate distance
clr_euc_temp <- phyloseq::distance(ps.pruned.temp.clr, method = "euclidean")
# Extract sample data
sampledf_temp <- data.frame(sample_data(ps.pruned.temp.clr))

### No treatment vs antibiotics
selected_treatments <- c("No Treatment", "Antibiotics")
meta_pair <- sampledf_temp[sampledf_temp$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_temp)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_temp_NoTreatAbx <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature challenge", Comp = "No Treatment vs Antibiotics")
#### dispersion stats
disp_temp_NoTreatAbx <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_temp_NoTreatAbx)
#### dispersion dataframe
disp_temp_NoTreatAbx <- data.frame(disp_temp_NoTreatAbx$distances)
colnames(disp_temp_NoTreatAbx)[colnames(disp_temp_NoTreatAbx) == 'disp_temp_NoTreatAbx.distances'] <- "Distances"
disp_df_temp_NoTreatAbx <- cbind(disp_temp_NoTreatAbx, subset(sampledf_temp, Pretreatment == "No Treatment" | Pretreatment == "Antibiotics"))

### No treatment vs temperature 
selected_treatments <- c("No Treatment", "Temperature")
meta_pair <- sampledf_temp[sampledf_temp$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_temp)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_temp_NoTreatTemp <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature challenge", Comp = "No Treatment vs Temperature")
#### dispersion stats
disp_temp_NoTreatTemp <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_temp_NoTreatTemp)

stats_disp_NoTreatTemp <- broom::tidy(p.adjust(permutest(betadisper(dist_pair, meta_pair$Pretreatment, type = "centroid", 
                                                                    bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
#### dispersion dataframe 
disp_temp_NoTreatTemp <- data.frame(disp_temp_NoTreatTemp$distances)
colnames(disp_temp_NoTreatTemp)[colnames(disp_temp_NoTreatTemp) == 'disp_temp_NoTreatTemp.distances'] <- "Distances"
disp_df_temp_NoTreatTemp <- cbind(disp_temp_NoTreatTemp, subset(sampledf_temp, Pretreatment == "No Treatment" | Pretreatment == "Temperature"))

### No treatment vs abx + temp
selected_treatments <- c("No Treatment", "Antibiotics + Temperature")
meta_pair <- sampledf_temp[sampledf_temp$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_temp)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_temp_NoTreatAbxTemp <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature challenge", Comp = "No Treatment vs Antibiotics vs Temperature")
#### dispersion stats
disp_temp_NoTreatAbxTemp <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_temp_NoTreatAbxTemp)

stats_disp_NoTreatAbxTemp <- broom::tidy(p.adjust(permutest(betadisper(dist_pair, meta_pair$Pretreatment, type = "centroid", 
                                                                       bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))

#### dispersion dataframe
disp_temp_NoTreatAbxTemp <- data.frame(disp_temp_NoTreatAbxTemp$distances)
colnames(disp_temp_NoTreatAbxTemp)[colnames(disp_temp_NoTreatAbxTemp) == 'disp_temp_NoTreatAbxTemp.distances'] <- "Distances"
disp_df_temp_NoTreatAbxTemp <- cbind(disp_temp_NoTreatAbxTemp, subset(sampledf_temp, Pretreatment == "No Treatment" | Pretreatment == "Antibiotics + Temperature"))

### Antibiotics vs temp
selected_treatments <- c("Antibiotics", "Temperature")
meta_pair <- sampledf_temp[sampledf_temp$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_temp)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_temp_AbxTemp <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature challenge", Comp = "Antibiotics vs Temperature")
#### dispersion stats
disp_temp_AbxTemp <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_temp_AbxTemp)
#### dispersion dataframe 
disp_temp_AbxTemp <- data.frame(disp_temp_AbxTemp$distances)
colnames(disp_temp_AbxTemp)[colnames(disp_temp_AbxTemp) == 'disp_temp_AbxTemp.distances'] <- "Distances"
disp_df_temp_AbxTemp <- cbind(disp_temp_AbxTemp, subset(sampledf_temp, Pretreatment == "Temperature" | Pretreatment == "Antibiotics"))

### Antibiotics vs temp + abx 
selected_treatments <- c("Antibiotics", "Antibiotics + Temperature")
meta_pair <- sampledf_temp[sampledf_temp$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_temp)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_temp_AbxAbxTemp <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature challenge", Comp = "Antibiotics vs Antibiotics + Temperature")
#### dispersion stats
disp_temp_AbxAbxTemp <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_temp_AbxAbxTemp)
#### dispersion dataframe
disp_temp_AbxAbxTemp <- data.frame(disp_temp_AbxAbxTemp$distances)
colnames(disp_temp_AbxAbxTemp)[colnames(disp_temp_AbxAbxTemp) == 'disp_temp_AbxAbxTemp.distances'] <- "Distances"
disp_df_temp_AbxAbxTemp <- cbind(disp_temp_AbxAbxTemp, subset(sampledf_temp, Pretreatment == "Antibiotics" | Pretreatment == "Antibiotics + Temperature"))

### Antibiotics + Temp vs temp 
selected_treatments <- c("Temperature", "Antibiotics + Temperature")
meta_pair <- sampledf_temp[sampledf_temp$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_temp)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_temp_TempAbxTemp <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature challenge", Comp = "Antibiotics + Temperature vs Temperature")
#### dispersion stats
disp_temp_TempAbxTemp <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_temp_TempAbxTemp)
#### dispersion dataframe
disp_temp_TempAbxTemp <- data.frame(disp_temp_TempAbxTemp$distances)
colnames(disp_temp_TempAbxTemp)[colnames(disp_temp_TempAbxTemp) == 'disp_temp_TempAbxTemp.distances'] <- "Distances"
disp_df_temp_TempAbxTemp <- cbind(disp_temp_TempAbxTemp, subset(sampledf_temp, Pretreatment == "Temperature" | Pretreatment == "Antibiotics + Temperature"))

## combine
adonis_output_temp <- rbind(adonis_output_temp_NoTreatAbx, adonis_output_temp_NoTreatTemp,
                            adonis_output_temp_NoTreatAbxTemp, adonis_output_temp_AbxTemp, 
                            adonis_output_temp_AbxAbxTemp, adonis_output_temp_TempAbxTemp)

adonis_output_all <- rbind(adonis_output_ramp, adonis_output_temp)

disp_output_temp <- rbind(disp_df_temp_AbxTemp, disp_df_temp_NoTreatAbx, disp_df_temp_NoTreatTemp,
                          disp_df_temp_NoTreatAbxTemp, disp_df_temp_AbxAbxTemp, disp_df_temp_TempAbxTemp)
disp_output_all <- rbind(disp_output_ramp, disp_output_temp)

write.csv(adonis_output_all, here::here("Output Files/07 - Beta Diversity - Output/adonis_output_all.csv"))

# plot
## Temperature challenge
disp_output_temp$Pretreatment <- factor(disp_output_temp$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
plot_disp_temp <-  ggplot(disp_output_temp, aes(x = Pretreatment, y = Distances, color = Pretreatment, linewidth = 1.5))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + 
  scale_color_manual(values = c("black", "#672998","#EF8737","darkred"))

plot_disp_temp <- plot_disp_temp + geom_point(aes(color = Pretreatment), alpha = 0.5, 
                                              position = position_jitterdodge(jitter.width = 0.1)) +theme_bw()
plot_disp_temp <- plot_disp_temp + xlab("") + ylab("Distance to Centroid") + ylim(10, 35)
plot_disp_temp <- plot_disp_temp + theme(axis.text.x = element_blank(), 
                                         axis.text.y = element_text(face = "bold", size = 12),
                                         axis.title = element_text(face = "bold", size = 12), 
                                         title = element_text(face = "bold"),
                                         panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
                                         strip.background = element_rect(
                                           color = "black", # Border color
                                           fill = "white", # Background fill color of the strip
                                           size = 1)) + 
  guides(color = guide_legend(title = "Treatment")) + theme(strip.text = element_text(face = "bold", size = 12)) +
  facet_wrap(~Phase)

## Combine plots ====
x <- ggarrange(plot_disp_abx, plot_disp_ramp, plot_disp_temp, legend = "none", ncol = 3)

ggplot2::ggsave(here::here("Output Files/07 - Beta Diversity - Output/plot_legend_dispersion.png"), plot_disp_temp,
                height = 300, width = 400, units = "mm",
                scale = 0.5, dpi = 1000)
ggplot2::ggsave(here::here("Output Files/07 - Beta Diversity - Output/plot_disp_phases_individ.png"), x,
                height = 250, width = 400, units = "mm",
                scale = 0.5, dpi = 1000)



## Repeating analysis without A. rohweri ====
ps.pruned <- readRDS(here::here("Output Files/04 - Phyloseq - Output/ps.pruned.RDS")) ## pruned data

## Edit data and subset
ps.pruned@sam_data$Phase <- paste0(ps.pruned@sam_data$Time_TotalDays)

ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "0"] <- "Antibiotic challenge"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "2"] <- "Antibiotic challenge"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "4"] <- "Antibiotic challenge"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "9"] <- "Temperature ramp"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "14"] <- "Temperature ramp"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "19"] <- "Temperature challenge"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "24"] <- "Temperature challenge"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "29"] <- "Temperature challenge"  
ps.pruned@sam_data$Phase[ps.pruned@sam_data$Phase == "34"] <- "Temperature challenge"  


tax_fix(
  ps.pruned,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
) -> ps.pruned

### remove a. rohweri
subset_taxa(ps.pruned, ASV != "ASV1") -> ps.pruned.noASV1

### remove the samples that aren't in the ps.rare object
subset_samples(ps.pruned.noASV1, Samples != "G50_T15_1_c") -> ps.pruned.noASV1
subset_samples(ps.pruned.noASV1, Samples != "G50_T10_3_a")-> ps.pruned.noASV1
subset_samples(ps.pruned.noASV1, Samples != "G50_T10_3_d")-> ps.pruned.noASV1
subset_samples(ps.pruned.noASV1, Samples != "G50_T10_3_c")-> ps.pruned.noASV1
subset_samples(ps.pruned.noASV1, Samples != "G50_T15_1_d")-> ps.pruned.noASV1
subset_samples(ps.pruned.noASV1, Samples != "G50_T35_3_d")-> ps.pruned.noASV1
subset_samples(ps.pruned.noASV1, Samples != "G50_T15_3_c")-> ps.pruned.noASV1
subset_samples(ps.pruned.noASV1, Samples != "G50_T20_3_b")-> ps.pruned.noASV1
subset_samples(ps.pruned.noASV1, Samples != "G50_T20_1_b")-> ps.pruned.noASV1

### remove any taxa no longer present 
ps.pruned.noASV1 <- prune_taxa(taxa_sums(ps.pruned.noASV1@otu_table) > 0, ps.pruned.noASV1)

# subset ps object by phase
subset_samples(ps.pruned.noASV1, Phase == "Antibiotic challenge") -> ps.pruned.noASV1.abx
subset_samples(ps.pruned.noASV1, Phase == "Temperature ramp") -> ps.pruned.noASV1.ramp
subset_samples(ps.pruned.noASV1, Phase == "Temperature challenge") -> ps.pruned.noASV1.temp

# clr transform
ps.pruned.noASV1.abx.clr <- ps.pruned.noASV1.abx %>%
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA")

ps.pruned.noASV1.ramp.clr <- ps.pruned.noASV1.ramp %>%
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA")

ps.pruned.noASV1.temp.clr <- ps.pruned.noASV1.temp %>%
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA")

# plot
ps.pruned.noASV1.abx.clr@sam_data$Pretreatment <- factor(ps.pruned.noASV1.abx.clr@sam_data$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
ps.pruned.noASV1.abx.clr@sam_data$Phase <- factor(ps.pruned.noASV1.abx.clr@sam_data$Phase, c("Antibiotic challenge", "Temperature ramp", "Temperature challenge"), ordered = TRUE)

ps.pruned.noASV1.ramp.clr@sam_data$Pretreatment <- factor(ps.pruned.noASV1.ramp.clr@sam_data$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
ps.pruned.noASV1.ramp.clr@sam_data$Phase <- factor(ps.pruned.noASV1.ramp.clr@sam_data$Phase, c("Antibiotic challenge", "Temperature ramp", "Temperature challenge"), ordered = TRUE)

ps.pruned.noASV1.temp.clr@sam_data$Pretreatment <- factor(ps.pruned.noASV1.temp.clr@sam_data$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
ps.pruned.noASV1.temp.clr@sam_data$Phase <- factor(ps.pruned.noASV1.temp.clr@sam_data$Phase, c("Antibiotic challenge", "Temperature ramp", "Temperature challenge"), ordered = TRUE)

## Antibiotic challenge
plot_abx_noASV1 <- ps.pruned.noASV1.abx.clr %>%
  ord_plot( 
    colour = "Pretreatment",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + stat_ellipse(aes(group = Pretreatment, color = Pretreatment), linewidth = 1.5) + facet_wrap(~Phase) + theme_bw() + 
  scale_color_manual(values = c("black", "#672998", "#EF8737", "darkred")) + theme_bw() + 
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1)) +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold"),
        legend.position = "bottom") + 
  theme(strip.text = element_text(face = "bold", size = 12)) + labs(color = "Treatment")

## Temperature ramp
plot_ramp_noASV1 <- ps.pruned.noASV1.ramp.clr %>%
  ord_plot( 
    colour = "Pretreatment",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + stat_ellipse(aes(group = Pretreatment, color = Pretreatment), linewidth = 1.5) + facet_wrap(~Phase) + theme_bw() + 
  scale_color_manual(values = c("black", "#672998", "#EF8737", "darkred")) + theme_bw() + 
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1)) +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold"),
        legend.position = "bottom") + 
  theme(strip.text = element_text(face = "bold", size = 12)) + labs(color = "Treatment")

## Temperature challenge
plot_temp_noASV1 <- ps.pruned.noASV1.temp.clr %>%
  ord_plot( 
    colour = "Pretreatment",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + stat_ellipse(aes(group = Pretreatment, color = Pretreatment), linewidth = 1.5) + facet_wrap(~Phase) + theme_bw() + 
  scale_color_manual(values = c("black", "#672998", "#EF8737", "darkred")) + theme_bw() + 
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1)) +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold"),
        legend.position = "bottom") + 
  theme(strip.text = element_text(face = "bold", size = 12)) + labs(color = "Treatment")

x <- ggarrange(plot_abx_noASV1, plot_ramp_noASV1, plot_temp_noASV1, ncol = 3, legend = "none")

ggplot2::ggsave(here::here("Output Files/07 - Beta Diversity - Output/plot_beta_phases_noASV1.png"), x,
                height = 250, width = 400, units = "mm",
                scale = 0.5, dpi = 1000)

## PERMANOVA and beta dispersion - Antibiotic challenge
# Calculate distance 
clr_euc_abx_noASV1 <- phyloseq::distance(ps.pruned.noASV1.abx.clr, method = "euclidean")
# Extract sample data
sampledf_abx_noASV1 <- data.frame(sample_data(ps.pruned.noASV1.abx.clr))
# adonis2
adonis_output_abx_noASV1 <- adonis2(clr_euc_abx_noASV1 ~ Pretreatment, data = sampledf_abx_noASV1, perm = 9999, strata = sampledf_abx$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Antibiotic challenge", Comp = "No Treatment vs Antibiotics") # add columns
# dispersion stats
disp_abx_noASV1 <- betadisper(clr_euc_abx_noASV1, sampledf_abx_noASV1$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_abx_noASV1) # significant 
stats_disp_abx_noASV1 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_abx_noASV1, sampledf_abx_noASV1$Pretreatment, type = "centroid", 
                                                            bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
# dispersion dataframe
disp_df_abx_noASV1 <- data.frame(disp_abx_noASV1$distances)
colnames(disp_df_abx_noASV1)[colnames(disp_df_abx_noASV1) == 'disp_abx_noASV1.distances'] <- "Distances"
disp_df_abx_noASV1 <- cbind(disp_df_abx_noASV1, sampledf_abx_noASV1)
# plot
disp_df_abx_noASV1$Pretreatment <- factor(disp_df_abx_noASV1$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
plot_disp_abx_noASV1 <-  ggplot(disp_df_abx_noASV1, aes(x = Pretreatment, y = Distances, color = Pretreatment, linewidth = 1.5))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + 
  scale_color_manual(values = c("black", "#672998","#EF8737","darkred"))

plot_disp_abx_noASV1 <- plot_disp_abx_noASV1 + geom_point(aes(color = Pretreatment), alpha = 0.5, 
                                            position = position_jitterdodge(jitter.width = 0.1)) +theme_bw()
plot_disp_abx_noASV1 <- plot_disp_abx_noASV1 + xlab("") + ylab("Distance to Centroid") #+ ylim(10, 35)
plot_disp_abx_noASV1 <- plot_disp_abx_noASV1 + theme(axis.text.x = element_blank(),
                                       axis.text.y = element_text(face = "bold", size = 12),
                                       axis.title = element_text(face = "bold", size = 12), 
                                       title = element_text(face = "bold"),
                                       panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
                                       strip.background = element_rect(
                                         color = "black", # Border color
                                         fill = "white", # Background fill color of the strip
                                         size = 1)) + 
  guides(color = guide_legend(title = "Treatment")) + theme(strip.text = element_text(face = "bold", size = 12)) +
  facet_wrap(~Phase)



## PERMANOVA and beta dispersion - Temperature ramp
# calculate distance
clr_euc_ramp_noASV1 <- phyloseq::distance(ps.pruned.noASV1.ramp.clr, method = "euclidean")
# Extract sample data
sampledf_ramp_noASV1 <- data.frame(sample_data(ps.pruned.noASV1.ramp.clr))

### No treatment vs antibiotics
selected_treatments <- c("No Treatment", "Antibiotics")
meta_pair <- sampledf_ramp_noASV1[sampledf_ramp_noASV1$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_ramp_noASV1)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_ramp_NoTreatAbx_noASV1 <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature ramp", Comp = "No Treatment vs Antibiotics") # add columns
#### dispersion stats
disp_ramp_NoTreatAbx_noASV1 <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_ramp_NoTreatAbx_noASV1) # significant 
stats_disp_ramp_NoTreatAbx_noASV1 <- broom::tidy(p.adjust(permutest(betadisper(dist_pair, meta_pair$Pretreatment, type = "centroid", 
                                                                   bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))

#### dispersion dataframe
disp_ramp_NoTreatAbx_noASV1 <- data.frame(disp_ramp_NoTreatAbx_noASV1$distances)
colnames(disp_ramp_NoTreatAbx_noASV1)[colnames(disp_ramp_NoTreatAbx_noASV1) == 'disp_ramp_NoTreatAbx_noASV1.distances'] <- "Distances"
disp_df_ramp_NoTreatAbx_noASV1 <- cbind(disp_ramp_NoTreatAbx_noASV1, subset(sampledf_ramp_noASV1, Pretreatment == "No Treatment" | Pretreatment == "Antibiotics"))

### No treatment vs temperature
selected_treatments <- c("No Treatment", "Temperature")
meta_pair <- sampledf_ramp_noASV1[sampledf_ramp_noASV1$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_ramp_noASV1)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_ramp_NoTreatTemp_noASV1 <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature ramp", Comp = "No Treatment vs Temperature")
#### dispersion stats
disp_ramp_NoTreatTemp_noASV1 <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_ramp_NoTreatTemp_noASV1) # not sig
#### dispersion dataframe
disp_ramp_NoTreatTemp_noASV1 <- data.frame(disp_ramp_NoTreatTemp_noASV1$distances)
colnames(disp_ramp_NoTreatTemp_noASV1)[colnames(disp_ramp_NoTreatTemp_noASV1) == 'disp_ramp_NoTreatTemp_noASV1.distances'] <- "Distances"
disp_df_ramp_NoTreatTemp_noASV1 <- cbind(disp_ramp_NoTreatTemp_noASV1, subset(sampledf_ramp_noASV1, Pretreatment == "No Treatment" | Pretreatment == "Temperature"))

### No treatment vs abx + temp
selected_treatments <- c("No Treatment", "Antibiotics + Temperature")
meta_pair <- sampledf_ramp_noASV1[sampledf_ramp_noASV1$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_ramp_noASV1)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_ramp_NoTreatAbxTemp_noASV1 <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature ramp", Comp = "No Treatment vs Antibiotics + Temperature")
#### dispersion stats
disp_ramp_NoTreatAbxTemp_noASV1 <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_ramp_NoTreatAbxTemp_noASV1) # sig
stats_disp_ramp_NoTreatAbxTemp_noASV1 <- broom::tidy(p.adjust(permutest(betadisper(dist_pair, meta_pair$Pretreatment, type = "centroid", 
                                                                               bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))

#### dispersion dataframe
disp_ramp_NoTreatAbxTemp_noASV1 <- data.frame(disp_ramp_NoTreatAbxTemp_noASV1$distances)
colnames(disp_ramp_NoTreatAbxTemp_noASV1)[colnames(disp_ramp_NoTreatAbxTemp_noASV1) == 'disp_ramp_NoTreatAbxTemp_noASV1.distances'] <- "Distances"
disp_df_ramp_NoTreatAbxTemp_noASV1 <- cbind(disp_ramp_NoTreatAbxTemp_noASV1, subset(sampledf_ramp_noASV1, Pretreatment == "No Treatment" | Pretreatment == "Antibiotics + Temperature"))

### Antibiotics vs temp
selected_treatments <- c("Antibiotics", "Temperature")
meta_pair <- sampledf_ramp_noASV1[sampledf_ramp_noASV1$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_ramp_noASV1)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_ramp_AbxTemp_noASV1 <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature ramp", Comp = "Antibiotics vs Temperature")
#### dispersion stats
disp_ramp_AbxTemp_noASV1 <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_ramp_AbxTemp_noASV1) # sig
stats_disp_ramp_AbxTemp_noASV1 <- broom::tidy(p.adjust(permutest(betadisper(dist_pair, meta_pair$Pretreatment, type = "centroid", 
                                                                                   bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))

#### dispersion dataframe
disp_ramp_AbxTemp_noASV1 <- data.frame(disp_ramp_AbxTemp_noASV1$distances)
colnames(disp_ramp_AbxTemp_noASV1)[colnames(disp_ramp_AbxTemp_noASV1) == 'disp_ramp_AbxTemp_noASV1.distances'] <- "Distances"
disp_df_ramp_AbxTemp_noASV1 <- cbind(disp_ramp_AbxTemp_noASV1, subset(sampledf_ramp_noASV1, Pretreatment == "Antibiotics" | Pretreatment == "Temperature"))

### Antibiotics vs temp + abx
selected_treatments <- c("Antibiotics", "Antibiotics + Temperature")
meta_pair <- sampledf_ramp_noASV1[sampledf_ramp_noASV1$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_ramp_noASV1)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_ramp_AbxAbxTemp_noASV1 <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature ramp", Comp = "Antibiotics vs Antibiotics + Temperature")
#### dispersion stats
disp_ramp_AbxAbxTemp_noASV1 <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_ramp_AbxAbxTemp_noASV1) # not significant 
#### dispersion dataframe
disp_ramp_AbxAbxTemp_noASV1 <- data.frame(disp_ramp_AbxAbxTemp_noASV1$distances)
colnames(disp_ramp_AbxAbxTemp_noASV1)[colnames(disp_ramp_AbxAbxTemp_noASV1) == 'disp_ramp_AbxAbxTemp_noASV1.distances'] <- "Distances"
disp_df_ramp_AbxAbxTemp_noASV1 <- cbind(disp_ramp_AbxAbxTemp_noASV1, subset(sampledf_ramp_noASV1, Pretreatment == "Antibiotics" | Pretreatment == "Antibiotics + Temperature"))

### Antibiotics + Temp vs temp
selected_treatments <- c("Temperature", "Antibiotics + Temperature")
meta_pair <- sampledf_ramp_noASV1[sampledf_ramp_noASV1$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_ramp_noASV1)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_ramp_TempAbxTemp_noASV1 <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature ramp", Comp = "Temperature vs Antibiotics + Temperature")
#### dispersion stats
disp_ramp_TempAbxTemp_noASV1 <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_ramp_TempAbxTemp_noASV1) # sig
stats_disp_ramp_TempAbxTemp_noASV1 <- broom::tidy(p.adjust(permutest(betadisper(dist_pair, meta_pair$Pretreatment, type = "centroid", 
                                                                            bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))

#### dispersion dataframe
disp_ramp_TempAbxTemp_noASV1 <- data.frame(disp_ramp_TempAbxTemp_noASV1$distances)
colnames(disp_ramp_TempAbxTemp_noASV1)[colnames(disp_ramp_TempAbxTemp_noASV1) == 'disp_ramp_TempAbxTemp_noASV1.distances'] <- "Distances"
disp_df_ramp_TempAbxTemp_noASV1 <- cbind(disp_ramp_TempAbxTemp_noASV1, subset(sampledf_ramp_noASV1, Pretreatment == "Temperature" | Pretreatment == "Antibiotics + Temperature"))

## combine files 
adonis_output_ramp_noASV1 <- rbind(adonis_output_abx_noASV1, adonis_output_ramp_NoTreatAbx_noASV1, adonis_output_ramp_NoTreatTemp_noASV1,
                                   adonis_output_ramp_NoTreatAbxTemp_noASV1, adonis_output_ramp_AbxTemp_noASV1, adonis_output_ramp_AbxAbxTemp_noASV1,
                                   adonis_output_ramp_TempAbxTemp_noASV1)

disp_output_ramp_noASV1 <- rbind(disp_df_ramp_AbxTemp_noASV1, disp_df_ramp_NoTreatAbx_noASV1, disp_df_ramp_NoTreatTemp_noASV1,
                                 disp_df_ramp_NoTreatAbxTemp_noASV1, disp_df_ramp_AbxAbxTemp_noASV1, disp_df_ramp_TempAbxTemp_noASV1)

# plot
disp_output_ramp_noASV1$Pretreatment <- factor(disp_output_ramp_noASV1$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
plot_disp_ramp_noASV1 <-  ggplot(disp_output_ramp_noASV1, aes(x = Pretreatment, y = Distances, color = Pretreatment, linewidth = 1.5))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + 
  scale_color_manual(values = c("black", "#672998","#EF8737","darkred"))

plot_disp_ramp_noASV1 <- plot_disp_ramp_noASV1 + geom_point(aes(color = Pretreatment), alpha = 0.5, 
                                                            position = position_jitterdodge(jitter.width = 0.1)) +theme_bw()
plot_disp_ramp_noASV1 <- plot_disp_ramp_noASV1 + xlab("") + ylab("Distance to Centroid") + ylim(10, 35)
plot_disp_ramp_noASV1 <- plot_disp_ramp_noASV1 + theme(axis.text.x = element_blank(),
                                                       axis.text.y = element_text(face = "bold", size = 12),
                                                       axis.title = element_text(face = "bold", size = 12), 
                                                       panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
                                                       strip.background = element_rect(
                                                         color = "black", # Border color
                                                         fill = "white", # Background fill color of the strip
                                                         size = 1),
                                                       title = element_text(face = "bold")) + 
  guides(color = guide_legend(title = "Treatment")) + theme(strip.text = element_text(face = "bold", size = 12)) +
  facet_wrap(~Phase)
## PERMANOVA and beta dispersion - Temperature challenge
# calculate distance
clr_euc_temp_noASV1 <- phyloseq::distance(ps.pruned.noASV1.temp.clr, method = "euclidean")
# Extract sample data
sampledf_temp_noASV1 <- data.frame(sample_data(ps.pruned.noASV1.temp.clr))

### No treatment vs antibiotics
selected_treatments <- c("No Treatment", "Antibiotics")
meta_pair <- sampledf_temp_noASV1[sampledf_temp_noASV1$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_temp_noASV1)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_temp_NoTreatAbx_noASV1 <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature challenge", Comp = "No Treatment vs Antibiotics")
#### dispersion stats
disp_temp_NoTreatAbx_noASV1 <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_temp_NoTreatAbx_noASV1) # not sig
#### dispersion dataframe
disp_temp_NoTreatAbx_noASV1 <- data.frame(disp_temp_NoTreatAbx_noASV1$distances)
colnames(disp_temp_NoTreatAbx_noASV1)[colnames(disp_temp_NoTreatAbx_noASV1) == 'disp_temp_NoTreatAbx_noASV1.distances'] <- "Distances"
disp_df_temp_NoTreatAbx_noASV1 <- cbind(disp_temp_NoTreatAbx_noASV1, subset(sampledf_temp_noASV1, Pretreatment == "No Treatment" | Pretreatment == "Antibiotics"))

### No treatment vs temperature
selected_treatments <- c("No Treatment", "Temperature")
meta_pair <- sampledf_temp_noASV1[sampledf_temp_noASV1$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_temp_noASV1)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_temp_NoTreatTemp_noASV1 <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature challenge", Comp = "No Treatment vs Temperature")
#### dispersion stats
disp_temp_NoTreatTemp_noASV1 <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_temp_NoTreatTemp_noASV1)
#### dispersion dataframe 
disp_temp_NoTreatTemp_noASV1 <- data.frame(disp_temp_NoTreatTemp_noASV1$distances)
colnames(disp_temp_NoTreatTemp_noASV1)[colnames(disp_temp_NoTreatTemp_noASV1) == 'disp_temp_NoTreatTemp_noASV1.distances'] <- "Distances"
disp_df_temp_NoTreatTemp_noASV1 <- cbind(disp_temp_NoTreatTemp_noASV1, subset(sampledf_temp_noASV1, Pretreatment == "No Treatment" | Pretreatment == "Temperature"))

### No treatment vs abx + temp
selected_treatments <- c("No Treatment", "Antibiotics + Temperature")
meta_pair <- sampledf_temp_noASV1[sampledf_temp_noASV1$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_temp_noASV1)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_temp_NoTreatAbxTemp_noASV1 <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature challenge", Comp = "No Treatment vs Antibiotics vs Temperature")
#### dispersion stats
disp_temp_NoTreatAbxTemp_noASV1 <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_temp_NoTreatAbxTemp_noASV1)
#### dispersion dataframe
disp_temp_NoTreatAbxTemp_noASV1 <- data.frame(disp_temp_NoTreatAbxTemp_noASV1$distances)
colnames(disp_temp_NoTreatAbxTemp_noASV1)[colnames(disp_temp_NoTreatAbxTemp_noASV1) == 'disp_temp_NoTreatAbxTemp_noASV1.distances'] <- "Distances"
disp_df_temp_NoTreatAbxTemp_noASV1 <- cbind(disp_temp_NoTreatAbxTemp_noASV1, subset(sampledf_temp_noASV1, Pretreatment == "No Treatment" | Pretreatment == "Antibiotics + Temperature"))

### Antibiotics vs temp
selected_treatments <- c("Antibiotics", "Temperature")
meta_pair <- sampledf_temp_noASV1[sampledf_temp_noASV1$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_temp_noASV1)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_temp_AbxTemp_noASV1 <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature challenge", Comp = "Antibiotics vs Temperature")
#### dispersion stats
disp_temp_AbxTemp_noASV1 <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_temp_AbxTemp_noASV1)
stats_disp_temp_AbxTemp_noASV1 <- broom::tidy(p.adjust(permutest(betadisper(dist_pair, meta_pair$Pretreatment, type = "centroid", 
                                                                            bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
#### dispersion dataframe 
disp_temp_AbxTemp_noASV1 <- data.frame(disp_temp_AbxTemp_noASV1$distances)
colnames(disp_temp_AbxTemp_noASV1)[colnames(disp_temp_AbxTemp_noASV1) == 'disp_temp_AbxTemp_noASV1.distances'] <- "Distances"
disp_df_temp_AbxTemp_noASV1 <- cbind(disp_temp_AbxTemp_noASV1, subset(sampledf_temp_noASV1, Pretreatment == "Temperature" | Pretreatment == "Antibiotics"))


### Antibiotics vs temp + abx
selected_treatments <- c("Antibiotics", "Antibiotics + Temperature")
meta_pair <- sampledf_temp_noASV1[sampledf_temp_noASV1$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_temp_noASV1)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_temp_AbxAbxTemp_noASV1 <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature challenge", Comp = "Antibiotics vs Antibiotics + Temperature")
#### dispersion stats
disp_temp_AbxAbxTemp_noASV1 <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_temp_AbxAbxTemp_noASV1)
#### dispersion dataframe
disp_temp_AbxAbxTemp_noASV1 <- data.frame(disp_temp_AbxAbxTemp_noASV1$distances)
colnames(disp_temp_AbxAbxTemp_noASV1)[colnames(disp_temp_AbxAbxTemp_noASV1) == 'disp_temp_AbxAbxTemp_noASV1.distances'] <- "Distances"
disp_df_temp_AbxAbxTemp_noASV1 <- cbind(disp_temp_AbxAbxTemp_noASV1, subset(sampledf_temp_noASV1, Pretreatment == "Antibiotics" | Pretreatment == "Antibiotics + Temperature"))

### Antibiotics + Temp vs temp
selected_treatments <- c("Temperature", "Antibiotics + Temperature")
meta_pair <- sampledf_temp_noASV1[sampledf_temp_noASV1$Pretreatment %in% selected_treatments, ]
dist_pair <- as.matrix(clr_euc_temp_noASV1)[meta_pair$Samples, meta_pair$Samples]
dist_pair <- as.dist(dist_pair)
#### adonis2
adonis_output_temp_TempAbxTemp_noASV1 <- adonis2(dist_pair ~ Pretreatment, data = meta_pair, perm = 9999, strata = meta_pair$Time_TotalDays) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "adj.p", method = "fdr") %>% # adjust p value
  broom::tidy() %>%
  dplyr::mutate(Phase = "Temperature challenge", Comp = "Antibiotics + Temperature vs Temperature")
#### dispersion stats
disp_temp_TempAbxTemp_noASV1 <- betadisper(dist_pair, meta_pair$Pretreatment, bias.adjust = TRUE, type = "centroid")
anova(disp_temp_TempAbxTemp_noASV1)
stats_disp_temp_TempAbxTemp_noASV1 <- broom::tidy(p.adjust(permutest(betadisper(dist_pair, meta_pair$Pretreatment, type = "centroid", 
                                                                                bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
#### dispersion dataframe
disp_temp_TempAbxTemp_noASV1 <- data.frame(disp_temp_TempAbxTemp_noASV1$distances)
colnames(disp_temp_TempAbxTemp_noASV1)[colnames(disp_temp_TempAbxTemp_noASV1) == 'disp_temp_TempAbxTemp_noASV1.distances'] <- "Distances"
disp_df_temp_TempAbxTemp_noASV1 <- cbind(disp_temp_TempAbxTemp_noASV1, subset(sampledf_temp_noASV1, Pretreatment == "Temperature" | Pretreatment == "Antibiotics + Temperature"))

## combine
adonis_output_temp_noASV1 <- rbind(adonis_output_temp_NoTreatAbx_noASV1, adonis_output_temp_NoTreatTemp_noASV1,
                                   adonis_output_temp_NoTreatAbxTemp_noASV1, adonis_output_temp_AbxTemp_noASV1, 
                                   adonis_output_temp_AbxAbxTemp_noASV1, adonis_output_temp_TempAbxTemp_noASV1)

adonis_output_all <- rbind(adonis_output_ramp_noASV1, adonis_output_temp_noASV1)

disp_output_temp_noASV1 <- rbind(disp_df_temp_AbxTemp_noASV1, disp_df_temp_NoTreatAbx_noASV1, disp_df_temp_NoTreatTemp_noASV1,
                                 disp_df_temp_NoTreatAbxTemp_noASV1, disp_df_temp_AbxAbxTemp_noASV1, disp_df_temp_TempAbxTemp_noASV1)
disp_output_all_noASV1 <- rbind(disp_output_ramp_noASV1, disp_output_temp_noASV1)

write.csv(adonis_output_all, here::here("Output Files/07 - Beta Diversity - Output/adonis_output_all_noASV1.csv"))

# plot
## Temperature challenge
disp_output_temp_noASV1$Pretreatment <- factor(disp_output_temp_noASV1$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
plot_disp_temp_noASV1 <-  ggplot(disp_output_temp, aes(x = Pretreatment, y = Distances, color = Pretreatment, linewidth = 1.5))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + 
  scale_color_manual(values = c("black", "#672998","#EF8737","darkred"))

plot_disp_temp_noASV1 <- plot_disp_temp_noASV1 + geom_point(aes(color = Pretreatment), alpha = 0.5, 
                                                            position = position_jitterdodge(jitter.width = 0.1)) +theme_bw()
plot_disp_temp_noASV1 <- plot_disp_temp_noASV1 + xlab("") + ylab("Distance to Centroid") + ylim(10, 35)
plot_disp_temp_noASV1 <- plot_disp_temp_noASV1 + theme(axis.text.x = element_blank(), 
                                                       axis.text.y = element_text(face = "bold", size = 12),
                                                       axis.title = element_text(face = "bold", size = 12), 
                                                       title = element_text(face = "bold"),
                                                       panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
                                                       strip.background = element_rect(
                                                         color = "black", # Border color
                                                         fill = "white", # Background fill color of the strip
                                                         size = 1)) + 
  guides(color = guide_legend(title = "Treatment")) + theme(strip.text = element_text(face = "bold", size = 12)) +
  facet_wrap(~Phase)


## Combine plots
x <- ggarrange(plot_disp_abx_noASV1, plot_disp_ramp_noASV1, plot_disp_temp_noASV1, legend = "none", ncol = 3)


ggplot2::ggsave(here::here("Output Files/07 - Beta Diversity - Output/plot_disp_phases_individ_noASV1.png"), x,
                height = 250, width = 400, units = "mm",
                scale = 0.5, dpi = 1000)
