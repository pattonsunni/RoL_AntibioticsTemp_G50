# Author: Sunni Patton
# Last edited: 10/18/24 
# Title: Relative Abundance

## Set seed ====
set.seed(123)

## Load libraries ====
library(microViz)
library(dplyr)
library(phyloseq)
library(ggpubr)
library(svglite)
library(tidyverse)

## Load data ====
readRDS(here::here("Output Files/04 - Phyloseq - Output/ps.rare.rds")) -> ps.rare

## Validate phyloseq object ====
tax_fix(
  ps.rare,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
) -> ps.rare

ps.rare <- phyloseq_validate(ps.rare) 

## Change taxa names ====
# Agree with current literature
taxa_change <- data.frame(tax_table(ps.rare))

for (i in 1:nrow(taxa_change)){
  
  if (taxa_change[i,6] == "MD3-55"){
    taxa_change[i, 6:7] <- "Aquarickettsia"
  }
}

tax_table(ps.rare) <- as.matrix(taxa_change)

saveRDS(ps.rare, here::here("Output Files/05 - Relative Abundance - Output/ps.rare.rds"))

## Relative abundance ps object ====
ps.rare.trans <- transform_sample_counts(ps.rare, function(OTU) OTU/sum(OTU))

## Identify top 15 most abundant taxa
top15 <- names(sort(taxa_sums(ps.rare.trans), decreasing = TRUE))[1:15] 

## Prune top 15
ps.rare.top15 <- prune_taxa(top15, ps.rare.trans)

## Plot relative abundance (Supplemental) ====
as.data.frame(ps.rare.top15@tax_table) -> tax
tax$GenusASV <- paste0(tax$Genus, " (", tax$ASV, ")")
as.matrix(tax) -> tax
tax_table(ps.rare.top15) <- tax

colors <- c("#66A61EFF", "#f38400", "#7570B3FF", "#2D2651FF", "#C969A1FF", "#be0032",
            "#dcd300", "#654522","#e25822", "#5D7298FF",
            "#a1caf1", "#81B28DFF", "#A6761DFF", "#f38400", "#c2b280")

plot.relAbund.blank <- plot_bar(subset_samples(ps.rare.top15, ps.rare.top15@sam_data$Pretreatment == "No Treatment"), x="Sample", fill="GenusASV")

plot.relAbund.blank <- plot.relAbund.blank + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), 
        title = element_text(face = "bold"), 
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1.5),
        legend.position = "none") + 
  xlab("Samples") + facet_wrap(~Pretreatment) + scale_fill_manual(values = colors)

plot.relAbund.abx <- plot_bar(subset_samples(ps.rare.top15, ps.rare.top15@sam_data$Pretreatment == "Antibiotics"), x="Sample", fill="GenusASV")

plot.relAbund.abx <- plot.relAbund.abx + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), 
        title = element_text(face = "bold"), 
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1.5),
        legend.position = "none") + 
  xlab("Samples") + facet_wrap(~Pretreatment) + scale_fill_manual(values = colors)


plot.relAbund.temp <- plot_bar(subset_samples(ps.rare.top15, ps.rare.top15@sam_data$Pretreatment == "Temperature"), x="Sample", fill="GenusASV")

plot.relAbund.temp <- plot.relAbund.temp + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), 
        title = element_text(face = "bold"), 
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1.5),
        legend.position = "none") + 
  xlab("Samples") + facet_wrap(~Pretreatment) + scale_fill_manual(values = colors)



plot.relAbund.abxtemp <- plot_bar(subset_samples(ps.rare.top15, ps.rare.top15@sam_data$Pretreatment == "Antibiotics + Temperature"), x="Sample", fill="GenusASV")

plot.relAbund.abxtemp <- plot.relAbund.abxtemp + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), 
        title = element_text(face = "bold"), 
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1.5),
        legend.position = "none") + 
  xlab("Samples") + facet_wrap(~Pretreatment) + scale_fill_manual(values = colors)


plot_relabund_all <- ggarrange(plot.relAbund.blank, plot.relAbund.abx, plot.relAbund.temp, plot.relAbund.abxtemp)

ggplot2::ggsave(here::here("Output Files/05 - Relative Abundance - Output/plot_relabund_all_legend.svg"), plot.relAbund.blank,
                height = 500, width = 600, units = "mm",
                scale = 0.5, dpi = 1000)

## Plot mean relative abundance ====
colors <- c("#66A61EFF", "#f38400", "#7570B3FF", "#2D2651FF", "#C969A1FF", "#be0032",
            "#dcd300", "#654522","#e25822", "#5D7298FF",
            "#a1caf1", "#81B28DFF", "#A6761DFF", "#f38400", "#c2b280")
## Blank
melt.test <- psmelt(ps.rare.top15)

write.csv(melt.test, here::here("Output Files/05 - Relative Abundance - Output/relabund_all.csv"))

melt.test %>%
  filter(Time == "0" & Pretreatment == "No Treatment") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Blank.T0
melt.Blank.T0$Time <- "0"
melt.Blank.T0$Treatment <- "No Treatment"

melt.test %>%
  filter(Time == "2" & Pretreatment == "No Treatment") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Blank.T2
melt.Blank.T2$Time <- "2"
melt.Blank.T2$Treatment <- "No Treatment"

melt.test %>%
  filter(Time == "5" & Pretreatment == "No Treatment") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Blank.T4
melt.Blank.T4$Time <- "4"
melt.Blank.T4$Treatment <- "No Treatment"

melt.test %>%
  filter(Time == "10" & Pretreatment == "No Treatment") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Blank.T9
melt.Blank.T9$Time <- "9"
melt.Blank.T9$Treatment <- "No Treatment"

melt.test %>%
  filter(Time == "15" & Pretreatment == "No Treatment") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Blank.T14
melt.Blank.T14$Time <- "14"
melt.Blank.T14$Treatment <- "No Treatment"

melt.test %>%
  filter(Time == "20" & Pretreatment == "No Treatment") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Blank.T19
melt.Blank.T19$Time <- "19"
melt.Blank.T19$Treatment <- "No Treatment"

melt.test %>%
  filter(Time == "25" & Pretreatment == "No Treatment") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Blank.T24
melt.Blank.T24$Time <- "24"
melt.Blank.T24$Treatment <- "No Treatment"

melt.test %>%
  filter(Time == "30" & Pretreatment == "No Treatment") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Blank.T29
melt.Blank.T29$Time <- "29"
melt.Blank.T29$Treatment <- "No Treatment"

melt.test %>%
  filter(Time == "35" & Pretreatment == "No Treatment") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Blank.T34
melt.Blank.T34$Time <- "34"
melt.Blank.T34$Treatment <- "No Treatment"

rbind(melt.Blank.T0, melt.Blank.T2, melt.Blank.T4, melt.Blank.T9, melt.Blank.T14, melt.Blank.T19, melt.Blank.T24,
      melt.Blank.T29, melt.Blank.T34) -> blank.relabund

blank.relabund$mean <- (blank.relabund$mean)*100
blank.relabund$sd <- (blank.relabund$sd)*100

as.data.frame(ps.rare.top15@tax_table) -> GenusASV
data.frame(GenusASV$ASV, GenusASV$Genus) -> GenusASV
colnames(GenusASV) <- c("ASV", "Genus")
           
merge(blank.relabund, GenusASV, by = "ASV") -> blank.relabund
blank.relabund$GenusASV <- paste0(blank.relabund$Genus, " (", blank.relabund$ASV, ")")

## Plot
blank.relabund$Time <- factor(blank.relabund$Time, c("0", "2", "4", "9", "14", "19","24",
                                                 "29", "34"), ordered = TRUE)

blank_MRA <- ggplot(data = blank.relabund, aes(x = Time, y = mean, fill = GenusASV)) + 
  geom_bar(stat = "identity") + facet_wrap(~Treatment) 

blank_MRA <- blank_MRA + 
  xlab("Time (Days)") + ylab("Mean Relative Abundance (%)") + 
  guides(fill = guide_legend(title = "ASV")) + 
  theme_bw() + 
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1.5)) +
  scale_fill_manual(values = colors)
## Antibiotics
melt.test %>%
  filter(Time == "2" & Pretreatment == "Antibiotics") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Abx.T2
melt.Abx.T2$Time <- "2"
melt.Abx.T2$Treatment <- "Antibiotics"

melt.test %>%
  filter(Time == "5" & Pretreatment == "Antibiotics") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Abx.T4
melt.Abx.T4$Time <- "4"
melt.Abx.T4$Treatment <- "Antibiotics"

melt.test %>%
  filter(Time == "10" & Pretreatment == "Antibiotics") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Abx.T9
melt.Abx.T9$Time <- "9"
melt.Abx.T9$Treatment <- "Antibiotics"

melt.test %>%
  filter(Time == "15" & Pretreatment == "Antibiotics") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Abx.T14
melt.Abx.T14$Time <- "14"
melt.Abx.T14$Treatment <- "Antibiotics"

melt.test %>%
  filter(Time == "20" & Pretreatment == "Antibiotics") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Abx.T19
melt.Abx.T19$Time <- "19"
melt.Abx.T19$Treatment <- "Antibiotics"

melt.test %>%
  filter(Time == "25" & Pretreatment == "Antibiotics") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Abx.T24
melt.Abx.T24$Time <- "24"
melt.Abx.T24$Treatment <- "Antibiotics"

melt.test %>%
  filter(Time == "30" & Pretreatment == "Antibiotics") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Abx.T29
melt.Abx.T29$Time <- "29"
melt.Abx.T29$Treatment <- "Antibiotics"

melt.test %>%
  filter(Time == "35" & Pretreatment == "Antibiotics") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Abx.T34
melt.Abx.T34$Time <- "34"
melt.Abx.T34$Treatment <- "Antibiotics"

rbind(melt.Abx.T2, melt.Abx.T4, melt.Abx.T9, melt.Abx.T14, melt.Abx.T19, melt.Abx.T24,
      melt.Abx.T29, melt.Abx.T34) -> abx.relabund

abx.relabund$mean <- (abx.relabund$mean)*100
abx.relabund$sd <- (abx.relabund$sd)*100

as.data.frame(ps.rare.top15@tax_table) -> GenusASV
data.frame(GenusASV$ASV, GenusASV$Genus) -> GenusASV
colnames(GenusASV) <- c("ASV", "Genus")

merge(abx.relabund, GenusASV, by = "ASV") -> abx.relabund
abx.relabund$GenusASV <- paste0(abx.relabund$Genus, " (", abx.relabund$ASV, ")")

## Plot
abx.relabund$Time <- factor(abx.relabund$Time, c("2", "4", "9", "14", "19","24",
                                                 "29", "34"), ordered = TRUE)

abx_MRA <- ggplot(data = abx.relabund, aes(x = Time, y = mean, fill = GenusASV)) + 
  geom_bar(stat = "identity") + facet_wrap(~Treatment)

abx_MRA <- abx_MRA + 
  xlab("Time (Days)") + ylab("Mean Relative Abundance (%)") + guides(fill = guide_legend(title = "ASV")) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1.5)) +
  scale_fill_manual(values = colors)


## Temperature 
melt.test %>%
  filter(Time == "10" & Pretreatment == "Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Temp.T9
melt.Temp.T9$Time <- "9"
melt.Temp.T9$Treatment <- "Temperature"

melt.test %>%
  filter(Time == "15" & Pretreatment == "Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Temp.T14
melt.Temp.T14$Time <- "14"
melt.Temp.T14$Treatment <- "Temperature"

melt.test %>%
  filter(Time == "20" & Pretreatment == "Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Temp.T19
melt.Temp.T19$Time <- "19"
melt.Temp.T19$Treatment <- "Temperature"

melt.test %>%
  filter(Time == "25" & Pretreatment == "Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Temp.T24
melt.Temp.T24$Time <- "24"
melt.Temp.T24$Treatment <- "Temperature"

melt.test %>%
  filter(Time == "30" & Pretreatment == "Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Temp.T29
melt.Temp.T29$Time <- "29"
melt.Temp.T29$Treatment <- "Temperature"

melt.test %>%
  filter(Time == "35" & Pretreatment == "Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.Temp.T34
melt.Temp.T34$Time <- "34"
melt.Temp.T34$Treatment <- "Temperature"

rbind(melt.Temp.T9, melt.Temp.T14, melt.Temp.T19, melt.Temp.T24,
      melt.Temp.T29, melt.Temp.T34) -> temp.relabund


temp.relabund$mean <- (temp.relabund$mean)*100
temp.relabund$sd <- (temp.relabund$sd)*100

as.data.frame(ps.rare.top15@tax_table) -> GenusASV
data.frame(GenusASV$ASV, GenusASV$Genus) -> GenusASV
colnames(GenusASV) <- c("ASV", "Genus")

merge(temp.relabund, GenusASV, by = "ASV") -> temp.relabund
temp.relabund$GenusASV <- paste0(temp.relabund$Genus, " (", temp.relabund$ASV, ")")

## Plot
temp.relabund$Time <- factor(temp.relabund$Time, c("9", "14", "19","24",
                                                   "29", "34"), ordered = TRUE)

temp_MRA <- ggplot(data = temp.relabund, aes(x = Time, y = mean, fill = GenusASV)) + 
  geom_bar(stat = "identity") + facet_wrap(~Treatment)

temp_MRA <- temp_MRA + 
  xlab("Time (Days)") + ylab("Mean Relative Abundance (%)") + guides(fill = guide_legend(title = "Genus")) + 
  theme_bw() + theme(axis.text = element_text(face = "bold", size = 11.5), 
                     axis.title = element_text(face = "bold", size = 12), 
                     title = element_text(face = "bold"),
                     strip.text = element_text(face = "bold", size = 12),
                     strip.background = element_rect(
                       color = "black", # Border color
                       fill = "white", # Background fill color of the strip
                       size = 1.5)) +
  scale_fill_manual(values = colors)

## Antibiotics + Temp
melt.test %>%
  filter(Time == "2" & Treatment_Long == "Antibiotics + Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.AbxTemp.T2
melt.AbxTemp.T2$Time <- "2"
melt.AbxTemp.T2$Treatment <- "Antibiotics + Temperature"

melt.test %>%
  filter(Time == "5" & Treatment_Long == "Antibiotics + Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.AbxTemp.T4
melt.AbxTemp.T4$Time <- "4"
melt.AbxTemp.T4$Treatment <- "Antibiotics + Temperature"

melt.test %>%
  filter(Time == "10" & Treatment_Long == "Antibiotics + Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.AbxTemp.T9
melt.AbxTemp.T9$Time <- "9"
melt.AbxTemp.T9$Treatment <- "Antibiotics + Temperature"

melt.test %>%
  filter(Time == "15" & Treatment_Long == "Antibiotics + Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.AbxTemp.T14
melt.AbxTemp.T14$Time <- "14"
melt.AbxTemp.T14$Treatment <- "Antibiotics + Temperature"

melt.test %>%
  filter(Time == "20" & Treatment_Long == "Antibiotics + Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.AbxTemp.T19
melt.AbxTemp.T19$Time <- "19"
melt.AbxTemp.T19$Treatment <- "Antibiotics + Temperature"

melt.test %>%
  filter(Time == "25" & Treatment_Long == "Antibiotics + Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.AbxTemp.T24
melt.AbxTemp.T24$Time <- "24"
melt.AbxTemp.T24$Treatment <- "Antibiotics + Temperature"

melt.test %>%
  filter(Time == "30" & Treatment_Long == "Antibiotics + Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.AbxTemp.T29
melt.AbxTemp.T29$Time <- "29"
melt.AbxTemp.T29$Treatment <- "Antibiotics + Temperature"

melt.test %>%
  filter(Time == "35" & Treatment_Long == "Antibiotics + Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.AbxTemp.T34
melt.AbxTemp.T34$Time <- "34"
melt.AbxTemp.T34$Treatment <- "Antibiotics + Temperature"

rbind(melt.AbxTemp.T2, melt.AbxTemp.T4, melt.AbxTemp.T9, melt.AbxTemp.T14, melt.AbxTemp.T19, melt.AbxTemp.T24,
      melt.AbxTemp.T29, melt.AbxTemp.T34) -> abxtemp.relabund

abxtemp.relabund$mean <- (abxtemp.relabund$mean)*100
abxtemp.relabund$sd <- (abxtemp.relabund$sd)*100

as.data.frame(ps.rare.top15@tax_table) -> GenusASV
data.frame(GenusASV$ASV, GenusASV$Genus) -> GenusASV
colnames(GenusASV) <- c("ASV", "Genus")

merge(abxtemp.relabund, GenusASV, by = "ASV") -> abxtemp.relabund
abxtemp.relabund$GenusASV <- paste0(abxtemp.relabund$Genus, " (", abxtemp.relabund$ASV, ")")

## Plot
abxtemp.relabund$Time <- factor(abxtemp.relabund$Time, c("2", "4", "9", "14", "19","24",
                                                         "29", "34"), ordered = TRUE)

abxtemp_MRA <- ggplot(data = abxtemp.relabund, aes(x = Time, y = mean, fill = GenusASV)) + 
  geom_bar(stat = "identity") + facet_wrap(~Treatment)

abxtemp_MRA <- abxtemp_MRA + 
  xlab("Time (Days)") + ylab("Mean Relative Abundance (%)") + guides(fill = guide_legend(title = "Genus")) + 
  theme_bw() + theme(axis.text = element_text(face = "bold", size = 11.5), 
                     axis.title = element_text(face = "bold", size = 12), 
                     title = element_text(face = "bold"),
                     strip.text = element_text(face = "bold", size = 12),
                     strip.background = element_rect(
                       color = "black", # Border color
                       fill = "white", # Background fill color of the strip
                       size = 1.5)) +
  scale_fill_manual(values = colors)

## Combine data and plot
rbind(blank.relabund, abx.relabund, temp.relabund, abxtemp.relabund) -> relabund.All

write.csv(relabund.All, here::here("Output Files/05 - Relative Abundance - Output/mean_relabund.csv"))

relabund.All$Treatment <- factor(relabund.All$Treatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)

relabund.All$Time <- factor(relabund.All$Time, c("0","2", "4", "9", "14", "19","24",
                                                         "29", "34"), ordered = TRUE)


MRA.all <- ggplot(data = relabund.All, aes(x = Time, y = mean, fill = GenusASV)) + 
  geom_bar(stat = "identity") + facet_wrap(~Treatment)

MRA.all <- MRA.all + 
  xlab("Time (Days)") + ylab("Mean Relative Abundance (%)") + guides(fill = guide_legend(title = "ASV")) + 
  theme_bw() + theme(axis.text = element_text(face = "bold", size = 11.5), 
                     axis.title = element_text(face = "bold", size = 12), 
                     title = element_text(face = "bold"),
                     strip.text = element_text(face = "bold", size = 12),
                     strip.background = element_rect(
                       color = "black", # Border color
                       fill = "white", # Background fill color of the strip
                       size = 1.5)) +
  scale_fill_manual(values = colors)
ggplot2::ggsave(here::here("Output Files/05 - Relative Abundance - Output/01 - plot_meanRelAbund_all.png"), MRA.all,
                height = 450, width = 700, units = "mm",
                scale = 0.5, dpi = 1000)