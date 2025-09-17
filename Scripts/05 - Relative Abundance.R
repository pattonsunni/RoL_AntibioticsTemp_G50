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
library(vegan)

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
  filter(Time == "2" & Treatment_Long == "Antibiotics + Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.AbxTemp.T2
melt.AbxTemp.T2$Time <- "2"
melt.AbxTemp.T2$Treatment <- "Antibiotics"

melt.test %>%
  filter(Time == "5" & Treatment_Long == "Antibiotics + Temperature") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = ASV) %>%
  arrange() -> melt.AbxTemp.T4
melt.AbxTemp.T4$Time <- "4"
melt.AbxTemp.T4$Treatment <- "Antibiotics"


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

rbind(melt.AbxTemp.T9, melt.AbxTemp.T14, melt.AbxTemp.T19, melt.AbxTemp.T24,
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
  theme_bw(base_line_size = 1, base_rect_size = 1.5) + theme(axis.text = element_text(face = "bold", size = 11.5), 
                     axis.title = element_text(face = "bold", size = 12), 
                     title = element_text(face = "bold"),
                     strip.text = element_text(face = "bold", size = 12),
                     strip.background = element_rect(
                       color = "black", # Border color
                       fill = "white", # Background fill color of the strip
                       size = 1.5)) +
  scale_fill_manual(values = colors)
ggplot2::ggsave(here::here("Output Files/05 - Relative Abundance - Output/01 - plot_meanRelAbund_all.svg"), MRA.all,
                height = 450, width = 700, units = "mm",
                scale = 0.5, dpi = 1000)

## Plot relative abundance (Supplemental) ====
as.data.frame(ps.rare.top15.noASV1@tax_table) -> tax
tax$GenusASV <- paste0(tax$Genus, " (", tax$ASV, ")")
as.matrix(tax) -> tax
tax_table(ps.rare.top15.noASV1) <- tax

### Relabel sample names for clearer plotting/ordering
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T2_", "T02_", ps.rare.top15.noASV1@sam_data$CoralID) 
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T5_", "T04_", ps.rare.top15.noASV1@sam_data$CoralID) 
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T10_", "T09_", ps.rare.top15.noASV1@sam_data$CoralID) 
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T15_", "T14_", ps.rare.top15.noASV1@sam_data$CoralID) 
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T20_", "T19_", ps.rare.top15.noASV1@sam_data$CoralID) 
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T25_", "T24_", ps.rare.top15.noASV1@sam_data$CoralID) 
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T30_", "T29_", ps.rare.top15.noASV1@sam_data$CoralID) 
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T35_", "T34_", ps.rare.top15.noASV1@sam_data$CoralID) 

## Set colors
colors.noASV1 <- c("#66A61EFF", "#f38400", "#2D2651FF", "#C969A1FF", "#be0032",
                   "#dcd300", "#654522","#e25822", "#5D7298FF","#666666FF",
                   "#a1caf1", "#81B28DFF", "#A6761DFF", "#f38400", "#c2b280")


## Plot individual treatments
### no treatment
plot.relAbund.blank.noASV1 <- plot_bar(subset_samples(ps.rare.top15.noASV1, ps.rare.top15.noASV1@sam_data$Pretreatment == "No Treatment"), x="CoralID", fill="GenusASV")
plot.relAbund.blank.noASV1 <- plot.relAbund.blank.noASV1 + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 20), 
        axis.text.y = element_text(face = "bold", size = 20), 
        axis.title.y = element_text(face = "bold", size = 26), 
        axis.title.x = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1),
        strip.text = element_text(face = "bold", size = 26)) + xlab("Samples") +
  facet_wrap(~Pretreatment) + scale_fill_manual(values = colors.noASV1)  + scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + ylab("Relative Abundance (%)")
### antibiotics
plot.relAbund.abx.noASV1 <- plot_bar(subset_samples(ps.rare.top15.noASV1, ps.rare.top15.noASV1@sam_data$Pretreatment == "Antibiotics"), x="CoralID", fill="GenusASV")

plot.relAbund.abx.noASV1 <- plot.relAbund.abx.noASV1 + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 20), 
        axis.text.y = element_text(face = "bold", size = 20), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.text = element_text(face = "bold", size = 20),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1),
        strip.text = element_text(face = "bold", size = 26)) + xlab("Samples") +
  facet_wrap(~Pretreatment) + scale_fill_manual(values = colors.noASV1)  + scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""))
## temperature
plot.relAbund.temp.noASV1 <- plot_bar(subset_samples(ps.rare.top15.noASV1, ps.rare.top15.noASV1@sam_data$Pretreatment == "Temperature"), x="CoralID", fill="GenusASV")

plot.relAbund.temp.noASV1 <- plot.relAbund.temp.noASV1 + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 20), 
        axis.text.y = element_text(face = "bold", size = 20), 
        axis.title.y = element_text(face = "bold", size = 26), 
        axis.title.x = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1),
        strip.text = element_text(face = "bold", size = 26)) + xlab("Samples") +
  facet_wrap(~Pretreatment)+ scale_fill_manual(values = colors.noASV1)  + scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + ylab("Relative Abundance (%)")
## antibiotics + temperature
plot.relAbund.tempabx.noASV1 <- plot_bar(subset_samples(ps.rare.top15.noASV1, ps.rare.top15.noASV1@sam_data$Pretreatment == "Antibiotics + Temperature"), x="CoralID", fill="GenusASV")

plot.relAbund.tempabx.noASV1 <- plot.relAbund.tempabx.noASV1 + theme_bw() + guides(fill = guide_legend(title = "ASV", face = bold)) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 20), 
        axis.text.y = element_text(face = "bold", size = 20), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1),
        strip.text = element_text(face = "bold", size = 26)) + xlab("Samples") +
  facet_wrap(~Pretreatment)+ scale_fill_manual(values = colors.noASV1) + scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""), limits = c(0,1)) 

ggarrange(plot.relAbund.blank.noASV1, plot.relAbund.abx.noASV1, plot.relAbund.temp.noASV1, plot.relAbund.tempabx.noASV1, legend = "none") -> relabund_noASV1

ggplot2::ggsave(here::here("Output Files/05 - Relative Abundance - Output/plot_RelAbund_noASV1.png"), relabund_noASV1,
                height = 600, width = 1200, units = "mm",
                scale = 0.5, dpi = 1000)


## Relative abundance when aquarickettsia removed ====
### load data
ps.pruned <- readRDS(here::here("Output Files/04 - Phyloseq - Output/ps.pruned.RDS")) ## pruned data
### remove aquarickettsia
ps.pruned.noASV1 <- subset_taxa(ps.pruned, ASV != "ASV1") # Went from 12060495 reads to 781573 reads
sort(sample_sums(ps.pruned.noASV1)) # Lowest sample has 183 reads (G50_T15_1_c), highest sample has 25,795 reads (G50_T2_0_c)

### Rarefaction curve 
as.matrix(as.data.frame(ps.pruned.noASV1@otu_table)) -> data
rarecurve(data, step = 100, label = FALSE)

rarefied_df <- rrarefy(data, sample = 1145) ## This doesn't capture really any of the diversity; so I don't think it's going to be the best idea
sort(rowSums(rarefied_df)) 

### Make new phyloseq object after subsampling to even depth 
rare_samData <- ps.pruned.noASV1@sam_data
rare_taxTable <- ps.pruned.noASV1@tax_table
rare_otuTable <- rarefied_df

### Remove unrarefied samples from phyloseq obj
ps.noASV1.rare <- phyloseq(otu_table(rare_otuTable, taxa_are_rows=FALSE),sample_data(rare_samData),tax_table(rare_taxTable))

tax_fix(
  ps.noASV1.rare,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
) -> ps.noASV1.rare

subset_samples(ps.noASV1.rare, Samples != "G50_T15_1_c") -> ps.noASV1.rare
subset_samples(ps.noASV1.rare, Samples != "G50_T10_3_a") -> ps.noASV1.rare
subset_samples(ps.noASV1.rare, Samples != "G50_T10_3_d") -> ps.noASV1.rare
subset_samples(ps.noASV1.rare, Samples != "G50_T10_3_c") -> ps.noASV1.rare
subset_samples(ps.noASV1.rare, Samples != "G50_T15_1_d") -> ps.noASV1.rare
subset_samples(ps.noASV1.rare, Samples != "G50_T35_3_d") -> ps.noASV1.rare
subset_samples(ps.noASV1.rare, Samples != "G50_T15_3_c") -> ps.noASV1.rare
subset_samples(ps.noASV1.rare, Samples != "G50_T20_3_b") -> ps.noASV1.rare
subset_samples(ps.noASV1.rare, Samples != "G50_T20_1_b") -> ps.noASV1.rare

write_rds(ps.noASV1.rare, here::here("Output Files/05 - Relative Abundance - Output/ps.noASV1.rare.rds"))

ps.noASV1.rare.trans <- transform_sample_counts(ps.noASV1.rare, function(OTU) OTU/sum(OTU))

## Identify top 15 most abundant taxa
top15_noASV1 <- names(sort(taxa_sums(ps.noASV1.rare.trans), decreasing = TRUE))[1:15] 

## Prune top 15
ps.rare.top15.noASV1 <- prune_taxa(top15_noASV1, ps.noASV1.rare.trans)

## subset by phase
# Add phase column
ps.rare.top15.noASV1@sam_data$Phase <- paste0(ps.rare.top15.noASV1@sam_data$Time_TotalDays)

ps.rare.top15.noASV1@sam_data$Phase[ps.rare.top15.noASV1@sam_data$Phase == "0"] <- "Antibiotic challenge" 
ps.rare.top15.noASV1@sam_data$Phase[ps.rare.top15.noASV1@sam_data$Phase == "2"] <- "Antibiotic challenge"  
ps.rare.top15.noASV1@sam_data$Phase[ps.rare.top15.noASV1@sam_data$Phase == "4"] <- "Antibiotic challenge"  
ps.rare.top15.noASV1@sam_data$Phase[ps.rare.top15.noASV1@sam_data$Phase == "9"] <- "Temperature ramp"
ps.rare.top15.noASV1@sam_data$Phase[ps.rare.top15.noASV1@sam_data$Phase == "14"] <- "Temperature ramp"  
ps.rare.top15.noASV1@sam_data$Phase[ps.rare.top15.noASV1@sam_data$Phase == "19"] <- "Temperature challenge"  
ps.rare.top15.noASV1@sam_data$Phase[ps.rare.top15.noASV1@sam_data$Phase == "24"] <- "Temperature challenge"  
ps.rare.top15.noASV1@sam_data$Phase[ps.rare.top15.noASV1@sam_data$Phase == "29"] <- "Temperature challenge"  
ps.rare.top15.noASV1@sam_data$Phase[ps.rare.top15.noASV1@sam_data$Phase == "34"] <- "Temperature challenge" 

as.data.frame(ps.rare.top15.noASV1@tax_table) -> tax
tax$GenusASV <- paste0(tax$Genus, " (", tax$ASV, ")")
as.matrix(tax) -> tax
tax_table(ps.rare.top15.noASV1) <- tax

ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T35", "T34", ps.rare.top15.noASV1@sam_data$CoralID, ignore.case = FALSE)
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T30", "T29", ps.rare.top15.noASV1@sam_data$CoralID, ignore.case = FALSE)
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T25", "T24", ps.rare.top15.noASV1@sam_data$CoralID, ignore.case = FALSE)
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T20", "T19", ps.rare.top15.noASV1@sam_data$CoralID, ignore.case = FALSE)
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T15", "T14", ps.rare.top15.noASV1@sam_data$CoralID, ignore.case = FALSE)
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T10", "T9", ps.rare.top15.noASV1@sam_data$CoralID, ignore.case = FALSE)
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T5", "T4", ps.rare.top15.noASV1@sam_data$CoralID, ignore.case = FALSE)

ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T9", "T09", ps.rare.top15.noASV1@sam_data$CoralID, ignore.case = FALSE)
ps.rare.top15.noASV1@sam_data$CoralID <- gsub("T4", "T04", ps.rare.top15.noASV1@sam_data$CoralID, ignore.case = FALSE)


## antibiotic challenge
plot.abx.no.noASV1 <- plot_bar(subset_samples(ps.rare.top15.noASV1, (ps.rare.top15.noASV1@sam_data$Phase == "Antibiotic challenge") & ps.rare.top15.noASV1@sam_data$Pretreatment == "No Treatment"), x="CoralID", fill="GenusASV")
plot.abx.no.noASV1 <- plot.abx.no.noASV1 + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 16), 
        axis.text.y = element_text(face = "bold", size = 16), 
        axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1),
        strip.text = element_text(face = "bold", size = 16)) + xlab("Samples") +
  facet_wrap(~Pretreatment) + scale_fill_manual(values = colors.noASV1) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""), limits = c(0,1)) + ylab("Relative Abundance (%)")

plot.abx.abx.noASV1 <- plot_bar(subset_samples(ps.rare.top15.noASV1, (ps.rare.top15.noASV1@sam_data$Phase == "Antibiotic challenge") & ps.rare.top15.noASV1@sam_data$Pretreatment == "Antibiotics"), x="CoralID", fill="GenusASV")
plot.abx.abx.noASV1 <- plot.abx.abx.noASV1 + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 16), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 14),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1),
        strip.text = element_text(face = "bold", size = 16)) + xlab("Samples") +
  facet_grid(Phase~Pretreatment) + scale_fill_manual(values = colors.noASV1)  + scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + ylab("Relative Abundance (%)")


ggarrange(plot.abx.no.noASV1, plot.abx.abx.noASV1, legend = "none") -> plot.abx.challenge


### ramp
plot.ramp.no.noASV1 <- plot_bar(subset_samples(ps.rare.top15.noASV1, (ps.rare.top15.noASV1@sam_data$Phase == "Temperature ramp") & ps.rare.top15.noASV1@sam_data$Pretreatment == "No Treatment"), x="CoralID", fill="GenusASV")
plot.ramp.no.noASV1 <- plot.ramp.no.noASV1 + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 16), 
        axis.text.y = element_text(face = "bold", size = 16), 
        axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1),
        strip.text = element_text(face = "bold", size = 16)) + xlab("Samples") +
  facet_wrap(~Pretreatment) + scale_fill_manual(values = colors.noASV1) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""), limits = c(0,1)) + ylab("Relative Abundance (%)")

plot.ramp.abx.noASV1 <- plot_bar(subset_samples(ps.rare.top15.noASV1, (ps.rare.top15.noASV1@sam_data$Phase == "Temperature ramp") & ps.rare.top15.noASV1@sam_data$Pretreatment == "Antibiotics"), x="CoralID", fill="GenusASV")
plot.ramp.abx.noASV1 <- plot.ramp.abx.noASV1 + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 16), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1),
        strip.text = element_text(face = "bold", size = 16)) + xlab("Samples") +
  facet_wrap(~Pretreatment) + scale_fill_manual(values = colors.noASV1) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""), limits = c(0,1)) + ylab("Relative Abundance (%)")

plot.ramp.temp.noASV1 <- plot_bar(subset_samples(ps.rare.top15.noASV1, (ps.rare.top15.noASV1@sam_data$Phase == "Temperature ramp") & ps.rare.top15.noASV1@sam_data$Pretreatment == "Temperature"), x="CoralID", fill="GenusASV")
plot.ramp.temp.noASV1 <- plot.ramp.temp.noASV1 + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 16), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1),
        strip.text = element_text(face = "bold", size = 16)) + xlab("Samples") +
  facet_wrap(~Pretreatment) + scale_fill_manual(values = colors.noASV1) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""), limits = c(0,1)) + ylab("Relative Abundance (%)")


plot.ramp.abxtemp.noASV1 <- plot_bar(subset_samples(ps.rare.top15.noASV1, (ps.rare.top15.noASV1@sam_data$Phase == "Temperature ramp") & ps.rare.top15.noASV1@sam_data$Pretreatment == "Antibiotics + Temperature"), x="CoralID", fill="GenusASV")
plot.ramp.abxtemp.noASV1 <- plot.ramp.abxtemp.noASV1 + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 16), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 14),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1),
        strip.text = element_text(face = "bold", size = 16)) + xlab("Samples") +
  facet_grid(Phase~Pretreatment) + scale_fill_manual(values = colors.noASV1)  + scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + ylab("Relative Abundance (%)")

ggarrange(plot.ramp.no.noASV1, plot.ramp.abx.noASV1, plot.ramp.temp.noASV1, plot.ramp.abxtemp.noASV1, ncol = 4, legend = "none") -> plot.ramp


### temp
plot.temp.no.noASV1 <- plot_bar(subset_samples(ps.rare.top15.noASV1, (ps.rare.top15.noASV1@sam_data$Phase == "Temperature challenge") & ps.rare.top15.noASV1@sam_data$Pretreatment == "No Treatment"), x="CoralID", fill="GenusASV")
plot.temp.no.noASV1 <- plot.temp.no.noASV1 + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 16), 
        axis.text.y = element_text(face = "bold", size = 16), 
        axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1),
        strip.text = element_text(face = "bold", size = 16)) + xlab("Samples") +
  facet_wrap(~Pretreatment) + scale_fill_manual(values = colors.noASV1) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""), limits = c(0,1)) + ylab("Relative Abundance (%)")

plot.temp.abx.noASV1 <- plot_bar(subset_samples(ps.rare.top15.noASV1, (ps.rare.top15.noASV1@sam_data$Phase == "Temperature challenge") & ps.rare.top15.noASV1@sam_data$Pretreatment == "Antibiotics"), x="CoralID", fill="GenusASV")
plot.temp.abx.noASV1 <- plot.temp.abx.noASV1 + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 16), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1),
        strip.text = element_text(face = "bold", size = 16)) + xlab("Samples") +
  facet_wrap(~Pretreatment) + scale_fill_manual(values = colors.noASV1) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""), limits = c(0,1)) + ylab("Relative Abundance (%)")

plot.temp.temp.noASV1 <- plot_bar(subset_samples(ps.rare.top15.noASV1, (ps.rare.top15.noASV1@sam_data$Phase == "Temperature challenge") & ps.rare.top15.noASV1@sam_data$Pretreatment == "Temperature"), x="CoralID", fill="GenusASV")
plot.temp.temp.noASV1 <- plot.temp.temp.noASV1 + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 16), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1),
        strip.text = element_text(face = "bold", size = 16)) + xlab("Samples") +
  facet_wrap(~Pretreatment) + scale_fill_manual(values = colors.noASV1) + 
  scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""), limits = c(0,1)) + ylab("Relative Abundance (%)")


plot.temp.abxtemp.noASV1 <- plot_bar(subset_samples(ps.rare.top15.noASV1, (ps.rare.top15.noASV1@sam_data$Phase == "Temperature challenge") & ps.rare.top15.noASV1@sam_data$Pretreatment == "Antibiotics + Temperature"), x="CoralID", fill="GenusASV")
plot.temp.abxtemp.noASV1 <- plot.temp.abxtemp.noASV1 + theme_bw() + guides(fill = guide_legend(title = "ASV")) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 16), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 14),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1),
        strip.text = element_text(face = "bold", size = 16)) + xlab("Samples") +
  facet_grid(Phase~Pretreatment) + scale_fill_manual(values = colors.noASV1)  + scale_y_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) + ylab("Relative Abundance (%)")

ggarrange(plot.temp.no.noASV1, plot.temp.abx.noASV1, plot.temp.temp.noASV1, plot.temp.abxtemp.noASV1, ncol = 4, legend = "none") -> plot.temp

ggarrange(plot.abx.challenge, plot.ramp, plot.temp, nrow = 3, common.legend = TRUE) -> x


ggplot2::ggsave(here::here("Output Files/05 - Relative Abundance - Output/plot_relabund_noASV1_new.svg"), x,
                height = 610, width = 750, units = "mm",
                scale = 0.5, dpi = 1000)
