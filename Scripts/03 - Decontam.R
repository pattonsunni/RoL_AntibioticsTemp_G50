# Author: Sunni Patton
# Last edited: 10/11/24
# Title: Running decontam

# Code largely adapted from https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

## Set seed ====
set.seed(123)

## Load libraries ====
library(decontam)
library(ggplot2)
library(phyloseq)

## Load data ====
readRDS(here::here("Output files/02 - Phyloseq Preprocessing - Output/ps.All.rds")) -> ps.All
## Plot sample library size ====
df <- as.data.frame(sample_data(ps.All)) 
df$LibrarySize <- sample_sums(ps.All)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
plot_decontam <- ggplot(data=df, aes(x = Index, y = LibrarySize, color = Sample_control)) + geom_point()

## Save plot ====
ggplot2::ggsave(here::here("Output Files/03 - Decontam - Output/plot_decontam.png"), plot_decontam,
                height = 200, width = 350, units = "mm",
                scale = 0.5, dpi = 1000)

## Run decontam ====
sample_data(ps.All)$is.neg <- sample_data(ps.All)$Sample_control == "Negative control"

contamdf_combined <- isContaminant(ps.All, method="combined", neg="is.neg", conc="quant_reading", threshold=0.5)
table(contamdf_combined$contaminant) #identified 168 as contaminants 
head(which(contamdf_combined$contaminant)) 

## Make presence/absence plot ====

### Make phyloseq object of presence-absence in negative controls and true samples
ps.presAbs <- transform_sample_counts(ps.All, function(abund) 1*(abund>0))
ps.presAbs.neg <- prune_samples(sample_data(ps.presAbs)$Sample_control == "Negative control", ps.presAbs)
ps.presAbs.pos <- prune_samples(sample_data(ps.presAbs)$Sample_control == "True sample", ps.presAbs)

# Make data.frame of prevalence in positive and negative samples
df.presAbs <- data.frame(presAbs.pos=taxa_sums(ps.presAbs.pos), presAbs.neg=taxa_sums(ps.presAbs.neg),
                         contaminant=contamdf_combined$contaminant)
plot.presAbs <- ggplot(data=df.presAbs, aes(x=presAbs.neg, y=presAbs.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence/Frequency (Negative Controls)") + ylab("Prevalence/Frequency (True Samples)")
# Looks like there are two taxa with somewhat high prevalence in negative controls that aren't considered contaminants -> need to look into these

## Save plot ====
ggplot2::ggsave(here::here("Output Files/03 - Decontam - Output/plot_presAbs.png"), plot.presAbs,
                height = 200, width = 200, units = "mm",
                scale = 0.5, dpi = 1000)

## Remove contaminants ====
### Make phyloseq object for contaminants only 
ps.contam <- prune_taxa(contamdf_combined$contaminant, ps.All) # Contains 168 ASVs (as expected)
View(ps.contam@tax_table) # identified an Aquarickettsia ASV as a contaminant 
ps.contam <- prune_taxa(taxa_sums(ps.contam@otu_table) > 0, ps.contam)
summary(taxa_sums(ps.contam))

## Make phyloseq object with contaminants removed
ps.nocontam <- prune_taxa(!contamdf_combined$contaminant, ps.All) 
ps.nocontam <- prune_taxa(taxa_sums(ps.nocontam@otu_table) > 0, ps.nocontam) # Contains 4739 ASVs after contaminants removed
summary(taxa_sums(ps.nocontam)) # Least abundant taxon/taxa has 2 reads overall, highest abundance taxon has 11,278,922

sum(ps.nocontam@otu_table) # 12,065,466 reads
sum(ps.All@otu_table) - sum(ps.nocontam@otu_table) # 9,109 reads were contaminants
saveRDS(ps.nocontam, here::here("Output Files/03 - Decontam - Output/ps.nocontam.rds"))

## Make phyloseq object with NCs removed ====
### Make phyloseq object with just NCs (for ps.All and ps.nocontam)
subset_samples(ps.nocontam, Treatment != "Negative Control") -> ps.noneg

subset_samples(ps.nocontam, Treatment != "Negative Control") -> ps.noneg

### Phyloseq object summary
ps.noneg <- prune_taxa(taxa_sums(ps.noneg@otu_table) > 0, ps.noneg) 

summary(taxa_sums(ps.noneg@otu_table))
sort(sample_sums(ps.noneg@otu_table)) 

saveRDS(ps.noneg, here::here("Output Files/03 - Decontam - Output/ps.noneg.rds"))
