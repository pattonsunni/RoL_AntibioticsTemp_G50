# Author: Sunni Patton
# Last edited: 10/18/24 
# Title: Sample processing in phyloseq
# Overview: Adding phylogenetic tree to phyloseq object and randomly subsampling to an even depth

## Set seed ====
set.seed(123)

## Load libraries ====
library(speedyseq)
library(tidyr)
library(phytools)
library(seqateurs)
library(vegan)
library(microViz)
library(ips)
library(ggplot2)

## Load phyloseq object (contaminants and NCs removed) ====
readRDS(here::here("Output Files/03 - Decontam - Output/ps.noneg.rds")) -> ps.noneg

## Add ASV to taxonomy table ====
ps.noneg <- ps.noneg %>% mutate_tax_table(ASV = paste0("ASV", 1:4709),
                                          Sequence = paste0(rownames(ps.noneg@tax_table)))

## Change DNA sequence to ASV number ====
dna <- Biostrings::DNAStringSet(taxa_names(ps.noneg))
names(dna) <- taxa_names(ps.noneg)
ps.noneg <- merge_phyloseq(ps.noneg, dna)
taxa_names(ps.noneg) <- paste0("ASV", seq(ntaxa(ps.noneg)))

saveRDS(ps.noneg, here::here("Output Files/04 - Phyloseq - Output/ps.noneg.rds"))

## Make rarefaction curve ====
as.matrix(as.data.frame(ps.noneg@otu_table)) -> data
rarecurve(data, step = 100, label = FALSE)

## Remove low abundance taxa ====
### Look at abundance and frequency of each taxon
### Doing this makes it easy to see in how many samples each taxon appears (Prevalence), and at what abundance (TotalCounts)
pst <- fast_melt(ps.noneg) # Takes a phyloseq object and returns a melted dataframe with sample ID, variable name, and abundance value for each observation

prevdt <- pst[, list(Prevalence = sum(count > 0), 
                     TotalCounts = sum(count)),
              by = taxaID]
View(prevdt) # Handful of high abundance and frequency taxa, but there a lot of taxa that are only present in one sample and/or present at low reads

### Prune taxa that are only present in 1 sample and with total reads < 1st quartile (5) ====
summary(taxa_sums(ps.noneg@otu_table))
keepTaxa <- prevdt[((Prevalence >1 | TotalCounts >5)), taxaID]
ps.pruned <- prune_taxa(keepTaxa, ps.noneg) # 3230 taxa remain after pruning

summary(taxa_sums(ps.pruned@otu_table))
# min: 4, 1st quart: 8, median: 15, mean: 3743, 3rd quart: 38, max: 11278922

sum(ps.noneg@otu_table) - sum(ps.pruned@otu_table)

saveRDS(ps.pruned, here::here("Output Files/04 - Phyloseq - Output/ps.pruned.rds"))

## Check library size by sample ====
ps.pruned@sam_data$Treatment_Long <- factor(ps.pruned@sam_data$Treatment_Long, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)


plot_libSize <- ggplot(ps.pruned@sam_data, aes(x = as.factor(Time), y = LibrarySize_PostQC, color = Treatment_Long)) + 
  scale_color_manual(values = c("black", "#5386B3", "#EF8737", "darkred")) +
  geom_boxplot(lwd = 1.1, outlier.color = "NA") 

plot_libSize <- plot_libSize + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Treatment_Long, color = Pretreatment), size = 3,
                                            position = position_dodge(width = 0.9))
plot_libSize <- plot_libSize + geom_point(aes(color = Treatment_Long), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold"))
plot_libSize <- plot_libSize + xlab("Treatment") + ylab("Number of Reads") 
plot_libSize <- plot_libSize + guides(color = guide_legend(title = "Treatment"))


## Randomly subsample to even depth (rrarefy) ====
### rrarefy requires a dataframe, so save ps.pruned out_table as df
otu_table <- as.matrix(as.data.frame(ps.pruned@otu_table))
View(otu_table)

sort(sample_sums(ps.pruned))
rarefied_df <- rrarefy(otu_table, sample = 52686)
sort(rowSums(rarefied_df)) 

### Make new phyloseq object after subsampling to even depth 
rare_samData <- ps.pruned@sam_data
rare_taxTable <- ps.pruned@tax_table
rare_otuTable <- rarefied_df

ps.rare <- phyloseq(otu_table(rare_otuTable, taxa_are_rows=FALSE),sample_data(rare_samData),tax_table(rare_taxTable))

### Make sure there aren't any taxa present that aren't actually in any sample after randomly subsampling
ps.rare <- prune_taxa(taxa_sums(ps.rare@otu_table) > 0, ps.rare)

summary(taxa_sums(ps.rare))


saveRDS(ps.rare, here::here("Output Files/04 - Phyloseq - Output/ps.rare.rds"))