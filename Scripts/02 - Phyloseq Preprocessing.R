# Author: Sunni Patton
# Last edited: 7/10/24
# Title: Creating initial phyloseq object

## Set seed ====
set.seed(123)

## Load libraries ====
library(phyloseq)

## Load metadata file ====
read.csv(here::here("Output Files/01 - Metadata - Output/metadata_quant.csv")) -> metadata_quant

## Load sequence table and taxonomy table ====
readRDS(here::here("Output Files/00 - Read Preprocessing - Output/seq_table_noNA.rds")) -> seq_table_noNA
readRDS(here::here("Output Files/00 - Read Preprocessing - Output/taxonomy_noNA.rds")) -> taxonomy_noNA

## Make sample numbers into sample names to work with phyloseq ====
x <- metadata_quant$Samples
rownames(metadata_quant) <- x

## Create phyloseq object ====
ps.All <- phyloseq(otu_table(seq_table_noNA, taxa_are_rows=FALSE), sample_data(metadata_quant), tax_table(taxonomy_noNA))
saveRDS(ps.All, here::here("Output Files/02 - Phyloseq Preprocessing - Output/ps.All.rds"))

## Inspect phyloseq object ====
ps.All 
# Taxa distribution 
summary(taxa_sums(ps.All@otu_table)) 
# Sample read distribution
summary(sample_sums(ps.All@otu_table)) 

