# Author: Sunni Patton
# Last edited: 10/17/2024
# Title: Read preprocessing

## Load libraries ====
library(here)
library(dada2)
library(dplyr)

## Set seed for reproducibility ====
set.seed(123)

## Pre-Processing ====

## Set path to raw files
### Primers have been removed from these reads
path <- here::here("Data (post-cutadapt)")

### Set path for forward and reverse reads
fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

## Extract sample names from files
### Assume name follows: lane1-s0001-index--Findenx-Rindex-SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(basename(fnFs), "-"), `[`,7) 
sampleNames <- sapply(strsplit(basename(sampleNames), "\\."), `[`,1) 
sampleNames <- gsub("_S\\d+\\d?\\d?\\d?_R1_001", "", sampleNames) 

## Set file destination for files after quality filtering
filtFs <- file.path(path, "filterAndTrim", paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filterAndTrim", paste0(sampleNames, "_R_filt.fastq.gz"))

## Filter and Trim ====
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(250, 200),
                        maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=FALSE)
## Save output
saveRDS(out, here::here("Output Files/00 - Read Preprocessing - Output/out.rds"))

## Assess number of reads lost
sum(out[,1])-sum(out[,2]) # 13388056 initially, 12608343 remaining, 779713 lost

## Dereplicate reads
derep_F <- derepFastq(filtFs)
derep_R <- derepFastq(filtRs)

## Learn errors and infer sample sequence ====
errF <- learnErrors(derep_F, multithread = FALSE, nbases = 5e8)  
saveRDS(errF, here::here("Output Files/00 - Read Preprocessing - Output/errF.rds"))

errR <- learnErrors(derep_R,multithread = FALSE, nbases = 5e8) 
saveRDS(errR, here::here("Output Files/00 - Read Preprocessing - Output/errR.rds"))

## Ensure sample naming is consistent
names(derep_F)<-sampleNames
names(derep_R)<-sampleNames

## Infer sample sequence
dadaForward <- dada(derep_F, err=errF, multithread=FALSE)
saveRDS(dadaForward, here::here("Output Files/00 - Read Preprocessing - Output/dadaForward.rds"))

dadaReverse <- dada(derep_R, err=errR, multithread=FALSE)
saveRDS(dadaReverse, here::here("Output Files/00 - Read Preprocessing - Output/dadaReverse.rds"))

## Create contigs and sequence table ====
contigs <- mergePairs(dadaForward, derep_F, dadaReverse, derep_R)
saveRDS(contigs, here::here("Output Files/00 - Read Preprocessing - Output/contigs.rds"))

## Make sequence table and visualize contig length and frequency
seq_table <- makeSequenceTable(contigs) 
dim(seq_table)
table(nchar(getSequences(seq_table)))

### Keep contigs within desired size range
seq_table <- seq_table[,nchar(colnames(seq_table)) %in% 250:256]
table(nchar(getSequences(seq_table))) 
dim(seq_table) 
sum(seq_table) 

## Save output
saveRDS(seq_table, here::here("Output Files/00 - Read Preprocessing - Output/seq_table.rds"))

## Remove chimeras ====
seq_table_nochimeri <- removeBimeraDenovo(seq_table, method="consensus", multithread=TRUE, verbose=TRUE) 
dim(seq_table_nochimeri) 
sum(seq_table) - sum(seq_table_nochimeri) 

saveRDS(seq_table_nochimeri, here::here("Output Files/00 - Read Preprocessing - Output/seq_table_nochimeri.rds"))

## Assign taxonomy ====
taxa <- assignTaxonomy(seq_table_nochimeri, here::here("silva_nr99_v138.1_train_set.fa.gz"), multithread=FALSE, minBoot = 80)
taxa <- addSpecies(taxa, here::here("silva_species_assignment_v138.1.fa.gz"))

saveRDS(taxa, here::here("Output Files/00 - Read Preprocessing - Output/taxa.rds")) 

## Remove off-target sequences ====
### New sequence table (no chloroplast)
is.chloroplast <- taxa[,"Order"] %in% "Chloroplast"
seq_table_nochloro <- seq_table_nochimeri[,!is.chloroplast]
dim(seq_table_nochloro) #5131 taxa
sum(seq_table_nochimeri) - sum(seq_table_nochloro) 
saveRDS(seq_table_nochloro, here::here("Output Files/00 - Read Preprocessing - Output/seq_table_nochloro.rds"))

### New taxonomy table (no chloroplast)
taxonomy_nochloro <- taxa[!is.chloroplast,]
dim(taxonomy_nochloro)
View(taxonomy_nochloro)
saveRDS(taxonomy_nochloro, here::here("Output Files/00 - Read Preprocessing - Output/taxonomy_nochlor.rds"))

### New sequence table (no mitochondria or chloroplast)
is.mitochondria <- taxonomy_nochloro[,"Family"] %in% "Mitochondria"
seq_table_nomito <- seq_table_nochloro[,!is.mitochondria]
dim(seq_table_nomito) #5004 taxa
sum(seq_table_nochloro) - sum(seq_table_nomito) 
saveRDS(seq_table_nomito, here::here("Output Files/00 - Read Preprocessing - Output/seq_table_nomito.rds"))

### New taxonomy table (no mitochondria or chloroplast)
taxonomy_nomito <- taxonomy_nochloro[!is.mitochondria,]
dim(taxonomy_nomito)
View(taxonomy_nomito)
saveRDS(taxonomy_nomito, here::here("Output Files/00 - Read Preprocessing - Output/taxonomy_nomito.rds"))

## Remove sequences not annotated beyond the Kingdom level
is.NA <- taxonomy_nomito[,"Phylum"] %in% NA
seq_table_noNA <- seq_table_nomito[,!is.NA]
dim(seq_table_noNA) 
sum(seq_table_nomito) - sum(seq_table_noNA) 
saveRDS(seq_table_noNA, here::here("Output Files/00 - Read Preprocessing - Output/seq_table_noNA.rds"))

taxonomy_noNA <- taxonomy_nomito[!is.NA,]
dim(taxonomy_noNA)
saveRDS(taxonomy_noNA, here::here("Output Files/00 - Read Preprocessing - Output/taxonomy_noNA.rds"))

## Assess total number of reads removed and ASV frequency ====
sum(seq_table_noNA) 
dim(seq_table_noNA) 

summary(colSums(seq_table_noNA)) 
summary(rowSums(seq_table_noNA))  

sort(rowSums(seq_table_noNA)) 

## Tracking reads removed at each step ====
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaForward, getN), sapply(dadaReverse, getN), sapply(contigs, getN), rowSums(seq_table_nochimeri), rowSums(seq_table_nochloro), rowSums(seq_table_nomito), rowSums(seq_table_noNA))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "nochloro", "nomito", "noNA.phylum")

### Update sample names based on track_G7 rows (simply forcing column names to be sampleNames incorrectly assigns the names)
sampleNames_new <- sapply(strsplit(rownames(track), "-"), `[`,7) 
sampleNames_new <- sapply(strsplit(basename(sampleNames_new), "\\."), `[`,1) 
sampleNames_new <- gsub("_S\\d+\\d?\\d?\\d?_R1_001", "", sampleNames_new) 

### Assign new row names (shorter sample name) using sampleNames_new
rownames(track) <- sampleNames_new

write.csv(track, here::here("Output Files/00 - Read Preprocessing - Output/track.csv")) 