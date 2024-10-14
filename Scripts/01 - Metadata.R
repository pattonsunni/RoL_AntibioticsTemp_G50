# Author: Sunni Patton
# Last edited: 10/11/24
# Title: Creating metadata file

## Set seed ====
set.seed(123)

## Load libraries ====
library(ggplot2)
library(dplyr)
library(rstatix)
library(broom)

## Extract sample name data ====
### Read in track if not already loaded (contains sample names)
read.csv(here::here("Output Files/00 - Read Preprocessing - Output/track.csv")) -> track


### Extract information 
time <- as.character(sapply(strsplit(rownames(track), "_T"), `[`,2))
time <- as.character(sapply(strsplit(time, "_"), `[`, 1)) 

treatment <- as.character(sapply(strsplit(rownames(track), "G50_T\\d+\\_"), `[`, 2)) 
treatment <- as.character(sapply(strsplit(treatment, "_"), `[`, 1))

replicate <- as.character(sapply(strsplit(rownames(track), "G50_T\\d+\\_"), `[`, 2))
replicate <- as.character(sapply(strsplit(replicate, "_"), `[`, 2))


### Make dataframe
metadata <- data.frame(rownames(track))

metadata$Time <- time
metadata$Treatment <- treatment
metadata$TankReplicate <- replicate

colnames(metadata)[1] <- "Samples"

metadata$Treatment[metadata$Samples == "NC_1" | 
                     metadata$Samples =="NC_2" |
                     metadata$Samples == "NC_3" |
                     metadata$Samples == "NC_4"] <- "Negative Control"

metadata$Treatment2 <- paste0(metadata$Treatment)

metadata$Treatment2[metadata$Treatment == "0"] <- "Blank"
metadata$Treatment2[metadata$Treatment == "1"] <- "Antibiotic"
metadata$Treatment2[metadata$Treatment == "2"] <- "Temperature"
metadata$Treatment2[metadata$Treatment == "3"] <- "Antibiotic + Temperature"

metadata$Pretreatment <- paste0(metadata$Treatment2)
metadata$Pretreatment[metadata$Time == "0"] <- "Pretreatment"

# Fix time column (end of antibiotic experiment was 96 hours, i.e., 4 days -- not 5. So 5 should be 4, and all times after that should be one less)
metadata$Time2 <- paste0(metadata$Time)
metadata$Time2[metadata$Time == "5"] <- "4"
metadata$Time2[metadata$Time == "10"] <- "9"
metadata$Time2[metadata$Time == "15"] <- "14"
metadata$Time2[metadata$Time == "20"] <- "19"
metadata$Time2[metadata$Time == "25"] <- "24"
metadata$Time2[metadata$Time == "30"] <- "29"
metadata$Time2[metadata$Time == "35"] <- "34"

## Save metadata file 
readr::write_csv(metadata, here::here("Output Files/01 - Metadata - Output/metadata.csv"))

## Read in dna quant file and combine with metadata
read.csv(here::here("DNA_quant.csv")) -> dna_quant
merge(metadata, dna_quant, by = "Samples") -> meta_quant
readr::write_csv(meta_quant, here::here("Output Files/01 - Metadata - Output/metadata_quant.csv"))

## Will probably want to add the library size column?