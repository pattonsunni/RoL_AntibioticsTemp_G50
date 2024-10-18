# Author: Sunni Patton
# Last edited: 10/18/24
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
colnames(track) <- c("Samples", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "nochloro", "nomito", "noNA.phylum")


### Extract information 
time <- as.character(sapply(strsplit(track$Samples, "_T"), `[`,2))
time <- as.character(sapply(strsplit(time, "_"), `[`, 1)) 

treatment <- as.character(sapply(strsplit(track$Samples, "G50_T\\d+\\_"), `[`, 2)) 
treatment <- as.character(sapply(strsplit(treatment, "_"), `[`, 1))

replicate <- as.character(sapply(strsplit(track$Samples, "G50_T\\d+\\_"), `[`, 2))
replicate <- as.character(sapply(strsplit(replicate, "_"), `[`, 2))

coralID <- as.character(sapply(strsplit(track$Samples, "G50_"), `[`, 2))

tankID <- as.character(sapply(strsplit(coralID, "T\\d+\\_"), `[`, 2))

### Make dataframe
metadata <- data.frame(track$Samples)

metadata$CoralID <- coralID
metadata$TankID <- tankID
metadata$Time <- time
metadata$Treatment <- treatment
metadata$TankReplicate <- replicate

colnames(metadata)[1] <- "Samples"

metadata$Treatment[metadata$Samples == "NC_1" | 
                     metadata$Samples =="NC_2" |
                     metadata$Samples == "NC_3" |
                     metadata$Samples == "NC_4"] <- "Negative Control"

metadata$Treatment_Long <- paste0(metadata$Treatment)

metadata$Treatment_Long[metadata$Treatment == "0"] <- "Blank"
metadata$Treatment_Long[metadata$Treatment == "1"] <- "Antibiotic"
metadata$Treatment_Long[metadata$Treatment == "2"] <- "Temperature"
metadata$Treatment_Long[metadata$Treatment == "3"] <- "Antibiotic + Temperature"

## Add column to specify actual treatment start
metadata$Pretreatment <- paste0(metadata$Treatment_Long)

metadata$Pretreatment[metadata$Time == "0"] <- "No Treatment" # all T0
metadata$Pretreatment[metadata$Treatment == "0"] <- "No Treatment" # all blank
metadata$Pretreatment[metadata$Treatment == "2" & (metadata$Time == "2" | metadata$Time == "5")] <- "No Treatment"
metadata$Pretreatment[metadata$Treatment == "3" & (metadata$Time == "2" | metadata$Time == "5")] <- "Antibiotics Only"
metadata$Pretreatment[metadata$Treatment == "1" & metadata$Time != "0"] <- "Antibiotics Only"    
metadata$Pretreatment[metadata$Time != "0" & metadata$Treatment == "2"] <- "Temperature Only" 
metadata$Pretreatment[metadata$Pretreatment == "Antibiotic + Temperature"] <- "Antibiotics + Temperature" 

# Fix time column (end of antibiotic experiment was 96 hours, i.e., 4 days -- not 5. So 5 should be 4, and all times after that should be one less)
metadata$Time_TotalDays <- paste0(metadata$Time)

metadata$Time_TotalDays[metadata$Time == "5"] <- "4"
metadata$Time_TotalDays[metadata$Time == "10"] <- "9"
metadata$Time_TotalDays[metadata$Time == "15"] <- "14"
metadata$Time_TotalDays[metadata$Time == "20"] <- "19"
metadata$Time_TotalDays[metadata$Time == "25"] <- "24"
metadata$Time_TotalDays[metadata$Time == "30"] <- "29"
metadata$Time_TotalDays[metadata$Time == "35"] <- "34"

metadata$Time_TempDays <- paste0(metadata$Time_TotalDays)

metadata$Time_TempDays[metadata$Time_TotalDays == "0" |
                         metadata$Time_TotalDays == "2" |
                         metadata$Time_TotalDays == "4"] <- "0"
metadata$Time_TempDays[metadata$Time_TotalDays == "9"] <- "5"
metadata$Time_TempDays[metadata$Time_TotalDays == "14"] <- "10"
metadata$Time_TempDays[metadata$Time_TotalDays == "19"] <- "15"
metadata$Time_TempDays[metadata$Time_TotalDays == "24"] <- "20"
metadata$Time_TempDays[metadata$Time_TotalDays == "29"] <- "25"
metadata$Time_TempDays[metadata$Time_TotalDays == "34"] <- "30"

## Read in dna quant file and combine with metadata
read.csv(here::here("DNA_quant.csv")) -> dna_quant
merge(metadata, dna_quant, by = "Samples") -> metadata

## Add library size columns
as.data.frame(c(track[1], track[2], track[10])) -> libSize
colnames(libSize) <- c("Samples", "LibrarySize_PreQC", "LibrarySize_PostQC")
merge(metadata, libSize) -> metadata
readr::write_csv(metadata, here::here("Output Files/01 - Metadata - Output/metadata.csv"))
