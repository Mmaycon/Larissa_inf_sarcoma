### Description: from IDAT files to DNAmet matrix (beta-values matrix) using ST and public samples
### Technology: 850k DNAmet array


# .ppt presentation at
# /Users/mmarcao/Library/CloudStorage/OneDrive-St.JudeChildren'sResearchHospital/Jasmine_group/LarissaFurtado_DNAmet_RNAseq_integration

# Source to Differential Region methylation (DMR) at # Source at https://github.com/hamidghaedi/Methylation_Analysis?tab=readme-ov-file



# Process IDAT files into beta-value matrix 
# Load packages
library(BiocManager)
require(IlluminaHumanMethylationEPICmanifest)
require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) 
library(shinyMethyl)
library(bumphunter)
library(minfi)
library(dplyr)

### PROCESSING IDAT FILES
# Load metadata
load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/met_meta.rda")
met # Felipe's DNAmet matrix
head(meta) # Metadata/labsheet
dim(meta)
# Creating group information
met$Group <- c(rep('ETV6_NTRK3_fused', 22), rep('Kinase_fused', 16))
target <- data.frame(sample = colnames(met), group = met$Group, source = 'StJude')
target$source[c(10:22)] <- 'PublicData'
met <- NULL # Delete it. I'm processing it on my own

# 1. Load IDAT files -----------
# ETV6_NTRK3 samples
ETV6_NTRK3_idat_dir <- '/media/ResearchHome/plummgrp/home/common/LFurtado-colab/idat_files_ETV6-NTRK3_fused_tumors'
RGC_data_etv6 <- read.metharray.exp(base = ETV6_NTRK3_idat_dir, targets = NULL, force = TRUE, recursive = T) 

summary.idat <- shinySummarize(RGC_data_etv6)
runShinyMethyl(summary.idat) 

# kinase samples
kinase_idat_dir <- '/media/ResearchHome/plummgrp/home/common/LFurtado-colab/idat_files_kinase_fused_tumors'
RGC_data_kinase <- read.metharray.exp(base = kinase_idat_dir, targets = NULL, force = TRUE, recursive = T) 

summary.idat <- shinySummarize(RGC_data_kinase)
runShinyMethyl(summary.idat) 

# GEO samples (GSE140686)
geo_idat_dir <- '/media/ResearchHome/plummgrp/home/common/LFurtado-colab/idat_files_GSE140686'
RGC_data_geo <- read.metharray.exp(base = geo_idat_dir, targets = NULL, force = TRUE, recursive = T) 

summary.idat <- shinySummarize(RGC_data_geo)
runShinyMethyl(summary.idat) 
summary.idat@sampleNames

# all samples 
all_idat_dir <- '/media/ResearchHome/plummgrp/home/common/LFurtado-colab'
RGC_data_all <- read.metharray.exp(base = all_idat_dir, targets = NULL, force = TRUE, recursive = T) 

summary.idat <- shinySummarize(RGC_data_all)
runShinyMethyl(summary.idat) 


# 2. pvalue detection
detP <- detectionP(RGC_data_all, type = "m+u") 
table(detP > 0.05)


# 3. Preprocess the data
proc_data <- preprocessRaw(RGC_data_all) # that was the best processing option I got


# 4. Mask probes that failed p-value detection
proc_data_r <- ratioConvert(proc_data)
is.na(assays(proc_data_r)$Beta) <- (detP[rownames(proc_data_r), colnames(proc_data_r)] > 0.05)
beta <- getBeta(proc_data_r)
head(beta)
dim(beta) #865859     38

# 5. Remove mask probes
load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/hm450.anno.Rda") # this is 450k but later I'm using EPIC annotation so the non-wanted CpG will be taken off. It Should not be a concern here.
probes_remove <- subset(hm450.anno, chrm_A %in% c("chrX","chrY", "chrM") & MASK_general == FALSE)$probeID 
beta <- beta[!rownames(beta) %in% probes_remove, ]
dim(beta) #856801     38
head(beta)
beta_meta <- target
# save(beta, beta_meta, file = '/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/beta_and_meta_data.rda')


# QC PLOTS ------------------
qc <- getQC(proc_data)

plotQC(qc) # it gets good with the preprocessRaw() method

densityPlot(proc_data, sampGroups = target$group)
densityBeanPlot(proc_data, sampGroups = target$group)
controlStripPlot(RGC_data_all, controls="BISULFITE CONVERSION II")


