# https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts
# package recommended by Felipe
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rsubread")

library(Rsubread)

# ETV6 samples ---------------
# Define the paths to your BAM files and annotation file
annotation_file <- '/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/gencode.v19.chr_patch_hapl_scaff.annotation.gtf'
bam_files <- dir('/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/', patter = '.bam')
# Run featureCounts on the BAM files
setwd("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps")
bam_files <- bam_files[-c(8)] # the SJST032952_D4-S9 is not properly aligned

start <- Sys.time()

counts <- featureCounts(files = bam_files[],
                        annot.ext = annotation_file,
                        isGTFAnnotationFile = TRUE,
                        strandSpecific = "2",
                        isPairedEnd=TRUE,
                        nthreads = 32)
end <- Sys.time()
print(end - start) # 5.331839 mins on nthreads = 32

saveRDS(counts, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_counts.rds')

# Tidying count matrix up
ETV6_counts <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_counts.rds")
head(ETV6_counts$annotation)
head(ETV6_counts$counts)
ETV6_counts <- as.matrix(ETV6_counts$counts)
library(stringr)
first_split <- str_split_fixed(as.character(colnames(ETV6_counts)), "[Aligned]", 2)[,1]
colnames(ETV6_counts) <- first_split
head(ETV6_counts)
ETV6_smp_ID <- colnames(ETV6_counts)
metadata_1 <- data.frame(ETV6_smp_ID, 
           rep('ETV6', length(ETV6_smp_ID)))

names(metadata_1) <- c('smpID', 'group')


# Kinase samples ---------------
# Define the paths to your BAM files and annotation file
annotation_file <- '/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/gencode.v19.chr_patch_hapl_scaff.annotation.gtf'
bam_files <- dir('/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_smps/', patter = '.bam')
# Run featureCounts on the BAM files
setwd("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_smps")

start <- Sys.time()

counts <- featureCounts(files = bam_files[], 
                        annot.ext = annotation_file, 
                        isGTFAnnotationFile = TRUE,
                        strandSpecific = "2",
                        isPairedEnd=TRUE,
                        nthreads = 32)
end <- Sys.time()
print(end - start) #8.907089 mins on nthreads = 32

saveRDS(counts, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_counts.rds')

# Tidying count matrix up
kinase_counts <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_counts.rds")
head(kinase_counts$annotation)
head(kinase_counts$counts)
kinase_counts <- as.matrix(kinase_counts$counts)
library(stringr)
first_split <- str_split_fixed(as.character(colnames(kinase_counts)), "[Aligned]", 2)[,1]
colnames(kinase_counts) <- first_split
head(kinase_counts)
kinase_smp_ID <- colnames(kinase_counts)
metadata_2 <- data.frame(kinase_smp_ID, 
                       rep('kinase', length(kinase_smp_ID)))

names(metadata_2) <- c('smpID', 'group')

# Combine metadata (there's only group info here)
metadata <- plyr::rbind.fill(metadata_1, metadata_2)
rownames(metadata) <- metadata$smpID
metadata$group <- plyr::mapvalues(metadata$group, from = c("ETV6", "kinase"), to = c("ETV6_NTRK3_fused", "Kinase_fused")) #changing group names
saveRDS(metadata, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/metadata.rds')

dim(ETV6_counts) #63568     9
dim(kinase_counts) #63568    18
identical(rownames(ETV6_counts), rownames(kinase_counts)) # TRUE
count_matrix <- cbind(ETV6_counts, kinase_counts)
head(count_matrix)

saveRDS(count_matrix, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/count_matrix.rds')


