### Description: FASTQ alignment and removal of bad/outlier RNA samples
### Technology: bulk RNA-seq

#.ppt presentation at
#/Users/mmarcao/Library/CloudStorage/OneDrive-St.JudeChildren'sResearchHospital/Jasmine_group/LarissaFurtado_DNAmet_RNAseq_integration


### Align ETV6 smps
files <- dir("/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors/", pattern = "fastq")
path = "/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors/"
path.out = "/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/"
# #out.prefix = str_extract(files,"[^_]*_[^_]*")
# library(stringr)
# file_split_1 <- str_split_fixed(as.character(files), "[_][R][0-9][_][00]", 2)[,1]
# file_split_2 <-  str_split_fixed(as.character(file_split_1), "[_][L]", 2)[,1]
# unique_sample_ids <- unique(file_split_2)
# saveRDS(unique_sample_ids, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/unique_sample_ids.rds')
# unique_sample_ids <- readRDS('/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/unique_sample_ids.rds')


unique_sample_ids <- c("SJBT032721_D1-S3", # done
                       "SJIFS031085_D1-S9", # done
                       "SJST030433_D1-S3", # done 
                       "SJST030567_D1-S10", # done
                       "SJST032838_D1-S4", # done 
                       "SJST032952_D1-S9", # NEED TO REVIEW THE FILE PAIR. MULTIPLE FASTQ.
                       "SJST032952_D2-S5", # done
                       "SJST032952_D4-S4", # done
                       "SJST032952_D4-S9", # NEED TO REVIEW  - only one fastq
                       "SJST033312_D1-S13", # done
                       "SJST033312_D1-S2") # done - it's a way small file.BUT it's fine. It matches to the size of their fastq

# library(dplyr)
for(i in 1:length(unique_sample_ids)) {
  #i <- 1
  samples <- files[grep(unique_sample_ids[i], files)]
  samples_R1 <- samples[grep("R1_0", samples)]
  samples_R2 <- samples[grep("R2_0", samples)]
  
  if (length(samples) > 2) {
    # test code
    # print('4 fastqs')
    system(paste0("STAR --genomeDir /media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2 --runThreadN 20 --readFilesIn ",
                  path, samples_R1[1],',',path,samples_R1[2],' ',path,samples_R2[1],',',path,samples_R2[2], " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ",path.out, unique_sample_ids[i]))
    
  } else { 
    # test code
    # print(paste0(length(samples), ' fastqs'))
    system((paste0("STAR --genomeDir /media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2 --runThreadN 20 --readFilesIn ",
                   path, samples_R1[1],' ',path,samples_R2[1], " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ",path.out, unique_sample_ids[i])))
  }
}




# Check it if we ran all fastq samples ------------
# all the fastq samples
unique_sample_ids
# catch all .bam files 
library(stringr)
bam_out <- dir("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/", pattern = ".bam$")
file_split_1 <- str_split_fixed(as.character(bam_out), "[Aligned]", 2)[,1]
unique_sample_ids_check <- unique(file_split_1)
# compare them. If it's TRUE so you're all covered. All the samples has been aligned
identical(unique_sample_ids, unique_sample_ids_check) #FALSE 
# SJST032952_D1-S9 - NEED TO REVIEW THE FILE PAIR
# SJST032952_D4-S9 - NEED TO REVIEW - ONLY ONE FASTQ FILE. It has a .bam output but it's not right


### Align kinase smps
files <- dir("/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_kinase_fused_tumors/", pattern = "fastq")
path = "/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_kinase_fused_tumors/"
path.out = "/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_smps/"
#out.prefix = str_extract(files,"[^_]*_[^_]*")
unique_sample_ids <- c("SJST030375_R2-S1",
                       "SJST031362_D1-S21",
                       "SJST031362_D1-S8",
                       "SJST031690_D1-S1",
                       "SJST031792_D1-S19",
                       "SJST031920_D1-S3",
                       "SJST032767_D1-S1",
                       "SJST032767_D2-S8",
                       "SJST033308_D1-S8",
                       "SJST033389_D1-S3",
                       "SJST033491_D1-S2",
                       "SJST033791_D1-S3",
                       "SJST033835_D1-S10",
                       "SJST033983_D1-S6",
                       "SJST034036_D1-S5",
                       "SJST034397_D1-S6",
                       "SJST034534_D1-S3",
                       "SJST034815_D1-S3")



# library(dplyr)
for(i in 1:length(unique_sample_ids)) {
  
  samples <- files[grep(unique_sample_ids[i], files)]
  samples_R1 <- samples[grep("R1_0", samples)]
  samples_R2 <- samples[grep("R2_0", samples)]
  
  if (length(samples) > 2) {
    system(paste0("STAR --genomeDir /media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2 --runThreadN 20 --readFilesIn ", path, samples_R1[1],',',path,samples_R1[2],' ',path,samples_R2[1],',',path,samples_R2[2], " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ",path.out, unique_sample_ids[i]))
    
  } else { system(paste0("STAR --genomeDir /media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2 --runThreadN 20 --readFilesIn ", path, samples_R1[1],' ', path, samples_R2[1], " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ", path.out, unique_sample_ids[i]))
  }
}



# Check it if we ran all fastq samples ------------
# all the fastq samples
unique_sample_ids
# catch all .bam files 
library(stringr)
bam_out <- dir("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_smps/", pattern = ".bam$")
file_split_1 <- str_split_fixed(as.character(bam_out), "[Aligned]", 2)[,1]
unique_sample_ids_check <- unique(file_split_1)
# compare them. If it's TRUE so you're all covered. All the samples has been aligned
identical(unique_sample_ids, unique_sample_ids_check) #TRUE 


####### Obs: MultiQC screen shoots on ppt slides 


### From .Bam to Count Matrix  
# https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts
# package recommended by Felipe
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rsubread")

library(Rsubread)
library(DESeq2)
library(sva)
library(dplyr)
library(stringr)

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

# saveRDS(counts, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_counts.rds')

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

# saveRDS(counts, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_counts.rds')

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
# saveRDS(metadata, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/metadata.rds')

dim(ETV6_counts) #63568     9
dim(kinase_counts) #63568    18
identical(rownames(ETV6_counts), rownames(kinase_counts)) # TRUE
count_matrix <- cbind(ETV6_counts, kinase_counts)
head(count_matrix)

# saveRDS(count_matrix, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/count_matrix.rds')





### Resolving sample IDs 
library(stringr)
# Matching samples ID from aligment to Larissa's labsheet ------------------
# bam files/featureCounts to count matrix checkout --------------
# ETV6 samples
bam_files_ETV6 <- dir('/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/', patter = '.bam')
library(stringr)
first_split_ETV6 <- str_split_fixed(as.character(bam_files_ETV6), "[Aligned]", 2)[,1]
# bam_files_ETV6 <- str_split_fixed(as.character(first_split_ETV6), "[-]", 2)[,1]
# kinase samples
library(Rsubread)
bam_files_kinase <- dir('/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_smps/', patter = '.bam')
library(stringr)
first_split_kinase <- str_split_fixed(as.character(bam_files_kinase), "[Aligned]", 2)[,1]
# bam_files_kinase <- str_split_fixed(as.character(first_split_kinase), "[-]", 2)[,1]

# length(bam_files_ETV6) #10 samples
# length(bam_files_kinase) #18 samples
length(first_split_ETV6) #10
first_split_ETV6 <- first_split_ETV6[!first_split_ETV6 %in% 'SJST033312_D1-S13']
length(first_split_ETV6) #9
length(first_split_kinase)#18


# .bam files to count matrix
count_matrix <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/count_matrix.rds")
count_matrix <- as.data.frame(count_matrix)
dim(count_matrix) #27 samples
table(duplicated(colnames(count_matrix))) # no duplicates
count_matrix <- count_matrix[, !colnames(count_matrix) %in% 'SJST033312_D1-S13'] # multiqc fall out

# fixing the sample ID
# library(stringr)
# colnames(count_matrix) <- str_split_fixed(as.character(colnames(count_matrix)), "[-]", 2)[,1]
# from featureCounts/count matrix
# table(colnames(count_matrix))
# table(duplicated(colnames(count_matrix))) # 2 duplicates: SJST031362_D1 and SJST033312_D1

# from .bam files
# table(c(bam_files_ETV6,bam_files_kinase)) 
# table(duplicated(c(bam_files_ETV6,bam_files_kinase))) 
# 3 duplicates
# SJST032952_D4 (had only one fastq. So it started aligning but then stopped.), 
# SJST033312_D1 ()
# SJST031362_D1
table(duplicated(c(first_split_ETV6, first_split_kinase))) # no duplicates

setdiff(colnames(count_matrix), c(first_split_ETV6, first_split_kinase)) #nada
setdiff(c(first_split_ETV6, first_split_kinase), colnames(count_matrix)) #SJST032952_D4-S9
# SJST032952_D4-S9 has only one fastq. So it had been aligning but then stopped it. That's why we see this sample ID on the .bam files

#@@@ Conclusion: the count matrix is good. There're all sampleID as it should be. We don't need to change anything there.


# labsheet/metadata checkout -------------
library(readxl)
labsheet <- read_excel("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/metadata_lvf_sub.xlsx")
labsheet <- data.frame(labsheet)
labsheet <- labsheet[!labsheet$fastq_ID %in% "NA", ]
names(labsheet)[names(labsheet) == 'fastq_ID'] <- 'smpID'
length(labsheet$smpID) #27 samples
table(duplicated(labsheet$smpID)) # no duplicates


# Take note: Labsheet is shorten on sample ID. eg: they have "SJST030375_D1" but there's still missing a piece of it (eg: SJST030375_D1_XX). So it might get things a bit confusing in the future analysis. 


# Match count matrix to labsheet -----------
# First, make their sampleID match
library(stringr)
colnames(count_matrix) <- str_split_fixed(as.character(colnames(count_matrix)), "[-]", 2)[,1]
# Comparing both sampleID
length(colnames(count_matrix)) #27 samples
table(duplicated(colnames(count_matrix))) # 1 duplicate
table(colnames(count_matrix)) # SJST031362_D1 
# Solve duplicate (SJST031362_D1) in count matrix 
colnames(count_matrix) <- make.names(colnames(count_matrix), unique = TRUE)
colnames(count_matrix)

length(labsheet$smpID)#27 samples
table(duplicated(labsheet$smpID)) # 0 duplicates now

setdiff(colnames(count_matrix), labsheet$smpID) #"SJST030375_R2" "SJST031362_D1.1"
setdiff(labsheet$smpID, colnames(count_matrix)) #"STST030375_R2" "SJST030375_D1" "SJST032952_D1"


#@@@ Conclusion: STST030375_R2 is misspelled in labshet; labsheet has two sample ("SJST030375_D1", "SJST032952_D1") which count_matrix doesn't. Total: Those 3 samples must fall out of our analysis; we must duplicate SJST031362_D1 in the labsheet

# editing misspelling 
labsheet$smpID[10] <- 'SJST030375_R2'
# deleting samples there're only in labsheet
labsheet <- labsheet[!labsheet$smpID %in% c("SJST030375_D1", "SJST032952_D1"), ]
# duplicate SJST033312_D1 in the labsheet
row1 <- labsheet[labsheet$smpID %in% 'SJST031362_D1', ]
row1$smpID <- 'SJST031362_D1.1'
labsheet <- rbind(labsheet, row1)



length(labsheet$smpID) #26 samples (1 uplicate) but still unique smpIDs
length(colnames(count_matrix))#26 samples (1 duplicate) but still unique smpIDs
length(intersect(labsheet$smpID, colnames(count_matrix))) #25
setdiff(labsheet$smpID, colnames(count_matrix))

# save(count_matrix, labsheet, file = '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/reviewed_smpID_countmtx_labsheet.rda')


