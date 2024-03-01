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
length(first_split_kinase)#18


# .bam files to count matrix
count_matrix <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/count_matrix.rds")
dim(count_matrix) #27 samples
table(duplicated(colnames(count_matrix))) # no duplicates
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
table(duplicated(colnames(count_matrix))) # 2 duplicates now
table(colnames(count_matrix)) # SJST031362_D1 and SJST033312_D1

length(labsheet$smpID)#27 samples
table(duplicated(labsheet$smpID)) # 0 duplicates now

setdiff(colnames(count_matrix), labsheet$smpID) #"SJST030375_R2"
setdiff(labsheet$smpID, colnames(count_matrix)) #"STST030375_R2" "SJST030375_D1" "SJST032952_D1"


#@@@ Conclusion: count_matrix has one sample (STST030375_R2) which labsheet doesn't; labsheet has two sample ("SJST030375_D1", "SJST032952_D1") which count_matrix doesn't. Total: Those 3 samples must fall out of our analysis 

count_matrix <- count_matrix[, c(1:9, 11:27), drop = FALSE] #it's a matrix. It requeries it's own subsetting
dim(count_matrix) #26 samples 
labsheet <- labsheet[!labsheet$smpID %in% c("STST030375_R2", "SJST030375_D1", "SJST032952_D1"), ]
dim(labsheet) #24 samples

setdiff(colnames(count_matrix), labsheet$smpID) #nada
setdiff(labsheet$smpID, colnames(count_matrix)) # nada


#@@@ Conclusion: We need to duplicate SJST031362_D1 and SJST033312_D1 in the labsheet
# Solve it in labsheet
row1 <- labsheet[labsheet$smpID %in% 'SJST031362_D1', ]
row1$smpID <- 'SJST031362_D1.1'
row2 <- labsheet[labsheet$smpID %in% 'SJST033312_D1', ]
row2$smpID <- 'SJST033312_D1.1'
labsheet <- rbind(labsheet, row1, row2)
# Solve it in count_matrix
colnames(count_matrix) <- make.names(colnames(count_matrix), unique = TRUE)


length(labsheet$smpID) #26 samples (2 duplicates) but still unique smpIDs
length(colnames(count_matrix))#26 samples (2 duplicates) but still unique smpIDs


save(count_matrix, labsheet, file = '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/reviewed_smpID_countmtx_labsheet.rda')



