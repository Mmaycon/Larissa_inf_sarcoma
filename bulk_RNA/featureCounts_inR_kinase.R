gtf_path <- '/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/gencode.v19.chr_patch_hapl_scaff.annotation.gtf'
out_dir <- '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_smps/counts/'
bam_path <- '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_smps/'
bam_file <- dir('/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_smps/', patter = '.bam')


# Getting sample ID
# library(stringr)
# file_split_1 <- str_split_fixed(as.character(bam_file), "[Aligned]", 2)[,1]
# unique_sample_ids <- unique(file_split_1)
# 
# unique_sample_ids_2 <- c("SJST030375_R2-S1", "SJST031362_D1-S21", "SJST031362_D1-S8", "SJST031690_D1-S1", "SJST031792_D1-S19",
#                 "SJST031920_D1-S3", "SJST032767_D1-S1", "SJST032767_D2-S8", "SJST033308_D1-S8", "SJST033389_D1-S3",
#                 "SJST033491_D1-S2", "SJST033791_D1-S3", "SJST033835_D1-S10", "SJST033983_D1-S6", "SJST034036_D1-S5",
#                 "SJST034397_D1-S6", "SJST034534_D1-S3", "SJST034815_D1-S3")
# identical(unique_sample_ids, unique_sample_ids_2) # TRUE

unique_sample_ids <- c("SJST030375_R2-S1", "SJST031362_D1-S21", "SJST031362_D1-S8", "SJST031690_D1-S1", "SJST031792_D1-S19",
                                         "SJST031920_D1-S3", "SJST032767_D1-S1", "SJST032767_D2-S8", "SJST033308_D1-S8", "SJST033389_D1-S3",
                                         "SJST033491_D1-S2", "SJST033791_D1-S3", "SJST033835_D1-S10", "SJST033983_D1-S6", "SJST034036_D1-S5",
                                         "SJST034397_D1-S6", "SJST034534_D1-S3", "SJST034815_D1-S3")


for (i in 1:length(unique_sample_ids)) {
  
  system((paste0("featureCounts -p -a ", gtf_path, ' -T 20 -o ', out_dir, unique_sample_ids[i],'kinase_counts.txt', ' ', bam_path, bam_file[i])))
  
}



# # Check it if we ran all bam samples ------------
# # all the bam samples
# unique_sample_ids
# # catch all .bam files 
# library(stringr)
# bam_out <- dir("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_smps/", pattern = ".bam$")
# file_split_1 <- str_split_fixed(as.character(bam_out), "[Aligned]", 2)[,1]
# unique_sample_ids_check <- unique(file_split_1)
# # compare them. If it's TRUE so you're all covered. All the samples has been aligned
# identical(unique_sample_ids, unique_sample_ids_check) #TRUE 

