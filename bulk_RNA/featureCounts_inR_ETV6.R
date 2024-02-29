gtf_path <- '/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/gencode.v19.chr_patch_hapl_scaff.annotation.gtf'
out_dir <- '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/counts/'
bam_path <- '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/'
bam_file <- dir('/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/', patter = '.bam')

# system(paste0("featureCounts -a -s 2 ", gtf_path, ' -T 20 -o ', out_dir, 'ETV6_counts.txt', ' ', bam_path, bam_file[1]))   
# system(paste0("featureCounts -p -a ", gtf_path, ' -T 20 -o ', out_dir, 'ETV6_counts.txt', ' ', bam_path, '*.bam'))   


# this worked for one .bam
# system(paste0("featureCounts -p -a ", gtf_path, ' -T 20 -o ', out_dir, ,'ETV6_counts.txt', ' ', bam_path, bam_file[1]))   

# library(stringr)
# file_split_1 <- str_split_fixed(as.character(bam_file), "[Aligned]", 2)[,1]
# unique_sample_ids_1 <- unique(file_split_1)
# 
# unique_sample_ids_2 <- c("SJBT032721_D1-S3", "SJIFS031085_D1-S9", "SJST030433_D1-S3", "SJST030567_D1-S10", "SJST032838_D1-S4", "SJST032952_D2-S5",
#   "SJST032952_D4-S4", "SJST032952_D4-S9", "SJST033312_D1-S13", "SJST033312_D1-S2")
# 
# identical(unique_sample_ids_1,unique_sample_ids_2 ) #TRUE


unique_sample_ids <- c("SJBT032721_D1-S3", "SJIFS031085_D1-S9", "SJST030433_D1-S3", "SJST030567_D1-S10", "SJST032838_D1-S4", "SJST032952_D2-S5",
  "SJST032952_D4-S4", "SJST032952_D4-S9", "SJST033312_D1-S13", "SJST033312_D1-S2")

for (i in 1:length(unique_sample_ids)) {
  
system((paste0("featureCounts -p -a ", gtf_path, ' -T 20 -o ', out_dir, unique_sample_ids[i],'ETV6_counts.txt', ' ', bam_path, bam_file[i])))

}

# Check it if we ran all bam samples ------------
# all the bam samples
unique_sample_ids
# catch all .bam files 
library(stringr)
bam_out <- dir("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/", pattern = ".bam$")
file_split_1 <- str_split_fixed(as.character(bam_out), "[Aligned]", 2)[,1]
unique_sample_ids_check <- unique(file_split_1)
# compare them. If it's TRUE so you're all covered. All the samples has been aligned
identical(unique_sample_ids, unique_sample_ids_check) #TRUE 








