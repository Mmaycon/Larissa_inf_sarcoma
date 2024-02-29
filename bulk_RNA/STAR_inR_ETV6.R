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
  "SJST032952_D1-S9", # NEED TO REVIEW THE FILE PAIR
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







