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




