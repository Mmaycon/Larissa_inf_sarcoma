# #!/bin/bash
# ERROR - I gotta change the way I'm inputing the fastq files

# #!/bin/bash
# THREADS=44
# GENOME_DIR="/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2" 
# FASTQ_DIR="/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors"
# read1="${FASTQ_DIR}/SJBT032721_D1-S3_L001_R1_001.fastq.gz ${FASTQ_DIR}/SJBT032721_D1-S3_L002_R1_001.fastq.gz"
# read2="${FASTQ_DIR}/SJBT032721_D1-S3_L001_R2_001.fastq.gz ${FASTQ_DIR}/SJBT032721_D1-S3_L002_R2_001.fastq.gz"
# #OUTPUT_DIR='/mnt/scratch1/maycon/LFurtado_colab_RNA_processed/ETV6_NTRK3_fused_tumors/code_test'
# OUTPUT_DIR='/media/ResearchHome/plummgrp/home/common/LFurtado-colab/STAR_out/ETV6_NTRK3_fused_tumors/code_test'
# sample_name=$(basename "${FASTQ_DIR}")
# 
# # Run STAR for the single sample
# STAR --runThreadN ${THREADS} \
#      --genomeDir ${GENOME_DIR} \
#      --readFilesIn ${read1} ${read2} \
#      --readFilesCommand zcat \
#      --outSAMtype BAM SortedByCoordinate \
#      #--quantMode GeneCounts \
#      --outFileNamePrefix ${OUTPUT_DIR}/${sample_name}
#      
     
     
     
     
     
     
     
    
# To merge files 
# https://knowledge.illumina.com/software/cloud-software/software-cloud-software-reference_material-list/000002035


# # Worked at my local machine 
# #!/bin/bash
# THREADS=8
# GENOME_DIR="/Users/mmarcao/Documents/RNA_process/HS_STAR_dir_2" 
# FASTQ_DIR="/Users/mmarcao/Documents/RNA_process/FASTAQ"
# read1=("${FASTQ_DIR}/SJBT032721_D1-S3_L001_R1_001.fastq" "${FASTQ_DIR}/SJBT032721_D1-S3_L002_R1_001.fastq")
# read2=("${FASTQ_DIR}/SJBT032721_D1-S3_L001_R2_001.fastq" "${FASTQ_DIR}/SJBT032721_D1-S3_L002_R2_001.fastq")
# OUTPUT_DIR='/Users/mmarcao/Documents/RNA_process/FASTAQ/star_output'
# sample_name=$(basename "${FASTQ_DIR}")
# 
# # Merge the input files
# cat "${read1[@]}" > "${FASTQ_DIR}/merged_R1.fastq"
# cat "${read2[@]}" > "${FASTQ_DIR}/merged_R2.fastq"
# 
# # Run STAR for the single sample
# STAR --runThreadN ${THREADS} \
#      --genomeDir ${GENOME_DIR} \
#      --readFilesIn "${FASTQ_DIR}/merged_R1.fastq" "${FASTQ_DIR}/merged_R2.fastq" \
#      --outSAMtype BAM SortedByCoordinate \
#      --outFileNamePrefix ${OUTPUT_DIR}/${sample_name}



# # Translating it to make it run on server (1 lane)
# THREADS=20
# GENOME_DIR="/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2" 
# FASTQ_DIR="/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors"
# read1="${FASTQ_DIR}/SJBT032721_D1-S3_L001_R1_001.fastq.gz"
# read2="${FASTQ_DIR}/SJBT032721_D1-S3_L001_R2_001.fastq.gz" 
# OUTPUT_DIR='/media/ResearchHome/plummgrp/home/common/LFurtado-colab/STAR_out/ETV6_NTRK3_fused_tumors/code_test'
# sample_name=$(basename "${OUTPUT_DIR}")
# 
# # Run STAR for the single sample
# STAR --runThreadN ${THREADS} \
#      --genomeDir ${GENOME_DIR} \
#      --readFilesCommand zcat \
#      --readFilesIn ${read1} ${read2} \
#      --outSAMtype BAM SortedByCoordinate \
#      --outFileNamePrefix ${OUTPUT_DIR}/${sample_name}


# Try it in one line pnly
# STAR --runThreadN 20 --genomeDir '/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2' --readFilesCommand zcat --readFilesIn '/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors/SJBT032721_D1-S3_L001_R1_001.fastq.gz' '/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors/SJBT032721_D1-S3_L001_R2_001.fastq.gz' --outSAMtype BAM SortedByCoordinate --outFileNamePrefix '/media/ResearchHome/plummgrp/home/common/LFurtado-colab/STAR_out/ETV6_NTRK3_fused_tumors/code_test/SJBT032721_D1-S3'




# # Translating it to make it run on server (2 Lanes)
# THREADS=20
# GENOME_DIR="/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2"
# FASTQ_DIR="/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors"
# read1="${FASTQ_DIR}/SJBT032721_D1-S3_L001_R1_001.fastq.gz","${FASTQ_DIR}/SJBT032721_D1-S3_L002_R1_001.fastq.gz"
# read2="${FASTQ_DIR}/SJBT032721_D1-S3_L001_R2_001.fastq.gz" "${FASTQ_DIR}/SJBT032721_D1-S3_L002_R2_001.fastq.gz"
# OUTPUT_DIR='/mnt/scratch1/maycon/LFurtado_colab_RNA_processed/ETV6_NTRK3_fused_tumors/code_test'
# MERGED_FQ_DRIR='/mnt/scratch1/maycon/LFurtado_colab_RNA_processed/ETV6_NTRK3_fused_tumors/code_test'
# sample_name=$(basename "${OUTPUT_DIR}")
# 
# # Merge the input files
# cat "${read1[@]}" > "${MERGED_FQ_DRIR}/merged_R1.fastq.gz"
# cat "${read2[@]}" > "${MERGED_FQ_DRIR}/merged_R2.fastq.gz"
# 
# # Run STAR for the single sample
# STAR --runThreadN ${THREADS} \
#      --genomeDir ${GENOME_DIR} \
#      --readFilesCommand zcat \
#      --readFilesIn "${MERGED_FQ_DRIR}/merged_R1.fastq.gz" "${MERGED_FQ_DRIR}/merged_R2.fastq.gz" \
#      --outSAMtype BAM SortedByCoordinate \
#      --outFileNamePrefix ${OUTPUT_DIR}/${sample_name}




# It works too (for 2 lane samples!)
THREADS=20
GENOME_DIR="/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2"
FASTQ_DIR="/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors"
read1_1="${FASTQ_DIR}/SJBT032721_D1-S3_L001_R1_001.fastq.gz"
read1_2="${FASTQ_DIR}/SJBT032721_D1-S3_L002_R1_001.fastq.gz"
read2_1="${FASTQ_DIR}/SJBT032721_D1-S3_L001_R2_001.fastq.gz"
read2_2="${FASTQ_DIR}/SJBT032721_D1-S3_L002_R2_001.fastq.gz"
OUTPUT_DIR='/mnt/scratch1/maycon/LFurtado_colab_RNA_processed/ETV6_NTRK3_fused_tumors/code_test'
sample_name=$(basename "${OUTPUT_DIR}")


# Run STAR for the single sample
STAR --runThreadN ${THREADS} \
     --genomeDir ${GENOME_DIR} \
     --readFilesCommand zcat \
     --readFilesIn ${read1_1},${read1_2} ${read2_1},${read2_2} \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix ${OUTPUT_DIR}/${sample_name}



# It worked !!!!
# STAR --runThreadN 20 --genomeDir '/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2' --readFilesCommand zcat --readFilesIn '/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors/SJBT032721_D1-S3_L001_R1_001.fastq.gz' '/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors/SJBT032721_D1-S3_L001_R2_001.fastq.gz' --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /mnt/scratch1/maycon/LFurtado_colab_RNA_processed/ETV6_NTRK3_fused_tumors/code_test




# THREADS=32
# GENOME_DIR="/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2"
# FASTQ_DIR="/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors"
# OUTPUT_DIR='/mnt/scratch1/maycon/LFurtado_colab_RNA_processed/ETV6_NTRK3_fused_tumors'
# sample_name=$(basename "${read1_1}" | cut -d'_' -f1)
# 
# # Check if files from the second lane exist
# if [ -f "${read1_1}" ] && [ -f "${read2_1}" ]; then
#     # Run STAR for samples from 2 lanes
#     STAR --runThreadN ${THREADS} \
#          --genomeDir ${GENOME_DIR} \
#          --readFilesCommand zcat \
#          --readFilesIn ${read1_1},${read1_2} ${read2_1},${read2_2} \
#          --outSAMtype BAM SortedByCoordinate \
#          --outFileNamePrefix ${OUTPUT_DIR}/${sample_name}
# else
#     # Run STAR for samples from 1 lane
#     STAR --runThreadN ${THREADS} \
#          --genomeDir ${GENOME_DIR} \
#          --readFilesCommand zcat \
#          --readFilesIn ${read1_1} ${read2_1} \
#          --outSAMtype BAM SortedByCoordinate \
#          --outFileNamePrefix ${OUTPUT_DIR}/${sample_name}
# fi




THREADS=32
GENOME_DIR="/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2"
FASTQ_DIR="/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors"
OUTPUT_DIR='/mnt/scratch1/maycon/LFurtado_colab_RNA_processed/ETV6_NTRK3_fused_tumors'

read1_1=${FASTQ_DIR}/'SJST030567_D1-S10_L004_R1_001.fastq.gz'
read1_2=${FASTQ_DIR}/'SJST030567_D1-S10_L004_R1_001.fastq.gz'
read2_1=${FASTQ_DIR}/'SJST030567_D1-S10_L003_R2_001.fastq.gz'
read2_2=${FASTQ_DIR}/'SJST030567_D1-S10_L004_R2_001.fastq.gz'

for read1_1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
    read1_2=${read1_1/_R1/_R2}
    read2_1=${read1_1}
    read2_2=${read1_2}

    sample_name=$(basename "${read1_1}" | cut -d'_' -f1)

    if [ -f "${read1_2}" ] && [ -f "${read2_2}" ]; then
        STAR --runThreadN ${THREADS} \
             --genomeDir ${GENOME_DIR} \
             --readFilesCommand zcat \
             --readFilesIn ${read1_1},${read1_2} ${read2_1},${read2_2} \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix ${OUTPUT_DIR}/${sample_name}
    else
        STAR --runThreadN ${THREADS} \
             --genomeDir ${GENOME_DIR} \
             --readFilesCommand zcat \
             --readFilesIn ${read1_1} ${read2_1} \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix ${OUTPUT_DIR}/${sample_name}
    fi
# done




  