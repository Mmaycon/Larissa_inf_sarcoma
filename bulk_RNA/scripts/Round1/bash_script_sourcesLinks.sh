#------------- Larissa bulk RNA-seq alignment -----------#
# tutorials  ------------
    # RNAseq processing from Sanbomics group at 
    # https://www.youtube.com/watch?v=oRC406tbB8w and
    # https://sanbomics.com/2022/01/08/complete-rnaseq-alignment-guide-from-fastq-to-count-table/
    
# directories  ------------
    # fasta files are living at /media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA
    # we should save temporarily at /mnt/scratch1/maycon/LFurtado_colab_RNA_processed

# activate conda commands  ------------
    # > nano ~/.bashrc
    # add to the last line the string below
    # export PATH="/mnt/scratch1/miniconda3/bin:$PATH"
    # exit and save it
    # > source ~/.bashrc
    # reopen the terminal
    # now we can see (base) before your usarname
    # check conda installation 
    # > conda --version
    

# install STAR - conda isn't working  ------------
    # I need permission to install it
    # > conda install -c bioconda star 
    # or
    # from https://github.com/alexdobin/STAR

# Download genome fasta GTF annotation (GRCh37 version)  ------------
    # from https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens
    # dna fasta - pick up "primary assembly"
    # GTF annoation - https://www.gencodegenes.org/human/release_19.html
    # I donwload it to my local machine and then I rsync them to the server - wget was not working that time
    # > rsync -avP ./gencode.v19.annotation.gtf.gz mmarcao@plummerlab01:/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens
    # > rsync -avP ./Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.gz mmarcao@plummerlab01:/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens
    
# Generating genome index ------------
    # at /media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens
    # > gunzip *.gz
    # > mkdir HS_STAR_dir # a directory for the index
    # run STAR at /GRCh37_homosapiens dir
cd /media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens
STAR --runMode genomeGenerate \
#--genomeDir HS_STAR_dir/ \
--genomeFastaFiles Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa \
#--sjdbGTFfile gencode.v19.annotation.gtf \
--sjdbGTFfile gencode.v19.chr_patch_hapl_scaff.annotation.gtf \
--runThreadN 8

# try more thereads
cd /media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens
STAR --runMode genomeGenerate \
--genomeDir HS_STAR_dir_2/ \
--genomeFastaFiles Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa \
--sjdbGTFfile gencode.v19.chr_patch_hapl_scaff.annotation.gtf \
--runThreadN 44



# command line tips  ------------
    # donwload from the internet 
    # > wget -P <preference dir>/ <url address>
    

# Align your bulk RNA-seq data ------------
# data Obs
    # we have samples with 2 files and samples with 4 files
    # it means there're samples that got sequenced in two lanes (2 files each)
    # we shall merge lanes within the same sample AFTER mapping them to the genome (aka after alignment)
    
    
# Code test
THREADS=44
GENOME_DIR="/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2" 
FASTQ_DIR="/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors"
read1="${FASTQ_DIR}/SJBT032721_D1-S3_L001_R1_001.fastq.gz"
read2="${FASTQ_DIR}/SJBT032721_D1-S3_L001_R2_001.fastq.gz"
#OUTPUT_DIR='/mnt/scratch1/maycon/LFurtado_colab_RNA_processed/ETV6_NTRK3_fused_tumors/code_test'
OUTPUT_DIR='/media/ResearchHome/plummgrp/home/common/LFurtado-colab/STAR_out/ETV6_NTRK3_fused_tumors/code_test'

# Run STAR for the single sample
STAR --runThreadN ${THREADS} \
     --genomeDir ${GENOME_DIR} \
     --readFilesIn ${read1} ${read2} \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --outFileNamePrefix ${OUTPUT_DIR}/sample_name_test


# Run qualimap on it
BAM_FILE='/media/ResearchHome/plummgrp/home/common/LFurtado-colab/STAR_out/ETV6_NTRK3_fused_tumors/code_test/sample_name_testAligned.sortedByCoord.out.bam'
OUTPUT_DIR='/media/ResearchHome/plummgrp/home/common/LFurtado-colab/STAR_out/ETV6_NTRK3_fused_tumors/code_test/map_report'
qualimap bamqc -bam ${BAM_FILE} -outdir ${OUTPUT_DIR} -outformat pdf





#!/bin/bash

# Path to the STAR genome index
GENOME_DIR="/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir" 

# Directory where the FASTQ files are located
FASTQ_DIR="/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors"

# Directory where you want to store the output files
OUTPUT_DIR="/mnt/scratch1/maycon/LFurtado_colab_RNA_processed/ETV6_NTRK3_fused_tumors"

# Number of threads to use for STAR
THREADS=8

# Loop through each sample directory
for sample_dir in ${FASTQ_DIR}/*; do
    if [ -d "${sample_dir}" ]; then
        # Extract the sample name from the directory path
        sample_name=$(basename "${sample_dir}")

        # Define the input read files (assuming two files per sample, one for each end of the pair)
        read1="${sample_dir}/${sample_name}_L*_R1_001.fastq.gz"
        read2="${sample_dir}/${sample_name}_L*_R2_001.fastq.gz"

        # Run STAR for each sample
        STAR --runThreadN ${THREADS} \
             --genomeDir ${GENOME_DIR} \
             --readFilesIn ${read1} ${read2} \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts \
             --outFileNamePrefix ${OUTPUT_DIR}/${sample_name}_
    fi
done