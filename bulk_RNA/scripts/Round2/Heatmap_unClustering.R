# ------------ Plot heatmap --------------- #
# Add "fusion” e “histology” labels. Remove “sample type” from the visualization

# Load packages 
library(dplyr)
library(DESeq2)

# Load expression matrix 
load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/bulk_RNA/objects/Round1/reviewed_smpID_countmtx_labsheet.rda")
labsheet
labsheet$group<- gsub(' ', '_', labsheet$group)

# Convert gene ID
rownames(count_matrix) <- str_split_fixed(as.character(rownames(count_matrix)), "[.]", 2)[,1] # good it doesn't have duplicates ENSEMBL here
count_matrix <- data.frame(count_matrix)
ens <- rownames(count_matrix) #genes from DE
library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')

symbols <- symbols[!is.na(symbols)]

symbols <- data.frame(ENSEMBL = names(symbols),
                      SYMBOL = as.vector(symbols))

count_matrix$ENSEMBL <- rownames(count_matrix)

count_matrix <- merge(count_matrix, symbols, by='ENSEMBL')
count_matrix <- count_matrix[!duplicated(count_matrix$SYMBOL), ]
rownames(count_matrix) <- count_matrix$SYMBOL
count_matrix$ENSEMBL <- NULL
count_matrix$SYMBOL <- NULL
head(count_matrix); dim(count_matrix) #33909    26

# remove outlier
count_matrix <- count_matrix[, !colnames(count_matrix) %in% "SJST034534_D1"]
dim(count_matrix) #25 smp
labsheet <- labsheet[!labsheet$smpID %in% "SJST034534_D1", ]
dim(labsheet)#25 smp
# make it identical (same sequence)
count_matrix <- count_matrix[, labsheet$smpID]
identical(colnames(count_matrix), labsheet$smpID) #TRUE
# add rownames into labsheet
rownames(labsheet) <- labsheet$smpID

# Add fusion variable from the labsheet (excel file) into our metadata (the "labsheet" object)
# from the current labsheet shared via One Drive:
smpID <- c("SJST033983_D1", "SJST032767_D1", "SJST032767_D2", "SJST034397_D1", "SJST031920_D1",
             "SJST033491_D1", "SJST033389_D1", "SJST031690_D1", "SJST031362_D1", "STST030375_R2",
             "SJST030375_D1", "SJST031792_D1", "SJST033791_D1", "SJST033835_D1", "SJST033308_D1",
             "SJST034534_D1", "SJST034036_D1", NA, NA, "SJST034815_D1", "SJST032952_D2",
             "SJST032952_D1", "SJST032952_D4", "SJST033312_D1", "SJST032838_D1", "SJBT032721_D1",
             "SJIFS031085_D1", "SJST030567_D1", "SJST030433_D1", NA)

fusion_type <- c("NRF1::BRAF", "PLEKHH2::ALK", "PLEKHH2::ALK", "TPR::NTRK1", "TPM3::NTRK1",
             "TPM3::NTRK1", "TPR::NTRK1", "LMNA::NTRK1", "LMNA::NTRK1", "PDE4DIP::NTRK1",
             "PDE4DIP::NTRK1", "EML4::NTRK3", "RBPMS::NTRK3", "RCC1::ALK", 
             "EML4::NTRK3","TTYH3::BRAF","TPM3::NTRK1","EML4::ALK Outside lab",
             'TPM3::NTRK1 Outside lab',"TPM3::NTRK1","ETV6:NTRK3","ETV6:NTRK3","ETV6:NTRK3",
             'ETV6::NTRK3','ETV6::NTRK3','ETV6::NTRK3','ETV6::NTRK3','ETV6::NTRK3',
             'ETV6::NTRK3','ETV6::NTRK3')

# Create a dataframe
fusion_df <- data.frame(smpID = smpID, fusion_type = fusion_type)
fusion_df <- fusion_df[!fusion_df$smpID %in% NA, ]
fusion_df$smpID[fusion_df$smpID %in% 'STST030375_R2'] <- 'SJST030375_R2'

table(duplicated(fusion_df$smpID))
setdiff(rownames(labsheet), fusion_df$smpID) 
length(intersect(rownames(labsheet), fusion_df$smpID)) #24 #SJST031362_D1.1 we have its annotation for SJST031362_D1. So it's okay move on with these 24 sample ID 

# Merge fusion_type variable into labsheet
labsheet_hm <- merge(labsheet, fusion_df, by = 'smpID')


# Load normalized expression matrix 
dds <- readRDS('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/bulk_RNA/objects/Round2/DESeq2_output.rds')
normalized_counts <- counts(dds, normalized=TRUE)
# Filter normalized_counts by labsheet
normalized_counts_hm <- normalized_counts[, colnames(normalized_counts) %in% labsheet_hm$smpID]

### PLOT - heatmap - whole RNAseq
library(pheatmap)
library(RColorBrewer)
normalized_counts_hm # exp matrix
labsheet_hm # metadata

my_sample_col <- data.frame(row.names = labsheet_hm$smpID, 
                            Group = labsheet_hm$group, 
                            Histology = labsheet_hm$histologic_diagnosis,
                            Fusion = labsheet_hm$fusion_type) # patient ID in row names; any columns are for sample information

pheatmap(normalized_counts_hm, 
         #annotation_row = my_probe_col, 
         annotation_col = my_sample_col,
         show_rownames = FALSE,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdYlBu")))(100),
         main = 'Whole RNAseq - normalized (18,777 features)')


### PLOT - heatmap - whole RNAseq - focus on NTRK label
NTRK <- my_sample_col[grep('NTRK', my_sample_col$Fusion), ]
NTRK_butETV6 <- rownames(NTRK3[!NTRK3$Fusion %in% c('ETV6::NTRK3', 'ETV6:NTRK3'), ])
NTRK_onlyETV6 <- rownames(NTRK3[!rownames(NTRK3) %in% NTRK3_butETV6, ])

my_sample_col$Fusion_concat <- NA
my_sample_col[rownames(my_sample_col) %in% NTRK_butETV6, ]$Fusion_concat <- 'Others_NTRK'
my_sample_col[rownames(my_sample_col) %in% NTRK_onlyETV6, ]$Fusion_concat <- 'ETV6_NTRK3'
my_sample_col[my_sample_col$Fusion_concat %in% NA, ]$Fusion_concat <- 'Not_NTRK'


ann_colors = list(
  Fusion_concat = c(ETV6_NTRK3 = 'black', Not_NTRK = 'blue', Others_NTRK = 'red'))


pheatmap(normalized_counts_hm, 
         #annotation_row = my_probe_col, 
         annotation_col = my_sample_col[, c(1,3,4)],
         show_rownames = FALSE,
         scale = "row",
         annotation_colors = ann_colors,
         main = 'Whole RNAseq - normalized (18,777 features)')



### PLOT - heatmap - only DEGs (ETV6 vs kinase) - focus on NTRK label
DEGs <- readRDS('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/bulk_RNA/objects/Round2/Volcano_DEGs_plot.rds')
head(DEGs)
dim(DEGs)

pheatmap(normalized_counts_hm[rownames(normalized_counts_hm) %in% rownames(DEGs), ], 
         #annotation_row = my_probe_col, 
         annotation_col = my_sample_col[, c(2,3,4)],
         show_rownames = FALSE,
         scale = "row",
         annotation_colors = ann_colors,
         main = 'DEGs - normalized (611 features)')

table(my_sample_col$Histology)

