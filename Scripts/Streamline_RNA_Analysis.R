library(DESeq2)
library(sva)
library(dplyr)
library(stringr)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(ggplot2)
library(UpSetR)

### Differential Expressed Genes -----------
load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/bulk_RNA/objects/Round1/reviewed_smpID_countmtx_labsheet.rda")
labsheet
count_matrix
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

# Create histologic_for_outlier column 
labsheet$histologic_for_outlier <- 'Others'
outlier_smp <- c('SJST033491_D1',
                 'SJST032952_D2',
                 'SJST032952_D4',
                 'SJST032767_D2',
                 'SJST033312_D1',
                 'SJST033835_D1',
                 'SJST031920_D1')
his_diag <- labsheet[labsheet$smpID %in% outlier_smp, ]$histologic_diagnosis

labsheet[labsheet$smpID %in% outlier_smp[1], ]$histologic_for_outlier <- his_diag[1]
labsheet[labsheet$smpID %in% outlier_smp[2], ]$histologic_for_outlier <- his_diag[2]
labsheet[labsheet$smpID %in% outlier_smp[3], ]$histologic_for_outlier <- his_diag[3]
labsheet[labsheet$smpID %in% outlier_smp[4], ]$histologic_for_outlier <- his_diag[4]
labsheet[labsheet$smpID %in% outlier_smp[5], ]$histologic_for_outlier <- his_diag[5]
labsheet[labsheet$smpID %in% outlier_smp[6], ]$histologic_for_outlier <- his_diag[6]
labsheet[labsheet$smpID %in% outlier_smp[7], ]$histologic_for_outlier <- his_diag[7]


# make it identical (same sequence)
count_matrix <- count_matrix[, labsheet$smpID]
identical(colnames(count_matrix), labsheet$smpID) #TRUE
# add rownames into labsheet
rownames(labsheet) <- labsheet$smpID


# Differential Expression -----
# make it identical (same sequence)
count_matrix <- count_matrix[, labsheet$smpID]
identical(colnames(count_matrix), labsheet$smpID) #TRUE
# add rownames into labsheet
rownames(labsheet) <- labsheet$smpID
# leave labsheet with only variables we want in the matrix design 
labsheet <- labsheet[, c("group", "smp_type") ]

# write.csv(labsheet, "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Objects/metadata_RNA_sample_ID.csv")

head(labsheet)
# Create DESEq object 
dds <- DESeqDataSetFromMatrix(count_matrix, colData = labsheet, design = ~ smp_type + group)

# Important step for further normalization during DESEq()
dds <- estimateSizeFactors(dds)
# sizeFactors(dds) #Idk what it does ... neither if we need to run it
# Pre-filtering 
smallestGroupSize <- min(as.vector(table(dds$group)))
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
# indicate direction for comparison
# dds$condition <- relevel(dds$group, ref = "kinase-fused_tumor")
# dds$group <- relevel(dds$group, ref = "ETV6-NTRK3_fused_tumor")
# dds$group <- relevel(dds$group, ref = "kinase-fused_tumor")

# Run Differential Expression analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", 
                                 "kinase-fused_tumor", #direction
                                 "ETV6-NTRK3_fused_tumor"))
res_kinaseDirection <- as.data.frame(res)

# Filter results for significant genes
# Subset stat. significant genes
res_pos <- subset(res_kinaseDirection, res_kinaseDirection$padj < 0.05 & res_kinaseDirection$log2FoldChange >= 0.5)
res_neg <- subset(res_kinaseDirection, res_kinaseDirection$padj < 0.05 & res_kinaseDirection$log2FoldChange <= -0.5)
res_sig <- rbind(res_pos, res_neg)
dim(res_sig) #617   6
# PLOT - FoldChange distribution - histogram
hist(res_sig$log2FoldChange, 
     main = "FoldChange dist. (s.sig DEGs)", 
     xlab = "FoldChange",
     ylab = "N of genes")

# saveRDS(res_kinaseDirection, "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Objects/DEGs_kinase_direction_noCutoffs.rds")
# saveRDS(res_sig, "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Objects/DEGs_kinase_direction_filtered_617genes.rds")
# saveRDS(dds, "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Objects/dds_RNA.rds")


## PCA 
dds <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Objects/dds_RNA.rds")
normalized_counts <- counts(dds, normalized=TRUE)
pca <- prcomp(t(normalized_counts)) 
aux <- as.data.frame(pca$x[, 1:3]) 
labsheet$smpID <- rownames(labsheet)
scores <- merge(labsheet, aux, by.y=0, by.x="smpID", all.x=T)
# Plot it 
library(ggplot2); theme_set(theme_classic())
ggplot(scores, aes(x=PC1, y=PC2, colour=factor(group), shape = smp_type)) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'), name="Group") +
  xlab(paste0("PC1 (", prettyNum(summary(pca)$importance[2,1]*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  #geom_text_repel(aes(label = smpID)) +
  ggtitle("Normalized / batch corrected (mtx design) - whole transcriptome") 


## Volcano plot 
res_pos <- subset(res_kinaseDirection, res_kinaseDirection$padj < 0.05 & res_kinaseDirection$log2FoldChange >= 0.5)
res_neg <- subset(res_kinaseDirection, res_kinaseDirection$padj < 0.05 & res_kinaseDirection$log2FoldChange <= -0.5)
res_sig <- rbind(res_pos, res_neg)

dim(res_pos) # 451 genes highly expressed in Kinase tumors
dim(res_neg) # 166 genes highly expressed in ETV6-NTRK3 tumors
dim(res_sig) # total of 617 genes differentially expressed 

top_Kinase <- rownames(res_pos[order(res_pos$log2FoldChange, decreasing = TRUE), ])[1:20]
top_ETV6 <- rownames(res_neg[order(res_neg$log2FoldChange, decreasing = FALSE), ])[1:20]
selected_genes <- c(top_Kinase, top_ETV6)

library(EnhancedVolcano) # ISSUE: somehow I can't add gene names onto volcano plot using EnhancedVolcano
res_kinaseDirection <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Objects/DEGs_kinase_direction_noCutoffs.rds")
#rownames(res_sig) <- res_sig$SYMBOL # don't run it if it doesn't need to
EnhancedVolcano(#res_sig,
  res_kinaseDirection, # DEG output not filtered yet
  #lab = rownames(res_kinaseDirection),
  #lab = as.character(rownames(res_kinaseDirection)),
  lab = '',
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = selected_genes,
  pCutoff = 10e-2, # p-value cutoff line
  FCcutoff = 0.5, # fold change cutoff line
  xlab = paste0("<---- ", "ETV6-NTRK3", "   ", "Log2 FoldChange", "   ", "Kinase", " ---->"),
  pointSize = 4.0,
  labSize = 6.0,
  labCol = 'black',
  labFace = 'bold',
  boxedLabels = TRUE,
  # colAlpha = 4/5,
  legendPosition = 'right',
  legendLabSize = 14,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 1.0,
  colConnectors = 'black',
  max.overlaps = Inf,
  caption = "", # Caption is empty as per your original code
  title = paste0("DEGs ", "Kinase", " vs ", "ETV6-NTRK3"), 
  subtitle = "", # Subtitle should be on the same line as title
)



### DEGs - Pathway Enrichmemnt (Gene Ontology) -----------
library(clusterProfiler)
library(org.Hs.eg.db)
res_sig <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Objects/DEGs_kinase_direction_filtered_617genes.rds")
DEGs <- res_sig #they're all > 0.5 FoldChange
DEGs$gene <- rownames(DEGs)
up_genes_kinase <- DEGs[DEGs$log2FoldChange >= 0.5 & DEGs$padj < 0.05, ]$gene
up_genes_etv6 <- DEGs[DEGs$log2FoldChange <= 0.5 & DEGs$padj < 0.05, ]$gene

# Kinase tumor regulated genes - Gene Ontology
ego2_kinase <- enrichGO(gene = up_genes_kinase,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'SYMBOL',
                    # ont	One of "BP", "MF", and "CC" subontologies,                                   # or "ALL" for all three
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)

ego2_kinase@result[order(ego2_kinase@result$Count, decreasing = TRUE),][1:10,]
barplot(ego2_kinase, drop=TRUE, main = "")
#clusterProfiler::dotplot(ego2_up, showCategory=30) + ggtitle("")

## Customized barplot 
ego2_kinase_plot <- ego2_kinase@result[ego2_kinase@result$Description %in% 
                                 c("cilium organization", "cilium assembly", "axoneme assembly", 
                                   "cilium movement", "microtubule-based movement", 
                                   "microtubule bundle formation", "axonemal dynein complex assembly",
                                   "cilium or flagellum-dependent cell motility"), ]

library(ggplot2); theme_set(theme_classic())
ggplot(ego2_kinase_plot, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_col() +
  labs(x = "Gene Count", y = "", title = "Tumor-Kinase Up Regulated Pathways") +
  scale_fill_gradient(low = "red", high = "blue", name = "p.adjust") 

## Network plot 
library(UpSetR)
table(unlist(strsplit(ego2_kinase_plot$geneID, "/")))
gene_list <- strsplit(ego2_kinase_plot$geneID, "/")
names(gene_list) <- ego2_kinase_plot$Description
UpSetR::upset(fromList(gene_list), 
              order.by = "freq", 
              nsets = 8) # Still need more time to interpret this ...
# # Not doing as commented below
# Remove c("cilium assembly", "axome assembly", "microtubele bundle formation", "axonemal dynein complex assembly", "cillium moviment")
# ego2_kinase_plot[!ego2_kinase_plot$Description %in% c("cilium assembly",), ]

ego2_up_sub <- filter(ego2_kinase, Description %in% ego2_kinase$Description)
edox <- setReadable(ego2_up_sub, 'org.Hs.eg.db', 'ENTREZID')
geneList <- up_genes_kinase

cnetplot(edox, foldChange=geneList, 
               #circular = TRUE, colorEdge = TRUE, 
               #showCategory = "", 
) + ggtitle("Tumor-Kinase Up Regulated Pathways ") 

## Save enrichment output table 
# write.csv(ego2_kinase_plot, '/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Objects/ego2_kinase_pathways_table_filtered.csv')



# ETV6 tumor regulated genes - Gene Ontology
ego2_etv6 <- enrichGO(gene = up_genes_etv6,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'SYMBOL',
                      # ont	One of "BP", "MF", and "CC" subontologies,                                   # or "ALL" for all three
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

ego2_etv6@result[order(ego2_etv6@result$Count, decreasing = TRUE),][1:10,]
barplot(ego2_etv6, drop=TRUE, main = "")
# clusterProfiler::dotplot(ego2_etv6, showCategory=30) + ggtitle("")

## Customized barplot 
ego2_etv6_plot <- ego2_etv6@result[ego2_etv6@result$Description %in% c("Wnt signaling pathway",
                                                                       "epithelial cell proliferation",
                                                                       "skeletal system morphogenesis",
                                                                       "embryonic skeletal system development",
                                                                       "bone development"), ]


library(ggplot2); theme_set(theme_classic())
ggplot(ego2_etv6_plot, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_col() +
  labs(x = "Gene Count", y = "", title = "Tumor-Kinase Down Regulated Pathways") +
  scale_fill_gradient(low = "red", high = "blue", name = "p.adjust") 


ego2_down_sub <- filter(ego2_etv6, Description %in% ego2_etv6_plot$Description) #trying to add wnt pathways into network plot
ego2_down_sub@result$p.adjust <- 0.001 # this is only to make it go to the plot !
ego2_down_sub@result$qvalue <- 0.001 # this is only to make it go to the plot !
edox <- setReadable(ego2_down_sub, 'org.Hs.eg.db', 'ENTREZID')
geneList <- up_genes_etv6

cnetplot(edox, foldChange=geneList, 
               #circular = TRUE, colorEdge = TRUE, 
               #showCategory = "", 
) 

## Save enrichment output table 
# write.csv(ego2_etv6_plot, '/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Objects/ego2_etv6_pathways_table_filtered.csv')



### Most Variable Genes (Heatmap - unsupervised clusters) -----------

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
dds <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Objects/dds_RNA.rds")
normalized_counts <- counts(dds, normalized=TRUE)
# Filter normalized_counts by labsheet
normalized_counts_hm <- normalized_counts[, colnames(normalized_counts) %in% labsheet_hm$smpID]

### PLOT - heatmap - whole RNAseq - focus on NTRK label
library(pheatmap)
library(RColorBrewer)
normalized_counts_hm # exp matrix
labsheet_hm # metadata

my_sample_col <- data.frame(row.names = labsheet_hm$smpID, 
                            Group = labsheet_hm$group, 
                            Histology = labsheet_hm$histologic_diagnosis,
                            Fusion = labsheet_hm$fusion_type) # patient ID in row names; any columns are for sample information
my_sample_col[my_sample_col$Fusion %in% "ETV6:NTRK3", ]$Fusion <- "ETV6::NTRK3" #correcting label typo

NTRK <- my_sample_col[grep('NTRK', my_sample_col$Fusion), ]
NTRK_butETV6 <- rownames(NTRK[!NTRK$Fusion %in% c('ETV6::NTRK3'), ])
NTRK_onlyETV6 <- rownames(NTRK[!rownames(NTRK) %in% NTRK_butETV6, ])

my_sample_col$Fusion_concat <- NA
my_sample_col[rownames(my_sample_col) %in% NTRK_butETV6, ]$Fusion_concat <- 'Others_NTRK'
my_sample_col[rownames(my_sample_col) %in% NTRK_onlyETV6, ]$Fusion_concat <- 'ETV6_NTRK3'
my_sample_col[my_sample_col$Fusion_concat %in% NA, ]$Fusion_concat <- 'Not_NTRK'


ann_colors = list(
  Fusion_concat = c(ETV6_NTRK3 = 'black', Not_NTRK = 'blue', Others_NTRK = 'red'))

my_sample_col$Histology <- NULL

# Most variable feature - from PC1 components 
pca_data <- na.omit(normalized_counts_hm)
pca_result <- prcomp(t(pca_data))
summary(pca_result) # "Proportion of Variance" gets in a platou at PC4. So Let's use 1000 CpGs from PC1, 2, 3, and 4. This we will call "Most Variable Feature"
loadings_PC1 <- pca_result$rotation[, 1]
abs_loadings_PC1 <- abs(loadings_PC1)
top_n <- 1000
top_cpg_indices <- order(abs_loadings_PC1, decreasing = TRUE)[1:top_n]
top_cpg_sites_PC1 <- rownames(pca_data)[top_cpg_indices]

loadings_PC2 <- pca_result$rotation[, 2]
abs_loadings_PC2 <- abs(loadings_PC2)
top_n <- 1000
top_cpg_indices <- order(abs_loadings_PC2, decreasing = TRUE)[1:top_n]
top_cpg_sites_PC2 <- rownames(pca_data)[top_cpg_indices]

loadings_PC3 <- pca_result$rotation[, 3]
abs_loadings_PC3 <- abs(loadings_PC3)
top_n <- 1000
top_cpg_indices <- order(abs_loadings_PC3, decreasing = TRUE)[1:top_n]
top_cpg_sites_PC3 <- rownames(pca_data)[top_cpg_indices]

MVC_PC1_2_3 <- c(top_cpg_sites_PC1,
                 top_cpg_sites_PC2,
                 top_cpg_sites_PC3)

MVC_PC1_2_3 <- MVC_PC1_2_3[!duplicated(MVC_PC1_2_3)]

# Correlation test
normalized_counts_hm_sub <- normalized_counts_hm[MVC_PC1_2_3, ]
library(ggplot2)
library(corrplot)

cor_test_mat <- cor.mtest(normalized_counts_hm_sub, conf.level = 0.95)
cor_test_mat <- cor_test_mat$p
cor_matrix <- cor(normalized_counts_hm_sub)

# To add statistics to the heatmap
convert_to_symbols <- function(p_value) {
  if (is.na(p_value) || p_value > 0.05) {
    return("NS")
  } else if (p_value <= 0.05 && p_value > 0.01) {
    return("*")
  } else if (p_value <= 0.01 && p_value > 0.001) {
    return("**")
  } else if (p_value <= 0.001) {
    return("***")
  }
}
# Assuming your matrix is named 'p_value_matrix'
# Apply the function to each element of the matrix
pvalue_symbol_matrix <- apply(cor_test_mat, c(1, 2), convert_to_symbols)

# Adding cor stats 
pheatmap(cor_matrix[,], #smp_ID on columns and rows
         #annotation_row = my_probe_col, 
         annotation_col = my_sample_col[, ],
         show_rownames = FALSE,
         cluster_columns = TRUE,
         #scale = "column",
         #scale = "row",
         scale = "none", #the actual beta-value 0-1.
         #annotation_colors = ann_colors,
         display_numbers = pvalue_symbol_matrix,
         fontsize_number = 10,
         annotation_colors = ann_colors,
         main = "Most Variable Genes from PC1,2,3 - Correlation - scale:none")

# No cor stats 
pheatmap(cor_matrix[,], #smp_ID on columns and rows
         #annotation_row = my_probe_col, 
         annotation_col = my_sample_col[, ],
         show_rownames = FALSE,
         cluster_columns = TRUE,
         #scale = "column",
         #scale = "row",
         scale = "none", #the actual beta-value 0-1.
         #annotation_colors = ann_colors,
         #display_numbers = pvalue_symbol_matrix,
         fontsize_number = 10,
         annotation_colors = ann_colors,
         main = "Most Variable Genes from PC1,2,3 - Correlation - scale:none")








