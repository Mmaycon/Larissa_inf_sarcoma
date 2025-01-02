### --- Differnetial Expression (DE) Analysis --- ###

# Load packages -------------
library(DESeq2)
library(sva)
library(dplyr)
library(stringr)


# 1. Load the data ----------------------
# Same data used on ~/DE_analysis_V2.R
# Obs: saving it under /maycon directory
# load("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/reviewed_smpID_countmtx_labsheet.rda")
# save(labsheet, count_matrix, file = '/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/bulk_RNA/objects/Round1/reviewed_smpID_countmtx_labsheet.rda')

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
# make it identical (same sequence)
count_matrix <- count_matrix[, labsheet$smpID]
identical(colnames(count_matrix), labsheet$smpID) #TRUE
# add rownames into labsheet
rownames(labsheet) <- labsheet$smpID



# 2. Run DE ----------------------
count_matrix <- count_matrix[, !colnames(count_matrix) %in% "SJST034534_D1"]
dim(count_matrix) #25 smp
labsheet <- labsheet[!labsheet$smpID %in% "SJST034534_D1", ]
dim(labsheet)#25 smp
# make it identical (same sequence)
count_matrix <- count_matrix[, labsheet$smpID]
identical(colnames(count_matrix), labsheet$smpID) #TRUE
# add rownames into labsheet
rownames(labsheet) <- labsheet$smpID
# Create DESEq object 
dds <- DESeqDataSetFromMatrix(count_matrix, colData = labsheet, design = ~ smp_type + group)
# Important step for further normalization during DESEq()
dds <- estimateSizeFactors(dds)
# sizeFactors(dds) #Idk what it does ... neither if we need to run it
# Pre-filtering 
table(dds$group)
smallestGroupSize <- 8
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
# indicate direction for comparison
dds$condition <- relevel(dds$group, ref = "kinase-fused_tumor")
# Run Differential Expression analysis
dds <- DESeq(dds)
# Filter results for significant genes
res <- results(dds)
res <- data.frame(res)
# Subset stat. significant genes
res_pos <- subset(res, res$padj < 0.05 & res$log2FoldChange >= 0.5)
res_neg <- subset(res, res$padj < 0.05 & res$log2FoldChange <= -0.5)
res_sig <- rbind(res_pos, res_neg)
dim(res_sig) #617   6
# PLOT - FoldChange distribution - histogram
hist(res_sig$log2FoldChange, 
     main = "FoldChange dist. (s.sig DEGs)", 
     xlab = "FoldChange",
     ylab = "N of genes")

# 3. Visualize DE results -------------
# Perform PCA
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_DEGs <- normalized_counts[rownames(res_sig), ]

pca <- prcomp(t(normalized_counts_DEGs)) 
aux <- as.data.frame(pca$x[, 1:3]) 
# Merge data
scores <- merge(labsheet, aux, by.y=0, by.x="smpID", all.x=T)
# Make a column for histologic diagnosis to be plotted over outilers samples only
scores$histologic_for_outlier <- 'Others'
outlier_smp <- c('SJST033491_D1',
                 'SJST032952_D2',
                 'SJST032952_D4',
                 'SJST032767_D2',
                 'SJST033312_D1',
                 'SJST033835_D1',
                 'SJST031920_D1')
his_diag <- scores[scores$smpID %in% outlier_smp, ]$histologic_diagnosis

scores[scores$smpID %in% outlier_smp[1], ]$histologic_for_outlier <- his_diag[1]
scores[scores$smpID %in% outlier_smp[2], ]$histologic_for_outlier <- his_diag[2]
scores[scores$smpID %in% outlier_smp[3], ]$histologic_for_outlier <- his_diag[3]
scores[scores$smpID %in% outlier_smp[4], ]$histologic_for_outlier <- his_diag[4]
scores[scores$smpID %in% outlier_smp[5], ]$histologic_for_outlier <- his_diag[5]
scores[scores$smpID %in% outlier_smp[6], ]$histologic_for_outlier <- his_diag[6]
scores[scores$smpID %in% outlier_smp[7], ]$histologic_for_outlier <- his_diag[7]

# PLOT - PCA within DEGs only 
library(ggplot2); theme_set(theme_classic())
# Create PCA plot for all samples separated by articles
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
  ggtitle("Normalized / batch corrected (mtx design) - only DEGs (617)") 


library(EnhancedVolcano)
# PLOT - volcano plot - DEGs
#rownames(res_sig) <- res_sig$SYMBOL # don't run it if it doesn't need to
EnhancedVolcano(res_sig,
                lab = "",
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c("FLG"  , "SLC6A15"    ,  "ALK"   ,  "CTXND1"  ,    "PAX3"),
                pCutoff = 10e-2, #pvalue cutoff line
                FCcutoff = 2.0, #foldChange cutoff line
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                #boxedLabels = TRUE,
                #colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                #drawConnectors = TRUE,
                #widthConnectors = 1.0,
                #colConnectors = 'black',
                title = "kinase-fused_tumor direction")


# 4. Saving outputs 
setwd("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/bulk_RNA/objects/Round2")

saveRDS(scores, './PCA_DEGs_plot.rds')
DEGs <- res_sig
saveRDS(DEGs, './Volcano_DEGs_plot.rds')
saveRDS(dds, './DESeq2_output.rds')


# 5. Go on enricher -> find kinase paths -> run it
library(writexl)
DEGs$genes <- rownames(DEGs)
hist(DEGs$log2FoldChange)
top_up_genes <- DEGs[DEGs$log2FoldChange >= 2, ]
write_xlsx(top_up_genes,"./top_up_genes_FC2_kinase_orient.xlsx")
top_down_genes <- DEGs[DEGs$log2FoldChange <= 2, ]
write_xlsx(top_down_genes,"./top_down_genes_FC2_kinase_orient.xlsx")

/DE_V3.R



