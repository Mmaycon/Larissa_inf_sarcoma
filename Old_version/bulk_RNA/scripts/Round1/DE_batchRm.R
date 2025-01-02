# Find/remove batch effects within DESEq workflow 

# Material sources -------------
# https://support.bioconductor.org/p/121408/
# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#using-sva-with-deseq2


# Load packages -------------
library(DESeq2)
library(sva)
library(dplyr)
library(stringr)
# Load the data  ------------
# Same data used on ~/DE_analysis_V2.R
load("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/reviewed_smpID_countmtx_labsheet.rda")
labsheet
count_matrix
labsheet$group<- gsub(' ', '_', labsheet$group)

# Convert gene ID -----------
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

# Create DESEq object 
dds <- DESeqDataSetFromMatrix(count_matrix, colData = labsheet, design = ~ group)
# Important step for further normalization during DESEq()
dds <- estimateSizeFactors(dds)
# sizeFactors(dds) #Idk what it does ... neither if we need to run it
# Pre-filtering 
table(dds$group)
smallestGroupSize <- 8
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Detect/add batch effect to the matrix design
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) >= 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ group, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)
print(svseq$sv)

# These are like boxplots
# Meant to visualize variation across your candidate batch effect variable. In this case it is smp_type 
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds$smp_type, #batch effect variable 
             vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + group
ddssva <- estimateSizeFactors(ddssva)


# DE analysis after batch removal -------
# indicate direction for comparison
ddssva$condition <- relevel(ddssva$group, ref = "kinase-fused_tumor")
# Run Differential Expression analysis
ddssva <- DESeq(ddssva)
# Filter results for significant genes
res <- results(ddssva)
res <- data.frame(res)
# Subset stat. significant genes
res_pos <- subset(res, res$padj < 0.05 & res$log2FoldChange >= 0.5)
res_neg <- subset(res, res$padj < 0.05 & res$log2FoldChange <= -0.5)
res_sig <- rbind(res_pos, res_neg)
dim(res_sig) #1069    6
hist(res_sig$log2FoldChange)
res_sig_batchrm <- data.frame(res_sig)
res_sig_batchrm


# DE analysis before batch removal -------
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

# Normalization  -----------
# Create DESEq object 
dds <- DESeqDataSetFromMatrix(count_matrix, colData = labsheet, design = ~ group)
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
dim(res_sig) #755    6
hist(res_sig$log2FoldChange)
res_sig_nobatchcor <- data.frame(res_sig)
res_sig_nobatchcor


# DE analysis after batch removal BY MATRIX DESIGN -------
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
hist(res_sig$log2FoldChange)
res_sig_batch_design <- data.frame(res_sig)
res_sig_batch_design


# Perform PCA
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_DEGs <- normalized_counts[rownames(res_sig_batch_design), ]
raw_counts_DEGs <- count_matrix[rownames(res_sig_batch_design), ]

# pca <- prcomp(t(normalized_counts)) 
pca <- prcomp(t(normalized_counts_DEGs)) 
# pca <- prcomp(t(count_matrix)) 
# pca <- prcomp(t(raw_counts_DEGs)) 
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
# Load necessary library and set theme
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
  ggtitle("Normalized / batch corrected data - only DEGs (617)") 



library(EnhancedVolcano)
# Plot volcano plot - basics 
#rownames(res_sig_batch_design) <- res_sig_batch_design$SYMBOL # don't run it if it doesn't need to
EnhancedVolcano(res_sig_batch_design,
                lab = rownames(res_sig_batch_design),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c("FLG"  , "SLC6A15"    ,  "ALK"   ,  "CTXND1"  ,    "PAX3"),
                pCutoff = 10e-2, #pvalue cutoff line
                FCcutoff = 2.0, #foldChange cutoff line
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                #colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')









# Create a function to compare DEGs from different normalization methods ---------
# *Using DESEq2 output

# DEGs
df1 <- res_sig_nobatchcor
dim(df1) # 950   6
df2 <- res_sig_batchrm
dim(df2) #1329    6
df3 <- res_sig_batch_design
dim(df3) #464   6

# length(rownames(df1)) #950
# length(rownames(df2)) #1329
# intersect_df1_per <- length(intersect(rownames(df1), rownames(df2)))/length(rownames(df1))
# intersect_df2_per <- length(intersect(rownames(df1), rownames(df2)))/length(rownames(df2))

### Better not use it ### 
# # Function to add "rank" by foldchange
# rank_by_logFC <- function(df) {
# positive_DEGs <- df[df$log2FoldChange > 0, ]
# positive_DEGs <- positive_DEGs[order(positive_DEGs$log2FoldChange, decreasing = TRUE),  ]
# positive_DEGs$rank <- 1:length(rownames(positive_DEGs))
# 
# negative_DEGs <- df[df$log2FoldChange < 0, ]
# negative_DEGs <- negative_DEGs[order(negative_DEGs$log2FoldChange, decreasing = FALSE),  ]
# negative_DEGs$rank <- 1:length(rownames(negative_DEGs))
# 
# return(rbind(positive_DEGs, negative_DEGs))
# }
# 
# df1_ranked <- rank_by_logFC(df1)
# df1_ranked[df1_ranked$rank < 10, ] 
# 
# df2_ranked <- rank_by_logFC(df2)
# df2_ranked[df2_ranked$rank < 10, ] 
# 
# df3_ranked <- rank_by_logFC(df3)
# df3_ranked[df3_ranked$rank < 10, ] 



df[df$log2FoldChange > 0, ]

# Upset plot 
install.packages('ComplexUpset')
library(ComplexUpset)

# # Prepareing the input #01
# listInput <- list(
#   DEG_nobatchcor = rownames(df1_ranked[df1_ranked$rank < 100, ] 
#   ),
#   DEG_batchrm_sva = rownames(df2_ranked[df2_ranked$rank < 100, ] 
#   ),
#   DEF_batch_design = rownames(df3_ranked[df3_ranked$rank < 100, ] ))

# list_names <- names(listInput)
# genes <- unlist(listInput)
# df <- data.frame(list_name = rep(list_names, sapply(listInput, length)),
#                  gene = genes,
#                  stringsAsFactors = FALSE)
# #rownames(df) <- 1:dim(df)[1]
# 
# df <- as.data.frame.matrix(table(df$list_name, df$gene))
# 
# for (col in names(df)) {
#   df[[col]] <- ifelse(df[[col]] == 1, TRUE, FALSE)
# }
# 
# genes <- colnames(df)
# df$DEG_list <- rownames(df)
# 
# # Run upset plot
# ComplexUpset::upset(df, genes,
#                     name='DEG_list',
#                     base_annotations=list(
#                       'Number of Genes'=intersection_size(
#                         counts=FALSE,
#                         mapping=aes(fill=DEG_list)
#                       )  + scale_fill_manual(values=c(
#                         rainbow(length(df$DEG_list))
#                       ))),
#                     width_ratio=0.1)
# upsetplot_pal <- rainbow(length(df$DEG_list))


# Prepareing the input #02
up_top_50_genes_df1 <- rownames(df1[order(df1$log2FoldChange, decreasing = TRUE), ])[1:50]
up_top_50_genes_df2 <- rownames(df2[order(df2$log2FoldChange, decreasing = TRUE), ])[1:50]
up_top_50_genes_df3 <- rownames(df3[order(df3$log2FoldChange, decreasing = TRUE), ])[1:50]


listInput <- list(
  DEG_nobatchcor = up_top_50_genes_df1,
  DEG_batchrm_sva = up_top_50_genes_df2,
  DEG_batch_design =  up_top_50_genes_df3)

length(listInput[[1]])
length(listInput[[2]])
length(listInput[[3]])

list_names <- names(listInput)
genes <- unlist(listInput)
df <- data.frame(list_name = rep(list_names, sapply(listInput, length)),
                 gene = genes,
                 stringsAsFactors = FALSE)
#rownames(df) <- 1:dim(df)[1]

df <- as.data.frame.matrix(table(df$gene, df$list_name))
df$genes <- rownames(df)
DEG_list <- colnames(df)[1:3]
# Run upset plot - @@@@ it is missing some genes IDK why @@@@#@#@#@#$%#$@#!@#@$#%$@#$%
length(df$genes) #why there's fer genes??
blue_palette <- colorRampPalette(c("lightblue", "darkblue"))(367)
ComplexUpset::upset(df, DEG_list, 
                    name='Genes', 
                    base_annotations=list(
                      'DEGenes'=intersection_size(
                        counts=FALSE,
                        #mapping=aes(fill=genes)
                        )  
                      + scale_fill_manual(values=blue_palette)
                      ), 
                    width_ratio=0.1)



install.packages('writexl')
library(writexl)
write_xlsx(data.frame(up_top_50_genes_df3),"/media/ResearchHome/plummgrp/home/common/LFurtado-colab/scripts_git/bulk_RNA/up_top_50_genes_kinase_orient.xlsx")

df3$genes <- rownames(df3)
write_xlsx(data.frame(df3[df3$log2FoldChange > 2, ]),"/media/ResearchHome/plummgrp/home/common/LFurtado-colab/scripts_git/bulk_RNA/DEG_UPgenes_batch_design.xlsx")
write_xlsx(data.frame(df3[df3$log2FoldChange < 2, ]),"/media/ResearchHome/plummgrp/home/common/LFurtado-colab/scripts_git/bulk_RNA/DEG_DOWNgenes_batch_design.xlsx")

saveRDS(df3, '/media/ResearchHome/plummgrp/home/common/LFurtado-colab/scripts_git/bulk_RNA/DEG_batch_design.rds')



# Heatmap 
install.packages("pheatmap")
library(pheatmap)
raw_counts_DEGs  # exp matrix
normalized_counts_DEGs # exp matrix
labsheet # metadata

my_sample_col <- data.frame(row.names = labsheet$smpID, 
                            Group = labsheet$group, 
                            Sample_type = labsheet$smp_type) # patient ID in row names; any columns are for sample information

pheatmap(normalized_counts_DEGs, 
             #annotation_row = my_probe_col, 
             annotation_col = my_sample_col,
             show_rownames = FALSE,
             scale = "row",
             main = 'gene expression visualization')






