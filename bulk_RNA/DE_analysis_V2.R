### --------------- DE analysis V. 2 --------------- ###

# Load the data ------------
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


# PCA of raw counts ------------------
# make it identical (same sequence)
count_matrix <- count_matrix[, labsheet$smpID]
identical(colnames(count_matrix), labsheet$smpID) #TRUE
# add rownames into labsheet
rownames(labsheet) <- labsheet$smpID

# Filtering low count genes
count_matrix <- as.matrix(count_matrix)
dds <- DESeqDataSetFromMatrix(count_matrix, colData = labsheet, design = ~ group)
# pre-filtering
table(metadata$group)
smallestGroupSize <-  9 #only 9 samples on the smallest group (ETV6)
keep <- rowSums(counts(dds) >= 1000) >= smallestGroupSize
dds <- dds[keep,]
raw_counts <- counts(dds, normalized=FALSE)
dim(raw_counts) #6446   26

# Perform PCA - removing outlier
count_matrix <- data.frame(count_matrix)
pca <- count_matrix[, !colnames(count_matrix) %in% "SJST034534_D1"]
pca <- prcomp(t(pca)) 
aux <- as.data.frame(pca$x[, 1:3]) 
# Merge data
scores <- merge(labsheet, aux, by.y=0, by.x="smpID", all.x=T)
# Load necessary library and set theme
library(ggplot2); theme_set(theme_classic())
# Create PCA plot for all samples separated by articles
ggplot(scores, aes(x=PC1, y=PC2, colour=factor(group), )) +
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
  geom_text_repel(aes(label = smpID)) +
  ggtitle("PCA raw data (6446 genes before DE and 25smps(23unique))") 


# PCA of normalized counts ------------------
# remove outlier
count_matrix <- count_matrix[, !colnames(count_matrix) %in% "SJST034534_D1"]
labsheet <- labsheet[!labsheet$smpID %in% "SJST034534_D1", ]
# make it identical (same sequence)
count_matrix <- count_matrix[, labsheet$smpID]
identical(colnames(count_matrix), labsheet$smpID) #TRUE
# add rownames into labsheet
rownames(labsheet) <- labsheet$smpID

# Filtering low count genes
count_matrix <- as.matrix(count_matrix)
dds <- DESeqDataSetFromMatrix(count_matrix, colData = labsheet, design = ~ group)
# pre-filtering
table(metadata$group)
smallestGroupSize <-  9 #only 9 samples on the smallest group (ETV6)
keep <- rowSums(counts(dds) >= 1000) >= smallestGroupSize
dds <- dds[keep,]
# Normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
dim(normalized_counts) #6082   25
# Perform PCA
pca <- prcomp(t(normalized_counts)) 
aux <- as.data.frame(pca$x[, 1:3]) 
# Merge data
scores <- merge(labsheet, aux, by.y=0, by.x="smpID", all.x=T)
# Load necessary library and set theme
library(ggplot2); theme_set(theme_classic())
# Create PCA plot for all samples separated by articles
ggplot(scores, aes(x=PC1, y=PC2, colour=factor(group), )) +
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
  geom_text_repel(aes(label = smpID)) +
  ggtitle("PCA normalized data (6082 genes before DE and 25smps(23unique))") 



# Differential Expression  ------------------
# remove outlier
count_matrix <- count_matrix[, !colnames(count_matrix) %in% "SJST034534_D1"]
labsheet <- labsheet[!labsheet$smpID %in% "SJST034534_D1", ]
# make it identical (same sequence)
count_matrix <- count_matrix[, labsheet$smpID]
identical(colnames(count_matrix), labsheet$smpID) #TRUE
# add rownames into labsheet
rownames(labsheet) <- labsheet$smpID
# Filtering low count genes
count_matrix <- as.matrix(count_matrix)
dds <- DESeqDataSetFromMatrix(count_matrix, colData = labsheet, design = ~ group)
# pre-filtering
table(metadata$group)
smallestGroupSize <-  9 #only 9 samples on the smallest group (ETV6)
keep <- rowSums(counts(dds) >= 1000) >= smallestGroupSize
dds <- dds[keep,]
# indicate direction for comparison
dds$condition <- relevel(dds$group, ref = "ETV6-NTRK3_fused_tumor")
# Run Differential Expression analysis
dds <- DESeq(dds)
# Filter results for significant genes
res <- results(dds)
res <- data.frame(res)
# Subset stat. significant genes
res_pos <- subset(res, res$padj < 0.05 & res$log2FoldChange >= 0.5)
res_neg <- subset(res, res$padj < 0.05 & res$log2FoldChange <= -0.5)
res_sig <- rbind(res_pos, res_neg)
dim(res_sig) #444   6
hist(res_sig$log2FoldChange)
res_sig <- data.frame(res_sig)

res_sig[order(res_sig$log2FoldChange, decreasing = TRUE), ]

library(EnhancedVolcano)
# Plot volcano plot - basics 
#rownames(res_sig) <- res_sig$SYMBOL # don't run it if it doesn't need to
EnhancedVolcano(res_sig,
                lab = rownames(res_sig),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c('SORCS1', 'VAT1L', 'RELN', 'PPL', 'ART3'),
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




# PCA of normalized counts DEGs only------------------
# remove outlier
count_matrix <- count_matrix[, !colnames(count_matrix) %in% "SJST034534_D1"]
labsheet <- labsheet[!labsheet$smpID %in% "SJST034534_D1", ]
# make it identical (same sequence)
count_matrix <- count_matrix[, labsheet$smpID]
identical(colnames(count_matrix), labsheet$smpID) #TRUE
# add rownames into labsheet
rownames(labsheet) <- labsheet$smpID

# Filtering low count genes
count_matrix <- as.matrix(count_matrix)
dds <- DESeqDataSetFromMatrix(count_matrix, colData = labsheet, design = ~ group)
# pre-filtering
table(metadata$group)
smallestGroupSize <-  9 #only 9 samples on the smallest group (ETV6)
keep <- rowSums(counts(dds) >= 1000) >= smallestGroupSize
dds <- dds[keep,]
# Normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
dim(normalized_counts) #6082   25
# Filter by DEGs
normalized_counts <- normalized_counts[rownames(normalized_counts) %in% rownames(res_sig), ]
dim(normalized_counts) #444  25

# Perform PCA
pca <- prcomp(t(normalized_counts)) 
aux <- as.data.frame(pca$x[, 1:3]) 
# Merge data
scores <- merge(labsheet, aux, by.y=0, by.x="smpID", all.x=T)
# Load necessary library and set theme
library(ggplot2); theme_set(theme_classic())
# Create PCA plot for all samples separated by articles
ggplot(scores, aes(x=PC1, y=PC2, colour=factor(group), )) +
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
  ggtitle("PCA normalized data (444 genes after DE and 25smps(23unique))") 

install.packages('writexl')
library(writexl)
res_sig$genes <- rownames(res_sig)
write_xlsx(res_sig[res_sig$log2FoldChange >= 0.5, ],"/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/FG05_high_176_genes.xlsx")
write_xlsx(res_sig[res_sig$log2FoldChange <= -0.5, ],"/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/FG05_low_268_genes.xlsx")

saveRDS(res_sig[res_sig$log2FoldChange >= 0.5, ], '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/FG05_high_176_genes.rds')
saveRDS(res_sig[res_sig$log2FoldChange <= -0.5, ], '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/FG05_low_268_genes.rds')


