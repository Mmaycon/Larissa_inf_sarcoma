# Tidying the data -----------
DE_input_matrix <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/DE_input_matrix.rds")
dim(DE_input_matrix) #24 smp
metadata <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/metadata_labsheet.rds")
dim(metadata) #25 smp
identical(colnames(DE_input_matrix), metadata$smpID) #FALSE



colnames(DE_input_matrix) <- make.names(colnames(DE_input_matrix), unique = TRUE)
metadata$smpID <- make.names(metadata$smpID, unique = TRUE)
setdiff(colnames(DE_input_matrix), metadata$smpID)
setdiff(metadata$smpID, colnames(DE_input_matrix))


metadata <- metadata[!metadata$smpID %in% c("SJST033312_D1.1", "SJST033491_D1"), ]
DE_input_matrix <- DE_input_matrix[, metadata$smpID]
identical(colnames(DE_input_matrix), metadata$smpID) #TRUE

# Diff. Expression --------
# Ordering the samples by group of comparison in the DE_input_matrix
etv6 <- metadata[metadata$group %in% 'ETV6_NTRK3_fused', ]$smpID
kinase <- metadata[metadata$group %in% 'Kinase_fused', ]$smpID
etv6_mt <- DE_input_matrix[, colnames(DE_input_matrix) %in% etv6]
length(colnames(etv6_mt)) #8 samples
kinase_mt <- DE_input_matrix[, colnames(DE_input_matrix) %in% kinase]
length(colnames(kinase_mt)) #16 samples
DE_input_matrix <- cbind(etv6_mt,kinase_mt)

# Comparing groups - we get pvalues from this step
#install.packages('exactRankTests')
require(exactRankTests)
require(parallel)
values <- t(DE_input_matrix) #transpose DE_input_matrix
values <- data.frame(values)
wpvalues <- unlist(mclapply(values,
                              function(gene) {
                                zz <- wilcox.exact(gene[1:  dim(etv6_mt)[2]],
                                                   gene[c(dim(etv6_mt)[2]+1) : dim(DE_input_matrix)[2]], exact=T) # excat = T nos da um resultado mais preciso, mas demora mais pra rodar o looping. Entao, exact = F é melhor custo-beneficio. . .
                                z <- zz$p.value
                                return(z)
                              }, mc.cores= 20))
wpvalues_adj <- p.adjust(wpvalues, method = "BH")
wpvalues_adj <- data.frame(wpvalues_adj)
wpvalues_adj$gene <- rownames(wpvalues_adj)
hist(wpvalues_adj$wpvalues_adj)
table(wpvalues_adj$wpvalues_adj < 0.05) # no genes  ....

# Compute FoldChange
# Calculate the mean expression for each gene across conditions
mean_expr_etv6 <- rowMeans(DE_input_matrix[, etv6]) #row = genes; columns = samples
mean_expr_kinase <- rowMeans(DE_input_matrix[, kinase])
# Calculate Fold Change for each gene
fold_change <- mean_expr_etv6 / mean_expr_kinase # interpret it as "etv6 group has n times more expression upon a given gene compared to kinase group
fold_change <- data.frame(fold_change) 
fold_change$gene <- rownames(fold_change)
hist(fold_change$fold_change, breaks = 1000)


FG_pvalue <- merge(fold_change, wpvalues_adj, by='gene')
volcano <- data.frame(DE_input_matrix)
volcano$gene <- rownames(volcano)
volcano <- merge(volcano, FG_pvalue, by='gene')
head(volcano)
table(volcano$wpvalues_adj < 0.05)

# Create labels for the volcano plot 
volcano$threshold <- "1" # threshold info will be used to define DMP and to color volcano plots
b <- volcano[volcano$wpvalues_adj < 0.05 & volcano$fold_change < -0.05,] #hypo NSC probes
volcano[rownames(b),"threshold"] <- "2"
c <- volcano[volcano$wpvalues_adj < 0.05 & volcano$fold_change > 0.05,] #hyper NSC probes 
volcano[rownames(c),"threshold"] <- "3"
print(table(volcano$threshold))

table(volcano$wpvalues_adj < 0.05) # none statistical significance =o make it again ...




### DESEq2 - let's proceed with the raw counts anyway ### ---------
library(DESeq2)
library(stringr)
library(dplyr)

setwd("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis")

# Load and tide the raw data
## Count matrix
count_matrix <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/count_matrix.rds")
colnames(count_matrix) <- str_split_fixed(as.character(colnames(count_matrix)), "[-]", 2)[,1]
## Metadata
metadata <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/metadata_labsheet.rds")
metadata_1 <- readRDS('/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/metadata.rds')
metadata_1$smpID <- str_split_fixed(as.character(metadata_1$smpID), "[-]", 2)[,1]
metadata <- merge(metadata_1, metadata[, !colnames(metadata) %in% ('group')], by='smpID')
## Match samples on both
identical(colnames(count_matrix), metadata$smpID) #FALSE
colnames(count_matrix) <- make.names(colnames(count_matrix), unique = TRUE)
metadata$smpID <- make.names(metadata$smpID, unique = TRUE)
rm_it_1 <- setdiff(colnames(count_matrix), metadata$smpID)
rm_it_2 <- setdiff(metadata$smpID, colnames(count_matrix))
remove_smp <- c(rm_it_1, rm_it_2)
metadata <- metadata[!metadata$smpID %in% remove_smp, ]
count_matrix <- count_matrix[, !colnames(count_matrix) %in% remove_smp]
## Just letting samples to be in the same sequence (optional. Not necessary)
count_matrix <- count_matrix[, metadata$smpID]
identical(colnames(count_matrix), metadata$smpID) #TRUE
## Remove a not well aligned sample  and an outlier 
count_matrix <- count_matrix[, !colnames(count_matrix) %in% c("SJST033312_D1-S2", "SJST033491_D1")]
metadata <- metadata[!metadata$smpID %in% c("SJST033312_D1-S2", "SJST033491_D1"), ]
dim(count_matrix) #63568    25
dim(metadata) # 25  6
 
# Run DESEq workflow 
## Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(count_matrix, colData = metadata, design = ~ group)
## pre-filtering
table(metadata$group)
smallestGroupSize <-  9 #only 9 samples on the smallest group (ETV6)
keep <- rowSums(counts(dds) >= 1000) >= smallestGroupSize
dds <- dds[keep,]
## indicate direction for comparison
dds$condition <- relevel(dds$group, ref = "ETV6_NTRK3_fused")

## Run Differential Expression analysis
dds <- DESeq(dds)

## Filter results for significant genes
res <- results(dds)
res <- data.frame(res)


res_pos <- subset(res, res$padj < 0.05 & res$log2FoldChange >= 0.5)
res_neg <- subset(res, res$padj < 0.05 & res$log2FoldChange <= -0.5)
res_sig <- rbind(res_pos, res_neg)
dim(res_sig) #1088    6
hist(res_sig$log2FoldChange)

## Removing "." from gene IDs
rownames(res_sig) <- str_split_fixed(as.character(rownames(res_sig)), "[.]", 2)[,1] # good it doesn't have duplicates ENSEMBL here

# Volcano plot ------------
# source from https://github.com/kevinblighe/EnhancedVolcano
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)

# Covert gene ID
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

ens <- rownames(res_sig) #genes from DE
library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')

symbols <- symbols[!is.na(symbols)]

symbols <- data.frame(ENSEMBL = names(symbols),
           SYMBOL = as.vector(symbols))

res_sig$ENSEMBL <- rownames(res_sig)

res_sig <- merge(res_sig, symbols, by='ENSEMBL')
dim(res_sig)# 385   8 

res_sig[order(res_sig$log2FoldChange, decreasing = TRUE), c('SYMBOL', 'padj', 'log2FoldChange')]


# Plot volcano plot - basics 
rownames(res_sig) <- res_sig$SYMBOL
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


# PCA - 6,620 feature data
# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(count_matrix, colData = metadata, design = ~ group)
# normalization
keep <- rowSums(counts(dds) >= 1000) >= smallestGroupSize
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
dim(normalized_counts) #6620   25

# Perform PCA
pca <- prcomp(t(normalized_counts)) 
aux <- as.data.frame(pca$x[, 1:3]) 
# Merge data
scores <- merge(metadata, aux, by.y=0, by.x="smpID", all.x=T)
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
  ggtitle("PCA normalized data (6620 genes before DE and 25smp)") 




# PCA - DE feature data
# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(count_matrix, colData = metadata, design = ~ group)
# normalization
keep <- rowSums(counts(dds) >= 1000) >= smallestGroupSize
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
dim(normalized_counts) #6620   25
normalized_counts

# Chanching ENSEMBL -> SYMBOL
rownames(normalized_counts) <- str_split_fixed(as.character(rownames(normalized_counts)), "[.]", 2)[,1] # good it doesn't have duplicates ENSEMBL here

ens <- rownames(normalized_counts) #genes from DE
library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')

symbols <- symbols[!is.na(symbols)]

symbols <- data.frame(ENSEMBL = names(symbols),
                      SYMBOL = as.vector(symbols))
normalized_counts <- data.frame(normalized_counts)
normalized_counts$ENSEMBL <- rownames(normalized_counts)

normalized_counts <- merge(normalized_counts, symbols, by='ENSEMBL')
dim(normalized_counts)#6351   27

normalized_counts_DE <- normalized_counts[normalized_counts$SYMBOL %in% res_sig$SYMBOL, ]


# Perform PCA
normalized_counts_DE$ENSEMBL <- NULL
normalized_counts_DE$SYMBOL <- NULL
dim(normalized_counts_DE) #374  25
pca2 <- prcomp(t(as.matrix(normalized_counts_DE))) 
aux <- as.data.frame(pca2$x[, 1:3]) 
# Merge data
scores2 <- merge(metadata, aux, by.y=0, by.x="smpID", all.x=T)
# Load necessary library and set theme
library(ggplot2); theme_set(theme_classic())
# Create PCA plot for all samples separated by articles
ggplot(scores2, aes(x=PC1, y=PC2, colour=factor(group), )) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'), name="Group") +
  xlab(paste0("PC1 (", prettyNum(summary(pca2)$importance[2,1]*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca2)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  ggtitle("PCA normalized data (374 genes after DE and 25smp)") 




