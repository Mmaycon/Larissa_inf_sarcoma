# Load DESeq2 library
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)

# Import count data
count_matrix <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/count_matrix.rds")
metadata <- readRDS('/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/metadata.rds')
setwd("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis")

# Remove a low quality sample - I know it because of multiqc from .bam files
count_matrix <- count_matrix[, !colnames(count_matrix) %in% "SJST033312_D1-S2"]
metadata <- metadata[!metadata$smpID %in% "SJST033312_D1-S2", ]

# Basic DESEq2 workflow -------------
# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(count_matrix, colData = metadata, design = ~ group)
# pre-filtering
smallestGroupSize <- 8 #only 8 samples on the smallest group (ETV6)
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
# indicate direction for comparison
dds$condition <- relevel(dds$group, ref = "ETV6_NTRK3_fused")

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds)

# Filter results for significant genes
resSig <- res[which(res$padj < 0.05), ]

# Visualize results
plotMA(res, main = "DESeq2 MA Plot")
plotCounts(dds, gene = rownames(res)[1])

plotMA(resSig, main = "DESeq2 MA Plot")
plotCounts(dds, gene = rownames(resSig)[1])

# Save results to a file
# write.csv(resSig, file = "deseq2_results.csv")




# PCA -------------------
setwd("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis")
# 1# rwa counts
# Perform PCA
pca1 <- prcomp(t(count_matrix)) 
aux <- as.data.frame(pca1$x[, 1:3]) 
# Merge data
scores1 <- merge(metadata, aux, by.y=0, by.x="smpID", all.x=T)
# Load necessary library and set theme
library(ggplot2); theme_set(theme_classic())
# Create PCA plot for all samples separated by articles
ggplot(scores1, aes(x=PC1, y=PC2, colour=factor(group), )) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'), name="Group") +
  xlab(paste0("PC1 (", prettyNum(summary(pca1)$importance[2,1]*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca1)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  ggtitle("PCA raw counts") 
# Save results
# saveRDS(scores1, file = "./PCA_raw_counts.rds")


# 2# normalized counts (from dds object)
# Get normalized data
# source from https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(count_matrix, colData = metadata, design = ~ group)
# normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
dim(normalized_counts)

# Perform PCA
pca2 <- prcomp(t(normalized_counts)) 
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
  ggtitle("PCA normalized data") 
# Save results
# saveRDS(scores2, file = "./PCA_normalized_data.rds")


# 3# normalized + pre-filtering 
# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(count_matrix, colData = metadata, design = ~ group)
# pre-filtering
smallestGroupSize <- 8 #only 8 samples on the smallest group (ETV6)
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
# indicate direction for comparison
dds$condition <- relevel(dds$group, ref = "ETV6_NTRK3_fused")
# normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
dim(normalized_counts) #24888    26
# Perform PCA
pca3 <- prcomp(t(normalized_counts)) 
aux <- as.data.frame(pca3$x[, 1:3]) 
# Merge data
scores3 <- merge(metadata, aux, by.y=0, by.x="smpID", all.x=T)
# Load necessary library and set theme
library(ggplot2); theme_set(theme_classic())
library(ggrepel)
# Create PCA plot for all samples separated by articles
ggplot(scores3, aes(x=PC1, y=PC2, colour=factor(group), )) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'),name="Group") +
  xlab(paste0("PC1 (", prettyNum(summary(pca3)$importance[2,1]*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca3)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  ggtitle("PCA normalized + filtering genes") +
  geom_text_repel(aes(label = smpID))
  
# Save results
# saveRDS(scores3, file = "./PCA_normalized_filteredgenes.rds")


#4 normalized + pre-filtering - SJST033491_D1-S2 (PCA outlier)
# Create a DESeqDataSet object
count_matrix_sub <- count_matrix[, !colnames(count_matrix) %in% 'SJST033491_D1-S2']
metadata_sub <- metadata[!metadata$smpID %in% 'SJST033491_D1-S2', ]

dds <- DESeqDataSetFromMatrix(count_matrix_sub, colData = metadata_sub, design = ~ group)
# pre-filtering
smallestGroupSize <- 8 #only 8 samples on the smallest group (ETV6)
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
# indicate direction for comparison
dds$condition <- relevel(dds$group, ref = "ETV6_NTRK3_fused")
# normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
dim(normalized_counts) #25270    25
# Perform PCA
pca4 <- prcomp(t(normalized_counts)) 
aux <- as.data.frame(pca4$x[, 1:3]) 
# Merge data
scores4 <- merge(metadata_sub, aux, by.y=0, by.x="smpID", all.x=T)
# Load necessary library and set theme
library(ggplot2); theme_set(theme_classic())
library(ggrepel)
# Create PCA plot for all samples separated by articles
ggplot(scores4, aes(x=PC1, y=PC2, colour=factor(group), )) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'),name="Group") +
  xlab(paste0("PC1 (", prettyNum(summary(pca4)$importance[2,1]*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca4)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  #geom_text_repel(aes(label = smpID)) +
  ggtitle("PCA normalized + filtering genes - SJST033491_D1-S2(kinase)") 

# Save results
# saveRDS(scores4, file = "./PCA_normalized_filteredgenes_rm1outlier.rds")



# Building up metadata from Larissa's lab sheet --------------
# labsheet is from Larissa
# metadata is from the aligned samples
library(readxl)
labsheet <- read_excel("./metadata_lvf_sub.xlsx")
labsheet <- data.frame(labsheet)
labsheet <- labsheet[!labsheet$fastq_ID %in% "NA", ]
names(labsheet)[names(labsheet) == 'fastq_ID'] <- 'smpID'
library(stringr)
metadata$smpID <- str_split_fixed(as.character(metadata$smpID), "[-]", 2)[,1]
metadata <- merge(labsheet, metadata, by.x="smpID", all.x=T)
head(metadata)
metadata[, c('smpID', 'group')]
table(duplicated(metadata$smpID)) # 1 duplicated sample 
metadata[duplicated(metadata$smpID), ] # 1 duplicated sample: SJST031362_D1
# checking on count_matrix samples name
library(stringr)
colnames(count_matrix) <- str_split_fixed(as.character(colnames(count_matrix)), "[-]", 2)[,1]
all(metadata$smpID %in% colnames(count_matrix)) # FALSE
setdiff(metadata$smpID, colnames(count_matrix)) # "SJST030375_D1" "SJST032952_D1" "STST030375_R2"
setdiff(colnames(count_matrix), metadata$smpID ) # "SJST030375_R2"

metadata <- metadata[!metadata$smpID %in% c("SJST030375_D1","SJST032952_D1","STST030375_R2"), ]
count_matrix <- count_matrix[, !colnames(count_matrix) %in%  c("SJST030375_D1","SJST032952_D1","STST030375_R2")]
# saveRDS(metadata, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/metadata_labsheet.rds')

all(metadata$smpID %in% colnames(count_matrix)) # FALSE
count_matrix <- count_matrix[, colnames(count_matrix) %in% metadata$smpID]
all(metadata$smpID %in% colnames(count_matrix)) # TRUE

# PCA plot after filtering "SJST030375_D1" "SJST032952_D1" "STST030375_R2" samples
# 3# normalized + pre-filtering 
# Create a DESeqDataSet object

dds <- DESeqDataSetFromMatrix(count_matrix, colData = metadata, design = ~ group)
# pre-filtering
smallestGroupSize <- 8 #only 8 samples on the smallest group (ETV6)
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
# indicate direction for comparison
dds$condition <- relevel(dds$group, ref = "ETV6_NTRK3_fused")
# normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
dim(normalized_counts) #25280    25
# saveRDS(normalized_counts, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/run_on_terminal/normalized_counts.rds')

# Perform PCA
pca3 <- prcomp(t(normalized_counts)) 
aux <- as.data.frame(pca3$x[, 1:3]) 
# Merge data
scores5 <- merge(metadata, aux, by.y=0, by.x="smpID", all.x=T)

# Make a column for histologic diagnosis to be plotted over outilers samples only
scores5$histologic_for_outlier <- 'Others'
outlier_smp <- c('SJST033491_D1',
     'SJST032952_D2',
     'SJST032952_D4',
     'SJST032767_D2',
     'SJST033312_D1',
     'SJST033835_D1',
     'SJST031920_D1')
his_diag <- scores5[scores5$smpID %in% outlier_smp, ]$histologic_diagnosis

scores5[scores5$smpID %in% outlier_smp[1], ]$histologic_for_outlier <- his_diag[1]
scores5[scores5$smpID %in% outlier_smp[2], ]$histologic_for_outlier <- his_diag[2]
scores5[scores5$smpID %in% outlier_smp[3], ]$histologic_for_outlier <- his_diag[3]
scores5[scores5$smpID %in% outlier_smp[4], ]$histologic_for_outlier <- his_diag[4]
scores5[scores5$smpID %in% outlier_smp[5], ]$histologic_for_outlier <- his_diag[5]
scores5[scores5$smpID %in% outlier_smp[6], ]$histologic_for_outlier <- his_diag[6]
scores5[scores5$smpID %in% outlier_smp[7], ]$histologic_for_outlier <- his_diag[7]




# Load necessary library and set theme
library(ggplot2); theme_set(theme_classic())
library(ggrepel)
# Create PCA plot for all samples separated by articles
ggplot(scores5, aes(x=PC1, y=PC2, colour=factor(group), shape=histologic_for_outlier, )) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'),name="Group") +
  
  xlab(paste0("PC1 (", prettyNum(summary(pca3)$importance[2,1]*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca3)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  ggtitle("PCA normalized + filtering genes - smps out of metadata") 
  #geom_text_repel(aes(label = smpID))
# saveRDS(scores5, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/run_on_terminal/metadata.rds')

# Fresh/frozen (smp_type) seems a batch effect 


  
# Remove batch effect ------- We don't need it until we find a likely confound variable 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("bladderbatch")

# Usage
library(bladderbatch)
data(bladderdata)
dat <- bladderEset[1:50,]
pheno <- pData(dat)
edata <- exprs(dat)
batch <- pheno$batch
mod <- model.matrix(~as.factor(cancer), data = pheno)
library(sva)
ComBat(dat = edata, batch = batch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)

# my data (running on terminal)
# order colnames to scores5$smpID
normalized_counts <- normalized_counts[, scores5$smpID]
identical(colnames(normalized_counts), scores5$smpID) #TRUE
# prepare covariete and batch columns
batch <- scores5$smp_type #batch
mod <- model.matrix(~as.factor(group), data = scores5) #covariate 
library(sva)
start <- Sys.time()
# run combat
norm_combat <- ComBat(dat=normalized_counts[, ], batch=batch, mod=mod, par.prior = FALSE)
# saveRDS(norm_combat, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/norm_combat.rds')
end <- Sys.time()
print(end - start)

norm_combat <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/norm_combat.rds")

norm_combat_NOcov_par <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/norm_combat_NOcov_par.rds")

# Perform PCA
pca6 <- prcomp(t(norm_combat_NOcov_par)) 
aux <- as.data.frame(pca6$x[, 1:3]) 
# Merge data
scores6 <- merge(metadata, aux, by.y=0, by.x="smpID", all.x=T)

# Load necessary library and set theme
library(ggplot2); theme_set(theme_classic())
library(ggrepel)
# Create PCA plot for all samples separated by articles
ggplot(scores6, aes(x=PC1, y=PC2, colour=factor(group), )) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'),name="Group") +
  
  xlab(paste0("PC1 (", prettyNum(summary(pca6)$importance[2,1]*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca6)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  ggtitle("PCA combat-normalized + filtering genes - smps out of metadata") +
  geom_text_repel(aes(label = smpID))
# saveRDS(scores6, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/run_on_terminal/combat_NOcov_par.rds')




# Parametrics batch correction
# Remove the SJST033491_D1 outlier
norm_combat_NOcov_par <- norm_combat_NOcov_par[, !colnames(norm_combat_NOcov_par) %in% 'SJST033491_D1']
metadata[metadata$smpID %in% 'SJST033491_D1', ] # just to see it

# Perform PCA
pca7 <- prcomp(t(norm_combat_NOcov_par)) 
aux <- as.data.frame(pca7$x[, 1:3]) 
# Merge data
scores7 <- merge(metadata, aux, by.y=0, by.x="smpID", all.x=T)

# Load necessary library and set theme
library(ggplot2); theme_set(theme_classic())
library(ggrepel)
# Create PCA plot for all samples separated by articles
ggplot(scores7, aes(x=PC1, y=PC2, colour=factor(group), )) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'),name="Group") +
  
  xlab(paste0("PC1 (", prettyNum(summary(pca7)$importance[2,1]*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca7)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  ggtitle("PCA combat-normalized + filtering genes - smps out of metadata - outliers") 
  #geom_text_repel(aes(label = smpID))
# saveRDS(scores7, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/PCA_combat_NOcov_par.rds')

ggplot(scores7, aes(x=PC2, y=PC3, colour=factor(group), )) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'),name="Group") +
  
  xlab(paste0("PC2 (", prettyNum(summary(pca7)$importance[2,2]*100, digits = 2), "%)")) +
  ylab(paste0("PC3 (", prettyNum(summary(pca7)$importance[2,3]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  ggtitle("PCA combat-normalized + filtering genes - smps out of metadata - outliers") #+
  #geom_text_repel(aes(label = smpID))

dim(norm_combat_NOcov_par)
# Just let not use this one norm_combat_NOcov_par. It looks like a lot to the NONparametric over a PCA but it's parametric. So as Idk the assumptions we should consider it may be better keep with the NONparametric bacth effect correction - it's right below 



# NonParametrics batch correction
# Remove the SJST033491_D1 outlier
norm_combat_NOcov_nonpar <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/norm_combat_NOcov_nonpar.rds")
norm_combat_NOcov_nonpar <- norm_combat_NOcov_nonpar[, !colnames(norm_combat_NOcov_nonpar) %in% 'SJST033491_D1']
metadata[metadata$smpID %in% 'SJST033491_D1', ] # just to see it

# Perform PCA
pca8 <- prcomp(t(norm_combat_NOcov_nonpar)) 
aux <- as.data.frame(pca8$x[, 1:3]) 
# Merge data
scores8 <- merge(metadata, aux, by.y=0, by.x="smpID", all.x=T)

# Load necessary library and set theme
library(ggplot2); theme_set(theme_classic())
library(ggrepel)
# Create PCA plot for all samples separated by articles
ggplot(scores8, aes(x=PC1, y=PC2, colour=factor(group), )) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'),name="Group") +
  
  xlab(paste0("PC1 (", prettyNum(summary(pca8)$importance[2,1]*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca8)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  ggtitle("PCA combat-normalized + filtering genes - smps out of metadata - outliers") 
#geom_text_repel(aes(label = smpID))


dim(norm_combat_NOcov_nonpar)
DE_input_matrix <- norm_combat_NOcov_nonpar
# saveRDS(DE_input_matrix, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/DE_input_matrix.rds')



# norm_combat_Cov_par batch correction (I've visualized)
# norm_combat_Cov_nonpar batch correction (I didn't view this)

# Conclusion: The Covariate mode didin't show any difference. So that is it. We're moving on with the norm_combat_NOcov_nonpar as input for the Differential Expression. 

