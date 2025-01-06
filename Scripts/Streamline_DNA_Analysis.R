### Description: from DNAmet matrix (beta-value matrix) to downsteram analysis 
### Technology: 850k DNAmet array

library(BiocManager)
require(IlluminaHumanMethylationEPICmanifest)
require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) 
library(shinyMethyl)
library(bumphunter)
library(minfi)
library(dplyr)


## Most Variable CpGs  --------------------
load('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/beta_and_meta_data.rda')
beta %>% dim() #856801     38
beta <- as.data.frame(beta)
beta <- na.omit(beta)
beta %>% dim() #836811     38
dim(beta_meta)
head(beta_meta)
names(beta_meta)[names(beta_meta) == "sample"] <- "smp_ID"

df <- data.frame(
  smp_ID = c("207558610121_R01C01", NA, "205715820168_R06C01", "207495040012_R05C01", NA, "206467000107_R05C01", "207558610121_R02C01", "207558610121_R03C01", "207558610121_R04C01", "207558610121_R05C01", NA, "207558610121_R06C01", "207558610121_R07C01", "207179240158_R05C01", NA, "207558610121_R08C01", "207179240091_R05C01", "207343240018_R08C01", "206154070072_R03C01", "207558610086_R02C01", "206250280180_R02C01", NA, "207339140028_R01C01", "207558610127_R01C01", "207558610127_R02C01", "207558610127_R03C01", "207558610127_R04C01", "207558610127_R05C01", "207558610127_R06C01", "207558610127_R07C01", "201869680190_R01C01", "201332340140_R08C01", "200848860136_R01C01", "200848860134_R08C01", "200928190007_R08C01", "200928190007_R07C01", "200928190007_R06C01", "200928190007_R04C01", "200928190007_R03C01", "200928190007_R02C01", "200928190007_R01C01", "200925700210_R06C01", "200848860099_R05C01", "9553932008_R01C02"),
  
  fusion = c("NRF1::BRAF", "PLEKHH2::ALK", "PLEKHH2::ALK", "TPR::NTRK1", "TPM3::NTRK1", "TPM3::NTRK1", "TPR::NTRK1", "LMNA::NTRK1", "LMNA::NTRK1", "PDE4DIP::NTRK1", "PDE4DIP::NTRK1", "EML4::NTRK3", "RBPMS::NTRK3", "RCC1::ALK", "EML4::NTRK3", "TTYH3::BRAF", "TPM3::NTRK1", "EML4::ALK\nOutside lab", "TPM3::NTRK1\nOutside lab", "TPM3::NTRK1", "ETV6:NTRK3", "ETV6:NTRK3", "ETV6:NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  
  histologic = c("Spindle cell neoplasm with NRF1::BRAF rearrangement", "Spindle cell neoplasm", "Spindle cell neoplasm", "Spindle cell neoplasm with evidence of TPR::NTRK1 rearrangement", "Spindle cell sarcoma with NTRK1 gene rearrangement", "NTRK-rearranged spindle cell neoplasm", "malignant spindle cell neoplasm with TPR::NTRK1 fusion, most in keeping with infantile fibrosarcoma", "low-grade sarcoma with LMNA-NTRK1 fusion and homozygous deletion of CDKN2A,", "Lipofibromatosis", "Recurrent/residual variant of infantile fibrosarcoma", "Infantile fibrosarcoma, variant with NTRK1 fusion, status post chemotherapy", "Infantile fibrosarcoma-like tumor", "Congenital NTRK-rearranged spindle cell neoplasm of the gastrointestinal tract with evidence of RBPMS::NTRK3 fusion transcript", "Spindle cell neoplasm with ALK gene rearrangement", "Low-grade fibroblastic neoplasm with neural differentiation suggestive of NTRK1-associated mesenchymal tumor.", "suspicious for infantile fibrosarcoma", "Spindle cell neoplasm with TPM3::NTRK1 fusion, most c/w infantile fibrosarcoma", "Spindle cell neoplasm with EML4::ALK fusion", "Spindled cell neoplasm with TPM3-NTRK1 fusion and meningioangiomatosis-like growth pattern;", "NTRK-rearranged sarcoma with evidence of TPM3::NTRK1 fusion", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  
  tumor_type = c("kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor"),
  
  smp_type = c("FFPE", "FFPE", "Fresh/frozen", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", "Fresh/frozen", "Fresh/frozen", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", NA, NA, "FFPE", "Fresh/frozen", "Fresh/frozen", "Fresh/frozen", "Fresh/frozen", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  
  stringsAsFactors = FALSE
)

beta_meta_st <- merge(beta_meta, df, by="smp_ID")
beta_meta_st[beta_meta_st$fusion %in% "ETV6:NTRK3", ]$fusion <- "ETV6::NTRK3"
beta_meta_st[beta_meta_st$fusion %in% "TPM3::NTRK1\nOutside lab", ]$fusion <- "TPM3::NTRK1"
beta_meta_st[beta_meta_st$fusion %in% "EML4::ALK\nOutside lab", ]$fusion <- "EML4::ALK"
beta_meta_GEO <- beta_meta[!beta_meta$smp_ID %in% beta_meta_st$smp_ID, ]
dim(beta_meta_GEO)#13 3
beta_meta <- plyr::rbind.fill(beta_meta_st, beta_meta_GEO)
dim(beta_meta) #38  7
head(beta_meta)


library(pheatmap)
library(RColorBrewer)
dim(beta) #
dim(beta_meta)#
my_sample_col <- data.frame(row.names = beta_meta$smp_ID, 
                            Group = beta_meta$group, 
                            #Histology = labsheet_hm$histologic_diagnosis,
                            Fusion = beta_meta$fusion) 

NTRK <- my_sample_col[grep('NTRK', my_sample_col$Fusion), ]
NTRK_butETV6 <- rownames(NTRK[!NTRK$Fusion %in% c('ETV6::NTRK3'), ])
NTRK_onlyETV6 <- rownames(NTRK[!rownames(NTRK) %in% NTRK_butETV6, ])

my_sample_col$Fusion_concat <- NA
my_sample_col[rownames(my_sample_col) %in% NTRK_butETV6, ]$Fusion_concat <- 'Others_NTRK'
my_sample_col[rownames(my_sample_col) %in% NTRK_onlyETV6, ]$Fusion_concat <- 'ETV6_NTRK3'

my_sample_col[my_sample_col$Fusion_concat %in% NA, ]$Fusion_concat[1:5] <- "Not_NTRK" 

my_sample_col[my_sample_col$Fusion_concat %in% NA, ]$Fusion_concat <- "NA"

my_sample_col[my_sample_col$Fusion %in% NA, ]$Fusion <- "NA"

my_sample_col$Source <- NA
my_sample_col[my_sample_col$Fusion %in% "NA", ]$Source <- "PublicData" 
my_sample_col[!my_sample_col$Fusion %in% "NA", ]$Source <- "StJude" 


ann_colors = list(
  Fusion_concat = c(ETV6_NTRK3 = 'black', Not_NTRK = 'blue', Others_NTRK = 'red', `NA` = "gray"))


# Most variable feature - from PC1 components 
pca_data <- na.omit(beta)
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

# Correlation test
beta_sub <- beta[MVC_PC1_2_3, ]
library(ggplot2)
library(corrplot)
length(colnames(beta_sub))
length(rownames(beta_meta))
intersect(colnames(beta_sub),
          rownames(beta_meta))
identical(colnames(beta_sub),
          rownames(beta_meta)) #FALSE
rownames(beta_meta) <- beta_meta$smp_ID
beta_meta <- beta_meta[colnames(beta_sub), ]
identical(colnames(beta_sub),
          rownames(beta_meta)) #TRUE

cor_test_mat <- cor.mtest(beta_sub, conf.level = 0.95)
cor_test_mat <- cor_test_mat$p
cor_matrix <- cor(beta_sub)

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
         main = "Most Variable CpGs from PC1,2,3 - Correlation - scale:none")


## Without displying pvalue_symbol_matrix

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Plots/DNAmet_cor_hm.png",  width = 8, height = 8, units = "in", res = 300)

p <- pheatmap(cor_matrix[,], #smp_ID on columns and rows
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
         main = "Most Variable CpGs from PC1,2,3 - Correlation - scale:none")

print(p)

dev.off()


## PCA and tSNE --------------
# Stjude + PublicData 
load('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/beta_and_meta_data.rda')
beta[1:3,1:3]
beta_meta[1:3,]
names(beta_meta)[names(beta_meta) == "sample"] <- "smp_ID"

df <- data.frame(
  smp_ID = c("207558610121_R01C01", NA, "205715820168_R06C01", "207495040012_R05C01", NA, "206467000107_R05C01", "207558610121_R02C01", "207558610121_R03C01", "207558610121_R04C01", "207558610121_R05C01", NA, "207558610121_R06C01", "207558610121_R07C01", "207179240158_R05C01", NA, "207558610121_R08C01", "207179240091_R05C01", "207343240018_R08C01", "206154070072_R03C01", "207558610086_R02C01", "206250280180_R02C01", NA, "207339140028_R01C01", "207558610127_R01C01", "207558610127_R02C01", "207558610127_R03C01", "207558610127_R04C01", "207558610127_R05C01", "207558610127_R06C01", "207558610127_R07C01", "201869680190_R01C01", "201332340140_R08C01", "200848860136_R01C01", "200848860134_R08C01", "200928190007_R08C01", "200928190007_R07C01", "200928190007_R06C01", "200928190007_R04C01", "200928190007_R03C01", "200928190007_R02C01", "200928190007_R01C01", "200925700210_R06C01", "200848860099_R05C01", "9553932008_R01C02"),
  
  fusion = c("NRF1::BRAF", "PLEKHH2::ALK", "PLEKHH2::ALK", "TPR::NTRK1", "TPM3::NTRK1", "TPM3::NTRK1", "TPR::NTRK1", "LMNA::NTRK1", "LMNA::NTRK1", "PDE4DIP::NTRK1", "PDE4DIP::NTRK1", "EML4::NTRK3", "RBPMS::NTRK3", "RCC1::ALK", "EML4::NTRK3", "TTYH3::BRAF", "TPM3::NTRK1", "EML4::ALK\nOutside lab", "TPM3::NTRK1\nOutside lab", "TPM3::NTRK1", "ETV6:NTRK3", "ETV6:NTRK3", "ETV6:NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  
  histologic = c("Spindle cell neoplasm with NRF1::BRAF rearrangement", "Spindle cell neoplasm", "Spindle cell neoplasm", "Spindle cell neoplasm with evidence of TPR::NTRK1 rearrangement", "Spindle cell sarcoma with NTRK1 gene rearrangement", "NTRK-rearranged spindle cell neoplasm", "malignant spindle cell neoplasm with TPR::NTRK1 fusion, most in keeping with infantile fibrosarcoma", "low-grade sarcoma with LMNA-NTRK1 fusion and homozygous deletion of CDKN2A,", "Lipofibromatosis", "Recurrent/residual variant of infantile fibrosarcoma", "Infantile fibrosarcoma, variant with NTRK1 fusion, status post chemotherapy", "Infantile fibrosarcoma-like tumor", "Congenital NTRK-rearranged spindle cell neoplasm of the gastrointestinal tract with evidence of RBPMS::NTRK3 fusion transcript", "Spindle cell neoplasm with ALK gene rearrangement", "Low-grade fibroblastic neoplasm with neural differentiation suggestive of NTRK1-associated mesenchymal tumor.", "suspicious for infantile fibrosarcoma", "Spindle cell neoplasm with TPM3::NTRK1 fusion, most c/w infantile fibrosarcoma", "Spindle cell neoplasm with EML4::ALK fusion", "Spindled cell neoplasm with TPM3-NTRK1 fusion and meningioangiomatosis-like growth pattern;", "NTRK-rearranged sarcoma with evidence of TPM3::NTRK1 fusion", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  
  tumor_type = c("kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor"),
  
  smp_type = c("FFPE", "FFPE", "Fresh/frozen", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", "Fresh/frozen", "Fresh/frozen", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", NA, NA, "FFPE", "Fresh/frozen", "Fresh/frozen", "Fresh/frozen", "Fresh/frozen", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  
  stringsAsFactors = FALSE
)

# Review $fusion labels 
beta_meta_st <- merge(beta_meta, df, by="smp_ID")
beta_meta_st[beta_meta_st$fusion %in% "ETV6:NTRK3", ]$fusion <- "ETV6::NTRK3"
beta_meta_st[beta_meta_st$fusion %in% "TPM3::NTRK1\nOutside lab", ]$fusion <- "TPM3::NTRK1"
beta_meta_st[beta_meta_st$fusion %in% "EML4::ALK\nOutside lab", ]$fusion <- "EML4::ALK"
beta_meta_GEO <- beta_meta[!beta_meta$smp_ID %in% beta_meta_st$smp_ID, ]
dim(beta_meta_GEO)#13 3
beta_meta <- plyr::rbind.fill(beta_meta_st, beta_meta_GEO)
dim(beta_meta) #38  7



# Run-plot PCA
beta <- na.omit(beta)
pca <- prcomp(t(beta)) 
aux <- as.data.frame(pca$x[, 1:3]) 
scores <- merge(beta_meta, aux, by.y=0, by.x="smp_ID", all.x=T)

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Plots/DNAmet_PCA_supplementary.png",  width = 8, height = 6, units = "in", res = 300)

library(ggplot2); theme_set(theme_classic())
p <- ggplot(scores, aes(x=PC1, y=PC2, colour=factor(group), 
                   #shape = smp_type
)) +
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
  ggtitle("Whole Array - StJude + PublicData")

print(p)

dev.off()



# Run-plot tSNE (V1)
library(Rtsne)
dim(beta)
set.seed(234234)
tsne_realData <- Rtsne(t(beta), perplexity=3, check_duplicates = FALSE) # #function to run t-sne
pdata.teste.tsne <- scores #renaming metadata
pdata.teste.tsne$tSNE1 <- tsne_realData$Y[,1]
pdata.teste.tsne$tSNE2 <- tsne_realData$Y[,2]

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Plots/DNAmet_tSNE.png",  width = 8, height = 6, units = "in", res = 300)

library(ggplot2); theme_set(theme_classic())
p <- ggplot(pdata.teste.tsne, aes(x=tSNE1, y=tSNE2, 
                             colour=factor(group), 
                             #shape = smp_type
)) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'), name="Group") +
  # xlab(paste0("PC1 (", prettyNum(summary(pca)$importance[2,1]*100, digits = 2), "%)")) +
  # ylab(paste0("PC2 (", prettyNum(summary(pca)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  #geom_text_repel(aes(label = smpID)) +
  ggtitle("Whole Array - StJude + PublicData") 

print(p)
dev.off()


# # Plot tSNE (V2)
# library(ggplot2); theme_set(theme_classic())
# ggplot(pdata.teste.tsne, aes(x=tSNE1, y=tSNE2, 
#                              colour=factor(group), 
#                              shape = source
# )) +
#   geom_point(size = 4) +
#   scale_color_manual(values=c('slateblue3', 'wheat3'), name="Group") +
#   # xlab(paste0("PC1 (", prettyNum(summary(pca)$importance[2,1]*100, digits = 2), "%)")) +
#   # ylab(paste0("PC2 (", prettyNum(summary(pca)$importance[2,2]*100, digits = 2), "%)")) +
#   scale_x_continuous(labels = scales::scientific_format()) +
#   scale_y_continuous(labels = scales::scientific_format()) +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14, face = "bold"),
#         legend.title = element_text(size = 16),  
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, color = "black", face = "bold")) +
#   #geom_text_repel(aes(label = smpID)) +
#   ggtitle("Whole Array - StJude + PublicData") 






## Differential Methylated Position (DMP) --------------
load('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/beta_and_meta_data.rda')
beta %>% dim() #856801     38
beta <- as.data.frame(beta)
beta <- beta[, colnames(beta) %in% beta_meta$sample]
identical(colnames(beta), beta_meta$sample) #TRUE
beta <- na.omit(beta)
dim(beta) #836811     38
# Run dmpFinder (method 1) -----
dmp_allspm <- dmpFinder(as.matrix(beta), pheno = beta_meta$group  , type = "categorical")
table(dmp_allspm$qval < 0.05) #99234 CpGs statistic significant  
dmp_sig <- dmp_allspm[dmp_allspm$qval < 0.05, ]
dmp_sig$CpG <- rownames(dmp_sig)
dmp_sig_method1 <- dmp_sig

# Wilcoxon test (method 2) ---------
identical(colnames(beta), beta_meta$sample) #TRUE
# Ordering the samples by group of comparison in beta
etv6 <- beta_meta[beta_meta$group %in% 'ETV6_NTRK3_fused', ]$sample
kinase <- beta_meta[beta_meta$group %in% 'Kinase_fused', ]$sample
etv6_mt <- beta[, colnames(beta) %in% etv6]
length(colnames(etv6_mt)) #22 samples
kinase_mt <- beta[, colnames(beta) %in% kinase]
length(colnames(kinase_mt)) #16 samples
beta <- cbind(kinase_mt, etv6_mt)

# Run Differential Methylation 
require(exactRankTests)
require(parallel)
values <- t(beta) #transpose beta
values <- data.frame(values)
# parallel processing ON for wilcoxon test
wpvalues <- unlist(mclapply(values,
                            function(CpG) {
                              zz <- wilcox.exact(CpG[1:  dim(kinase_mt)[2]],
                                                 CpG[c(dim(kinase_mt)[2]+1) : dim(beta)[2]], exact=T) 
                              z <- zz$p.value
                              return(z)
                            }, mc.cores= 20)) # set n of cores
# adjust pvalue 
wpvalues_adj <- p.adjust(wpvalues, method = "BH")
wpvalues_adj <- data.frame(wpvalues_adj)
wpvalues_adj$CpG <- rownames(wpvalues_adj)
hist(wpvalues_adj$wpvalues_adj)
table(wpvalues_adj$wpvalues_adj < 0.05) #75795 CpGs statistic significant 
wpvalues_adj_sig <- wpvalues_adj[wpvalues_adj$wpvalues_adj < 0.05, ]
dmp_sig_method2 <- wpvalues_adj_sig

# Merge both DMP results
DMP_CpGs <- merge(dmp_sig_method1, dmp_sig_method2, by = 'CpG') 
dim(DMP_CpGs) #71094 CpGs statistic significant 
# saveRDS(dmp_sig_merged, '/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/dmp_sig_merged.rds')

# It takes a long time to get DMP_CpGs. Start from it's load if you already ran it.


# Calculate Mean Differential DNA methylation ----------
DMP_CpGs <- readRDS('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/dmp_sig_merged.rds')
load('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/beta_and_meta_data.rda')

calculate_diff_mean <- function(data, 
                                metadata, 
                                condition_one, 
                                threshold) {
  
  list_diffmean_dfs <- list()
  data <- as.data.frame(data)
  condition_all_but_one <- group_names[!group_names %in% condition_one]
  ###condition_all_but_one <- 'Mesenchymal-like' #TESTING only. 
  for(i in 1:length(condition_all_but_one)) {
    print(paste0("Processing group: ", condition_all_but_one[i]))
    
    try({ #it ignores any error that might prevent the code of going forward
      #I did it because of the absence of hyper or hypo probes were making the loop stop
      loop_data <- data
      condition_1_sampleID <- metadata[metadata$group %in% condition_one, ]$sample
      condition_2_sampleID <- metadata[metadata$group %in% condition_all_but_one[i], ]$sample
      
      print(paste0("Length of condition_1_sampleID: ", length(condition_1_sampleID)))
      print(paste0("Length of condition_2_sampleID: ", length(condition_2_sampleID)))
      
      loop_data$meanM1 <- apply(loop_data[, condition_1_sampleID], 1, mean, na.rm = TRUE)
      loop_data$meanM2 <- apply(loop_data[, condition_2_sampleID], 1, mean, na.rm = TRUE)
      loop_data$DiffMean <- loop_data$meanM1 - loop_data$meanM2
      loop_data$Comparison <- paste0(condition_one, '_', condition_all_but_one[i], '_', condition_one, '_', 'Orientation')
      
      
      
      loop_data$DNAmet_orientation <- NA
      #loop_data$DNAmet_orientation <- as.character(loop_data$DNAmet_orientation) #this time it's not necessary
      if(dim(loop_data[loop_data$DiffMean > threshold, ])[1] > 0) {
        loop_data[loop_data$DiffMean > threshold, ]$DNAmet_orientation <- 'hyper'
      } else {
        # do nothing
      }
      
      if(dim(loop_data[loop_data$DiffMean < -threshold, ])[1] > 0) {
        loop_data[loop_data$DiffMean < -threshold, ]$DNAmet_orientation <- 'hypo'
      } else {
        # do nothing
      }
      
      loop_data[loop_data$DNAmet_orientation %in% NA, ]$DNAmet_orientation <- 'not_diff'
      
      
      
      loop_data$probeID <- rownames(loop_data)
      
      list_diffmean_dfs[[i]] <- loop_data[, c('DiffMean', 'Comparison', 'DNAmet_orientation', 'probeID')]
      print(paste0(list_diffmean_dfs[[i]]$Comparison[1], ' has been stored into the list.'))
    },  silent = FALSE)
  }
  return(list_diffmean_dfs) 
}

# Filtering beta to only sig. CpGs 
beta <- beta[rownames(beta) %in% DMP_CpGs$CpG, ]
dim(beta) #71094    38
group_names <- names(table(beta_meta$group))
list_diffmean_dfs <- calculate_diff_mean(data = beta,
                                         metadata = beta_meta,
                                         condition_one = group_names[2], #"kinese_fused_tumors" is the direction 
                                         threshold = 0.3)

dmp_diffmean <- do.call('rbind', list_diffmean_dfs)
head(dmp_diffmean)
length(names(table(dmp_diffmean$Comparison))) # it should be equal to the number if groups you are comparing - 1

names(DMP_CpGs)[names(DMP_CpGs) == "CpG"] <- "probeID"
DMP_CpGs <- DMP_CpGs[, c('probeID', 'wpvalues_adj')]
dmp_diffmean <- merge(dmp_diffmean, DMP_CpGs, by = 'probeID')
# saveRDS(dmp_diffmean_filtered, '/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/dmp_diffmean.rds')


## Plot Volcano 
load('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/beta_and_meta_data.rda')
volcano_plot_input <- dmp_diffmean

## Volcano V1
library(EnhancedVolcano) # ISSUE: somehow I can't add gene names onto volcano plot using EnhancedVolcano
rownames(volcano_plot_input) <- volcano_plot_input$probeID # don't run it if it doesn't need to

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Plots/DNAmet_DMP_Volcano.png",  width = 10, height = 6, units = "in", res = 300)

p <- EnhancedVolcano(
  volcano_plot_input, # DMP
  #lab = rownames(res_kinaseDirection),
  #lab = as.character(rownames(res_kinaseDirection)),
  lab = '',
  x = 'DiffMean',
  y = 'wpvalues_adj',
  #selectLab = selecected_features,
  pCutoff = 10e-2, # p-value cutoff line
  FCcutoff = 0.30, # fold change cutoff line
  xlab = paste0("<---- ", "ETV6-NTRK3", "   ", "Mean Methylation Difference", "   ", "Kinase", " ---->"),
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
  title = paste0("DMPs ", "Kinase", " vs ", "ETV6-NTRK3"), 
  subtitle = "", # Subtitle should be on the same line as title
  ylim = c(0, 5),
  xlim = c(-0.5, 0.5)
)

print(p)
dev.off()

## Volcano V2
# dim(beta) # 856801     38
# beta <- as.data.frame(beta)
# beta$probeID <- rownames(beta)
# volcano <- merge(dmp_diffmean, beta, by = "probeID")
# dim(volcano)
# table(volcano$DNAmet_orientation)
# # hyper   hypo  not_diff 
# # 75      580    70439
# 
# 
# volcano$Statistical_Sig <- "not_stat_sig"
# volcano[volcano$wpvalues_adj <= 0.05 &
#           volcano$DNAmet_orientation %in% c("hypo"), ]$Statistical_Sig <- "hypo"
# volcano[volcano$wpvalues_adj <= 0.05 &
#           volcano$DNAmet_orientation %in% c("hyper"), ]$Statistical_Sig <- "hyper"
# 
# 
# library(ggplot2); theme_set(theme_classic()) 
# ggplot(data=volcano[, ], aes(x=DiffMean, y=-1*log10(wpvalues_adj), colour=Statistical_Sig)) +
#   geom_point() +
#   xlab("Diff Mean Methylation") + ylab("-1 * log10 Significance") + 
#   scale_color_manual(breaks=c("hyper","hypo","not_stat_sig"), # color scale (for points) 
#                      values=c("red","blue","black"),
#                      labels=c('Hyper methylated',"Hypo methylated","< |0.3| Diff Mean"),
#                      name="Legend")  +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14, face = "bold"),
#         legend.title = element_text(size = 20),  
#         legend.text = element_text(size = 18),
#         plot.title = element_text(size = 16, color = "black", face = "bold")) +
#   ggtitle("Kinase_fused vs ETV6_NTRK3_fused")




## Differential Methylated Region (DMR) --------------
library(dplyr)
library(minfi)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DMRcate")
library(DMRcate)
library(Gviz)
library(RColorBrewer)
library(GenomicRanges)
library(rtracklayer)
library(HelpersMG)
library(data.table)


# Generate DMRs table (DNAmet orientation, CpGs regions, n of CpGs, CpGs probe_ID etc) 
load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/beta_and_meta_data.rda")
head(beta)
dim(beta) #856801     38
table(is.na(beta))
# FALSE     TRUE 
# 32506123    52315
beta <- na.omit(beta)
dim(beta) #36811     38
head(beta_meta)
dim(beta_meta) #38  3
rownames(beta_meta) <- beta_meta$sample

identical(rownames(beta_meta), colnames(beta)) #TRUE
# write.csv(beta_meta, "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/objects/objects_from_to_rmd/metadata_DNAmet_sample_ID.csv")

# Prepare metadata as "design" for the DMR
beta_meta$group <- as.factor(beta_meta$group)
beta_meta$group <- relevel(beta_meta$group, ref = "ETV6_NTRK3_fused") # this set Kinase_fused as the direction of DNAmeth  
design <- model.matrix(~ group, data = beta_meta)
# setting some annotation
myAnnotation <- cpg.annotate(object = beta, datatype = "array", 
                             what = "Beta", 
                             analysis.type = "differential", 
                             design = design, 
                             contrasts = FALSE, 
                             coef = 2, # groupKinase_fused design column - the orientation we want to see (same orientation we set for RNAseq)
                             arraytype = "EPIC",
                             fdr = 0.001)
str(myAnnotation)

# DMR analysis
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2) #pcutoff = defaut
results.ranges <- extractRanges(DMRs, genome = 'hg38')
print(results.ranges) # DMR results 
dmr.table <- as.data.frame(results.ranges)
dim(dmr.table) #572  13 
### Saving
# saveRDS(dmr.table, '/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/reviewed_DM_direction/dmr_table.rds')
# dmr.table <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/reviewed_DM_direction/dmr_table.rds")

## Retrieving CpGs from DMRs
# CpGs to promoters 
### Retrieving CpGs from DMRs ----------
# Turn illumina manifest (CpG-gene annotation) into GRange object
# So we can get info. from indiviudal CpGs and gene annotations from the manifest
library(readr)
EPIC.hg38.anno <- read_csv("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip = 7)
# Make it a genomic range object
EPIC.hg38.anno.gr <- makeGRangesFromDataFrame(EPIC.hg38.anno, keep.extra.columns=T, start.field = "Start_hg38", end.field = "End_hg38", seqnames.field = "CHR_hg38", strand.field="Strand_hg38", na.rm=TRUE) 

# Components from DMR 
chr <- dmr.table$seqnames
start <- dmr.table$start
end <- dmr.table$end
# Merging illumina manifest annotation to DMR output
filtered_gr <- subsetByOverlaps(EPIC.hg38.anno.gr, GRanges(seqnames = chr, ranges = IRanges(start = start, end = end), keep.extra.columns = TRUE))
df <- as.data.frame(filtered_gr)
length(df$Name) #426 CpGs from DMR output
df$UCSC_RefGene_Name
df_CpGs_DMR <- data.frame(df$Name)

#Saving 
# saveRDS(df_CpGs_DMR, '/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/df_CpGs_DMR.rds')


df_CpGs_DMR <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/reviewed_DM_direction/df_CpGs_DMR.rds")
names(df_CpGs_DMR) <- "CpG_ID"

dmr_table <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/reviewed_DM_direction/dmr_table.rds")
# transforming both back in a Genomic Range obj to retrieve all dmr_table columns but with the specific CpGs in each genomic interval 
library(readr)
EPIC.hg38.anno <- read_csv("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip = 7)
EPIC.hg38.anno <- as.data.frame(EPIC.hg38.anno)
df_CpGs_DMR <- EPIC.hg38.anno[EPIC.hg38.anno$Name %in% df_CpGs_DMR$CpG_ID, ]


# Reference Genome Annotation (illumina manifest annotation)
EPIC.hg38.anno.gr <- makeGRangesFromDataFrame(EPIC.hg38.anno, keep.extra.columns=T, start.field = "Start_hg38", end.field = "End_hg38", seqnames.field = "CHR_hg38", strand.field="Strand_hg38", na.rm=TRUE) 

# DMR anno output 
dmr_table.anno.gr <- makeGRangesFromDataFrame(dmr_table, keep.extra.columns=T, start.field = "start", end.field = "end", seqnames.field = "seqnames", strand.field="strand", na.rm=TRUE)

## Retrieving DMR CpGs - because in DMR ouput we don't have "probe_ID" information + Mapping CpGs top promoter regions 
library(ELMER)
promoter.gr <- get.feature.probe(promoter = TRUE, TSS.range = list(upstream = 2000, downstream = 500), rm.chr=c("chrX", "chrY"), met.platform = "EPIC")

# Interact to one genomic interval at a time
df_CpG_list <- list()
for(i in 1:length(dmr_table.anno.gr$no.cpgs)) {
  dmr_table_gr_1 <- dmr_table.anno.gr[i, ]
  # Merging them both - illumina manifest annotation to DMR output
  filtered_gr <- subsetByOverlaps(promoter.gr, dmr_table_gr_1)
  df_CpG <- as.data.frame(filtered_gr)
  if(nrow(df_CpG) != 0){ 
    df_dmr <- as.data.frame(dmr_table_gr_1)
    df_CpG$maxdiff <- NA
    df_CpG$meandiff <- NA
    df_CpG$meth_status <- NA
    df_CpG$maxdiff <- df_dmr$maxdiff
    df_CpG$meandiff <- df_dmr$meandiff
    
    table(df_CpG$meandiff > 0)
    if(any(df_CpG$meandiff > 0)){
      df_CpG[df_CpG$meandiff > 0,]$meth_status <- "hyper"
    }
    if(any(df_CpG$meandiff < 0)){
      df_CpG[df_CpG$meandiff < 0,]$meth_status <- "hypo"
    }
    df_CpG_list[[i]] <- df_CpG
  } 
  
}

df_CpG_promoter <- do.call(rbind, df_CpG_list) # basic way to get DMR CpGs + meth_status + promoter regions  


# Hyper meth CpGs (Kinase hyper methylated compared to ETV6) 
hyper_CpG_promoter_df <- df_CpG_promoter[df_CpG_promoter$meth_status %in% "hyper", ]

df <- table(hyper_CpG_promoter_df$gene) %>% data.frame()
df <- df[order(df$Freq, decreasing = TRUE), ]

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Plots/DNAmet_hist_Kinase_Hyper_CpGfreq_to_promoter.png",  width = 8, height = 8, units = "in", res = 300)

library(ggplot2); theme_set(theme_classic())
p <- ggplot(df[, ], aes(reorder(x = factor(Var1), Freq), y = Freq, )) +
  geom_bar(stat = "identity", color = "black") +
  labs(#title = "Frequency of Cells by Sample ID",
    x = "Genes",
    y = "Frequency CpG on putative promoters") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  labs(fill = "") +
  #facet_wrap(~ factor(Var1), ncol = 3) +
  
  ggtitle("Freq of Kinase HYPER CpG from DMR and the promoter-genes annotation") 

print(p)
dev.off()

# Hypo meth CpGs (Kinase hypo methylated compared to ETV6) 
hypo_CpG_promoter_df <- df_CpG_promoter[df_CpG_promoter$meth_status %in% "hypo", ]

df <- table(hypo_CpG_promoter_df$gene) %>% data.frame()
df <- df[order(df$Freq, decreasing = TRUE), ]

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Plots/DNAmet_hist_Kinase_Hypo_CpGfreq_to_promoter.png",  width = 8, height = 8, units = "in", res = 300)

library(ggplot2); theme_set(theme_classic())
p <- ggplot(df[, ], aes(reorder(x = factor(Var1), Freq), y = Freq, )) +
  geom_bar(stat = "identity", color = "black") +
  labs(#title = "Frequency of Cells by Sample ID",
    x = "Genes",
    y = "Frequency CpG on putative promoters") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  labs(fill = "") +
  #facet_wrap(~ factor(Var1), ncol = 3) +
  
  ggtitle("Freq of Kinase HYPO CpG from DMR and the promoter-genes annotation") 

print(p)
dev.off()


## Stemness score ------------
load('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/beta_and_meta_data.rda')
names(beta_meta)[names(beta_meta) == "sample"] <- "smp_ID"
load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/objects/DNAmet_stemness_model.Rda") 
w <- mm$w
w[1:5]
w_df <- as.data.frame(w)
w_df$probeID <- rownames(w_df)
X <- as.data.frame(beta)
X <- na.omit(X) 
# Make DNA mtx (X) and the stemness weigth (w_df) be the same dim and in the same order
X <- X[rownames(X) %in%  rownames(w_df) ,]
w_df <- w_df[rownames(X), ]
length(intersect(rownames(X),rownames(w_df))) # 206
identical(rownames(X),rownames(w_df)) # TRUE
X <- as.matrix(X)
dim(X) # 208  38
dim(w_df) #208   2
length(intersect(rownames(w_df),rownames(X))) # 206
identical(rownames(w_df),rownames(X)) #TRUE
w_df$probeID <- NULL 
# Apply stemness model into DNA mtx 
ss <- t(w_df) %*% X 
ss[1,1:3]
# Scale the scores into a ratio from 0 to 1 and store as data frame.
ss <- ss - min(ss)
ss <- ss / max(ss)
ss <- as.data.frame(t(ss))
colnames(ss) <- "Stemness_DNAmet" 
ss$smp_ID <- rownames(ss)

# Plot PCA 
beta_meta <- merge(beta_meta, ss, by="smp_ID")
beta <- na.omit(beta)
pca <- prcomp(t(beta)) 
aux <- as.data.frame(pca$x[, 1:3]) 
scores <- merge(beta_meta, aux, by.y=0, by.x="smp_ID", all.x=T)

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Plots/DNAmet_PCA_Stemness.png",  width = 8, height = 6, units = "in", res = 300)

library(ggplot2); theme_set(theme_classic())
p <- ggplot(scores, aes(PC1, PC2)) + geom_point(size=5, aes( fill=Stemness_DNAmet, color = group), pch = 21, stroke = 1.5) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_color_manual(values = c('slateblue3', 'wheat3')) +
  ylab(paste0('PC2 ', summary(pca)$importance[2, 2] * 100, '%')) +
  xlab(paste0('PC1 ', summary(pca)$importance[2, 1] * 100, '%') ) + theme_bw() +
  ggtitle('WholeArray - Stemness prediction')

print(p)
dev.off()


# Plot Boxplot 
# Statistical tests
t.test(scores[scores$group %in% "ETV6_NTRK3_fused", ]$Stemness_DNAmet,
       scores[scores$group %in% "Kinase_fused", ]$Stemness_DNAmet) # p-value = 0.01745
library(rstatix)
scores %>%
  wilcox_test(Stemness_DNAmet ~ group, p.adjust.method = "none") %>%
  add_significance() #p-value = 0.0224


png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Plots/DNAmet_Boxplot_Stemness.png",  width = 8, height = 6, units = "in", res = 300)

library(ggplot2); theme_set(theme_classic())
p <- ggplot(scores[, ], aes(x=group, y=Stemness_DNAmet)) + 
  geom_boxplot(fill= c('slateblue3', 'wheat3'), 
               outlier.color = NA) + 
  geom_jitter (alpha=0.2)  +
  xlab("Groups") + 
  ylab("Stemness") + 
  theme(axis.title.x = element_text(colour="black", size = 12, vjust=-1), 
        axis.title.y = element_text(colour="black", size = 12), 
        axis.text.x = element_text(angle = 45, vjust= 1, size = 12, hjust = 1)) +
  #facet_wrap(~ sample_char) +
  ggtitle("Stemness diff. beetween groups (p = 0.01745)") 

print(p)
dev.off()











