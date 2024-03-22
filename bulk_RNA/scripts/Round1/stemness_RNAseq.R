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



# Stemness index prediction script ---------
load("/media/ResearchHome/plummgrp/home/common/LFurtado-colab/scripts_git/bulk_RNA/model_RNA_MALTA.2018.Rda") # mm = modelo
w = mm$w
w[1:5]
length(w) #12953 (number of genes on the model)
# Filter matrix expression by the genes on stemenss model
matrix = count_matrix
# matrix = normalized_counts # same thing than raw count for this stemness prediction 
length(intersect(rownames(matrix), names(w))) #12773 
predict.DATA = matrix[rownames(matrix) %in% names(w) ,]
length(rownames(predict.DATA)) # 12773
w = w[ rownames(predict.DATA) ]
length(intersect(names(w),rownames(predict.DATA))) # 12773 
w[1:5]
length(names(w)) # 12773
is.vector(w) #TRUE

# Score the Matrix `X` using Spearman correlation
s = apply( predict.DATA, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
s[1:5]

# Scale the scores to be between 0 and 1
s = s - min(s)
s = s / max(s)
s[1:5]
s = as.data.frame(s)
names(s) = "stemness"
s$smpID <- rownames(s)


labsheet <- merge(labsheet, s, by='smpID')


# Boxplot 
library(ggplot2); theme_set(theme_classic()) 
ggplot(labsheet, aes(x=group,y=stemness, fill = group)) + 
  scale_fill_manual(values=c('slateblue3', 'wheat3'), name="Group") +
  geom_boxplot(outlier.color = NA)+ geom_jitter (alpha=0.5)  +
  xlab("") + 
  ylab("") + 
  theme(axis.title.x = element_text(colour="black", size = 12, vjust=-1), 
        axis.title.y = element_text(colour="black", size = 12), 
        axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1)) + 
  ggtitle("Stemness (RNAseq)") 


t.test(labsheet[labsheet$group %in% 'ETV6-NTRK3_fused_tumor', ]$stemness, 
       labsheet[labsheet$group %in% 'kinase-fused_tumor', ]$stemness)




