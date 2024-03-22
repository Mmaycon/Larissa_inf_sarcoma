### --- Select Most Variable Genes on Kinase groups for Reference Ontology --- ###

# Load packages -------------
library(DESeq2)
library(sva)
library(dplyr)
library(stringr)

# Load data 
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


# Load Normalized/raw counts 
dds <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/bulk_RNA/objects/Round2/DESeq2_output.rds")
kinase_smp <- labsheet[labsheet$group %in% 'kinase-fused_tumor', ]$smpID
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts <- normalized_counts[, colnames(normalized_counts) %in% kinase_smp]
dim(normalized_counts) #18777    17
raw_counts <- counts(dds, normalized=FALSE)
raw_counts <- raw_counts[, colnames(raw_counts) %in% kinase_smp]
dim(raw_counts) #18777    17



# hc clustering can't hand bulk RNAseq
# using only var function also can't capture the variance in the data
# NOW: DE 1 vs ALL -> keep low FC -> these are the most variable features








# Calculate dissimilarities
dist_matrix_nor <- dist(normalized_counts)
# Perform hierarchical clustering
hc_nor <- hclust(dist_matrix_nor)
# Plot the dendrogram
plot(hc_nor)

clusters <- cutree(hc_nor, k = 3)  # Adjust 'k' as needed for the desired number of clusters
table(clusters)

# Get the variable features based on the clustering
variable_features <- normalized_counts[clusters == 1, ]  # Adjust cluster number as needed








# Calculate dissimilarities
dist_matrix <- dist(raw_counts)
# Perform hierarchical clustering
hc <- hclust(dist_matrix)
# Plot the dendrogram
plot(hc)








# Find Most Variable Feature
# Calculate variance 
variance_values <- apply(raw_counts, 1, var)

# Create a data frame with CpG sites and their variance values
variance_data <- data.frame(genes = rownames(raw_counts), Variance = variance_values)
hist(variance_data$Variance)

# Visualize all probes on the elbow plot
library(ggplot2); theme_set(theme_classic())
elbow_plot <- variance_data %>%
  arrange(desc(Variance)) %>%
  mutate(Rank = row_number()) %>%
  ggplot(aes(x = Rank, y = Variance)) +
  geom_line() +
  geom_point() +
  labs(title = "Elbow Plot for Variance Cutoff",
       x = "Number of Genes",
       y = "Variance") +
  theme_minimal(); elbow_plot # XXXX seems a good cutoff 

variance_data[variance_data$Variance >=  XXXX, ] %>% dim() # Cutoff choosen by eye

# Select the top variance probes (n probes)
library(dplyr)
top_n <- n
selected_features <- variance_data %>%
  arrange(desc(Variance)) %>%
  slice(1:top_n)