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
                                                   gene[c(dim(etv6_mt)[2]+1) : dim(kinase_mt)[2]], exact=T) # excat = T nos da um resultado mais preciso, mas demora mais pra rodar o looping. Entao, exact = F é melhor custo-beneficio. . .
                                z <- zz$p.value
                                return(z)
                              }, mc.cores= 20))
wpvalues_adj <- p.adjust(wpvalues, method = "BH")
wpvalues_adj <- data.frame(wpvalues_adj)
wpvalues_adj$gene <- rownames(wpvalues_adj)

# Compute FoldChange

# Calculate the mean expression for each gene across conditions
mean_expr_etv6 <- rowMeans(DE_input_matrix[, etv6]) #row = genes; columns = samples
mean_expr_kinase <- rowMeans(DE_input_matrix[, kinase])
# Calculate Fold Change for each gene
fold_change <- mean_expr_etv6 / mean_expr_kinase # interpret it as "etv6 group has n times more expression upon a given gene compared to kinase group
fold_change <- data.frame(fold_change) 
fold_change$gene <- rownames(fold_change)

FG_pvalue <- merge(fold_change, wpvalues_adj, by='gene')
volcano <- data.frame(DE_input_matrix)
volcano$gene <- rownames(volcano)
volcano <- merge(volcano, FG_pvalue, by='gene')
head(volcano)

# Create labels for the volcano plot 
volcano$threshold <- "1" # threshold info will be used to define DMP and to color volcano plots
b <- volcano[volcano$wpvalues_adj < 0.05 & volcano$fold_change < -0.05,] #hypo NSC probes
volcano[rownames(b),"threshold"] <- "2"
c <- volcano[volcano$wpvalues_adj < 0.05 & volcano$fold_change > 0.05,] #hyper NSC probes 
volcano[rownames(c),"threshold"] <- "3"
print(table(volcano$threshold))

table(volcano$wpvalues_adj < 0.05) # none statistical significance =o make it again ...



