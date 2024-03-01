### --------------- Pathway Analysis ------------------- ###
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)
FG05_high_176_genes <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/FG05_high_176_genes.rds")
FG05_low_268_genes <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/FG05_low_268_genes.rds")



up_genes = rownames(FG05_high_176_genes)

ego2 <- enrichGO(gene = up_genes,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 # ont	One of "BP", "MF", and "CC" subontologies,                                   # or "ALL" for all three
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego2@result[order(ego2@result$Count, decreasing = TRUE),][1:10,]
barplot(ego2, drop=TRUE, showCategory=12, main = "up genes (EVT6_orientation)")
clusterProfiler::dotplot(ego2) + ggtitle("up genes (EVT6_orientation)")










down_genes = rownames(FG05_low_268_genes)

ego2 <- enrichGO(gene = down_genes,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 # ont	One of "BP", "MF", and "CC" subontologies,                                   # or "ALL" for all three
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego2@result[order(ego2@result$Count, decreasing = TRUE),][1:10,]
barplot(ego2, drop=TRUE, showCategory=12, main = "down genes (EVT6_orientation)")
clusterProfiler::dotplot(ego2) + ggtitle("down genes (EVT6_orientation)")

