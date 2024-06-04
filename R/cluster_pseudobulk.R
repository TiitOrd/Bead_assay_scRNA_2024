# cluster pseudobulk expression calculation and heatmaps

library(pheatmap)
library(Seurat) # v4.3
library(ggplot2)
library(dplyr)
set.seed(1234)


# bead assay Seurat object generated previously 
beadAssaySeurat
# An object of class Seurat 
# 36601 features across 6946 samples within 1 assay 
# Active assay: RNA (36601 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap


# cells are annotated as either fibroblast or one of 5 EC clusters (EC 1...5)
Idents(beadAssaySeurat) <- "FB_and_EC_clustering_5way"


# cluster average expression for all genes
cluster.averages <- AverageExpression(beadAssaySeurat,
                                      assays="RNA",
                                      return.seurat = F,
                                      use.scale = F,
                                      use.counts = F)
cluster.averages <- cluster.averages[["RNA"]]
colSums(cluster.averages)
# all are 10000
# rescale values to be per million UMI-s (transcripts)
cluster.averages <- cluster.averages*100
colSums(cluster.averages)
# all now sum to 1e6
cluster.averages.df <- data.frame(cluster.averages)
write.table(cluster.averages.df, 
            file=paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " cluster average expression TPM.txt"),
            col.names = T, row.names = T, sep="\t", dec=",")






# cluster pseudobulk expression heatmap

# example gene list: tube formation genes
genesToPlot <- c("KIF20B","PFN1","VASP","STIL","TMED2","MTHFD1","MTHFD1L","PRICKLE1","BMP4","BCL10","PRKACB","RALA")

# subset pseudobulk expression table
genesAvgExpr <- data.frame(cluster.averages[genesToPlot , ])

# keep only endothelial clusters
genesAvgExpr <- genesAvgExpr[ , c(1:5)]

# scale each gene by the max level of that gene (thus, gene expression will be from 0 to max observed and in linear scale)
genesAvgExpr <- t(apply(genesAvgExpr, MARGIN=1, function(x) { x / max(x) }))

myPheatmap <- pheatmap::pheatmap(
  mat=t(genesAvgExpr),
  scale = "none",
  color = colorRampPalette(c("white", RColorBrewer::brewer.pal(n = 9, name = "Reds")[1:5], "darkred"))(50),
  cluster_cols = F,
  cluster_rows = F,
  angle_col = 45,
  main="Gene expression (0...Max)",
  cellwidth=15,
  cellheight = 15,
  border_color = "lightgray"
)
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " selected genes pseudobulk heatmap.pdf"), width = 9, height = 3)
print(myPheatmap)
dev.off()





