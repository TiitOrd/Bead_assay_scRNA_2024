# gene ontology enrichment, term similiarity filtering and plotting

library(GOSemSim)
library(rrvgo)
library(pheatmap)
library(Seurat) # v4.3
library(ggplot2)
library(dplyr)
set.seed(1234)


# bead assay endothelial cells Seurat object generated previously 
beadAssaySeurat.ECs
# An object of class Seurat 
# 36601 features across 6482 samples within 1 assay 
# Active assay: RNA (36601 features, 2000 variable features)
#  3 layers present: counts, data, scale.data
#  2 dimensional reductions calculated: pca, umap

# cell annotations are the 5-way EC clustering
Idents(beadAssaySeurat.ECs) <- "RNA_snn_res.0.4"


# generate the set of background genes to be used in enrichment testing
# the background will be all genes expressed at >1 TPM in any EC cluster at the pseudobulk level

# cluster pseudobulk expression levels for all genes
cluster.averages <- AverageExpression(beadAssaySeurat.ECs,
                                      assays="RNA",
                                      return.seurat = F,
                                      use.scale = F,
                                      use.counts = F)
cluster.averages <- cluster.averages[["RNA"]]
colSums(cluster.averages)
# all are 10000
# scale to transcripts per million
cluster.averages <- cluster.averages*100
colSums(cluster.averages)
# all are 1e6
cluster.averages.df <- data.frame(cluster.averages)

# for each gene, max pseudobulk expression observed in any EC cluster:
cluster.averages.max <- apply(cluster.averages.df, 1, max)

# how many genes were >1 TPM in at least one cluster?
summary(cluster.averages.max > 1)
# Mode   FALSE    TRUE 
# logical   22203   14398

# keep TPM > 1 genes for use as the enrichment testing background
beadAssaySeurat.ECs.TPM1_genes <- sort(names(cluster.averages.max)[cluster.averages.max > 1])
write.table(beadAssaySeurat.ECs.TPM1_genes,
            paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " ECs TPM1 genes (in at least one cluster).txt"),
            sep="\t", dec=",", col.names = F, row.names = F, quote = F)





# enrichment analysis: was run on the gProfiler website using up to top 500 markers per cluster. Custom background (genes selected above).

# load a gProfiler multiquery result CSV downloaded from the website (all results = not filtered for pval)
# this includes GP:BP, MF, CC all in one table
myGO <- read.csv("./gProfiler.csv")



# term similarity filtering with rrvgo
# 
# based on https://bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html
# 
# GP:BP, MF, CC will be processed separately in rrvgo since term similarity matrices will be calculated within each GO parent.


###
# for GO:BP:

# pre-calculate similarity data (will be used multiple times)
semdata <- GOSemSim::godata("org.Hs.eg.db", ont="BP")

# collect categories to plot:
# for each cell cluster, select significant terms, remove very large terms, reduce terms with rrvgo, keep up to top 10 terms
pvalCols <- grep("adjusted_p_value__cluster_", colnames(myGO), value = T)
reducedTerms.collect <- data.frame()
for (currentCol in pvalCols) { 
  
  # select signif GO:BP terms for one cluster, remove large terms
  myGO.currentSignif <- myGO
  myGO.currentSignif$current_padj <- myGO.currentSignif[ , currentCol]
  myGO.currentSignif <- myGO.currentSignif %>% filter(source=="GO:BP" & current_padj < 0.05) %>% arrange(current_padj) %>% filter(term_size < 500)
  
  # pairwise similarities of the selected GO terms
  simMatrix <- calculateSimMatrix(myGO.currentSignif$term_id,
                                  semdata=semdata,
                                  orgdb="org.Hs.eg.db",
                                  ont="BP",
                                  method="Rel")
  
  # scores are interpreted in the direction that higher are better
  scores <- setNames(-log10(myGO.currentSignif$current_padj), myGO.currentSignif$term_id)
  
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.6,
                                  orgdb="org.Hs.eg.db")
  
  # keep the most significant term from each similarity group
  reducedTerms.rmRedundant <- reducedTerms[!duplicated(reducedTerms$parentTerm), ]
  
  reducedTerms.rmRedundant$pval_column <- currentCol
  reducedTerms.collect <- rbind(reducedTerms.collect, reducedTerms.rmRedundant)
}

# limit categories contributed by each cell cluster
reducedTerms.collect.limit <- reducedTerms.collect[reducedTerms.collect$cluster <= 10, ]

dataToPlot <- myGO[match(unique(reducedTerms.collect.limit$go), myGO$term_id), ]

rownames(dataToPlot) <- dataToPlot$term_name
dataToPlot <- dataToPlot[ , pvalCols]
colnames(dataToPlot)  <- sub("adjusted_p_value__cluster_", "", colnames(dataToPlot))
dataToPlot <- -log10(dataToPlot)

# use the 95% quantile as the most intense color cutoff:
temp.strongestColorCutoff <- quantile(as.matrix(dataToPlot), probs = c(0.95))
colorBreaks=c(seq(0, temp.strongestColorCutoff, length=100))
hmColors <- gplots::colorpanel(n=length(colorBreaks)-1, 
                               low="gray95", 
                               mid = "orange",
                               high="mediumvioletred")

pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " GOBP heatmap top500.pdf"), width = 6.6, height = 8)
pheatmap::pheatmap(
  dataToPlot,
  scale="none",
  breaks=colorBreaks,
  color=hmColors,
  cluster_rows = F,
  cluster_cols = F,
  angle_col=0,
  main=paste0("GO:BP enrichment in EC cluster markers\n(top 500 genes)\nscale: -log10(padj)")
)
dev.off()
###


###
# for GO:MF:

# term similarity filtering based on https://bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html

# pre-calculate similarity data (will be used multiple times)
semdata <- GOSemSim::godata("org.Hs.eg.db", ont="MF")

# collect categories to plot:
# for each cell cluster, select significant terms, remove very large terms, reduce terms with rrvgo, keep up to top 10 terms
pvalCols <- grep("adjusted_p_value__cluster_", colnames(myGO), value = T)
reducedTerms.collect <- data.frame()
for (currentCol in pvalCols) { 
  
  # select signif GO:MF terms for one cluster, remove large terms
  myGO.currentSignif <- myGO
  myGO.currentSignif$current_padj <- myGO.currentSignif[ , currentCol]
  myGO.currentSignif <- myGO.currentSignif %>% filter(source=="GO:MF" & current_padj < 0.05) %>% arrange(current_padj) %>% filter(term_size < 500)
  
  # pairwise similarities of the selected GO terms
  simMatrix <- calculateSimMatrix(myGO.currentSignif$term_id,
                                  semdata=semdata,
                                  orgdb="org.Hs.eg.db",
                                  ont="MF",
                                  method="Rel")
  
  # scores are interpreted in the direction that higher are better
  scores <- setNames(-log10(myGO.currentSignif$current_padj), myGO.currentSignif$term_id)
  
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.6,
                                  orgdb="org.Hs.eg.db")
  
  # keep the most significant term from each similarity group
  reducedTerms.rmRedundant <- reducedTerms[!duplicated(reducedTerms$parentTerm), ]
  
  reducedTerms.rmRedundant$pval_column <- currentCol
  reducedTerms.collect <- rbind(reducedTerms.collect, reducedTerms.rmRedundant)
}

# limit categories contributed by each cell cluster
reducedTerms.collect.limit <- reducedTerms.collect[reducedTerms.collect$cluster <= 10, ]

dataToPlot <- myGO[match(unique(reducedTerms.collect.limit$go), myGO$term_id), ]

rownames(dataToPlot) <- dataToPlot$term_name
dataToPlot <- dataToPlot[ , pvalCols]
colnames(dataToPlot)  <- sub("adjusted_p_value__cluster_", "", colnames(dataToPlot))
dataToPlot <- -log10(dataToPlot)

# use the 95% quantile as the most intense color cutoff:
temp.strongestColorCutoff <- quantile(as.matrix(dataToPlot), probs = c(0.95))
colorBreaks=c(seq(0, temp.strongestColorCutoff, length=100))
hmColors <- gplots::colorpanel(n=length(colorBreaks)-1, 
                               low="gray95", 
                               mid = "orange",
                               high="mediumvioletred")

pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " GOMF heatmap top500.pdf"), width = 6, height = 6)
pheatmap::pheatmap(
  dataToPlot,
  scale="none",
  breaks=colorBreaks,
  color=hmColors,
  cluster_rows = F,
  cluster_cols = F,
  angle_col=0,
  main=paste0("GO:MF enrichment in EC cluster markers\n(top 500 genes)\nscale: -log10(padj)")
)
dev.off()
###



###
# for GO:CC:

# term similarity filtering based on https://bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html

# pre-calculate similarity data (will be used multiple times)
semdata <- GOSemSim::godata("org.Hs.eg.db", ont="CC")

# collect categories to plot:
# for each cell cluster, select significant terms, remove very large terms, reduce terms with rrvgo, keep up to top 10 terms
pvalCols <- grep("adjusted_p_value__cluster_", colnames(myGO), value = T)
reducedTerms.collect <- data.frame()
for (currentCol in pvalCols) { 
  
  # select signif GO:CC terms for one cluster, remove large terms
  myGO.currentSignif <- myGO
  myGO.currentSignif$current_padj <- myGO.currentSignif[ , currentCol]
  myGO.currentSignif <- myGO.currentSignif %>% filter(source=="GO:CC" & current_padj < 0.05) %>% arrange(current_padj) %>% filter(term_size < 500)
  
  # pairwise similarities of the selected GO terms
  simMatrix <- calculateSimMatrix(myGO.currentSignif$term_id,
                                  semdata=semdata,
                                  orgdb="org.Hs.eg.db",
                                  ont="CC",
                                  method="Rel")
  
  # scores are interpreted in the direction that higher are better
  scores <- setNames(-log10(myGO.currentSignif$current_padj), myGO.currentSignif$term_id)
  
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.6,
                                  orgdb="org.Hs.eg.db")
  
  # keep the most significant term from each similarity group
  reducedTerms.rmRedundant <- reducedTerms[!duplicated(reducedTerms$parentTerm), ]
  
  reducedTerms.rmRedundant$pval_column <- currentCol
  reducedTerms.collect <- rbind(reducedTerms.collect, reducedTerms.rmRedundant)
}

# limit categories contributed by each cell cluster
reducedTerms.collect.limit <- reducedTerms.collect[reducedTerms.collect$cluster <= 10, ]

dataToPlot <- myGO[match(unique(reducedTerms.collect.limit$go), myGO$term_id), ]

rownames(dataToPlot) <- dataToPlot$term_name 
dataToPlot <- dataToPlot[ , pvalCols]
colnames(dataToPlot)  <- sub("adjusted_p_value__cluster_", "", colnames(dataToPlot))
dataToPlot <- -log10(dataToPlot)

# use the 95% quantile as the most intense color cutoff:
temp.strongestColorCutoff <- quantile(as.matrix(dataToPlot), probs = c(0.95))
colorBreaks=c(seq(0, temp.strongestColorCutoff, length=100))
hmColors <- gplots::colorpanel(n=length(colorBreaks)-1, 
                               low="gray95", 
                               mid = "orange",
                               high="mediumvioletred")

pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " GOCC heatmap top500.pdf"), width = 5.3, height = 7)
pheatmap::pheatmap(
  dataToPlot,
  scale="none",
  breaks=colorBreaks,
  color=hmColors,
  cluster_rows = F,
  cluster_cols = F,
  angle_col=0,
  main=paste0("GO:CC enrichment in EC cluster markers\n(top 500 genes)\nscale: -log10(padj)")
)
dev.off()
###

