# predict endothelial cell subtype using an in vivo reference dataset

library(Seurat) # v4.3
library(ggplot2)
library(dplyr)
set.seed(1234)


# previously published scRNA-Seq of choroidal neovascularization
# Rohlenova et al 2020, https://doi.org/10.1016/j.cmet.2020.03.009
# raw counts and cell annotations downloaded from https://carmelietlab.sites.vib.be/en/software-tools

mouse_eye <- as(as.matrix(data.table:::fread("./mouse eye - Raw counts/Data.csv", sep=",", stringsAsFactors = F),rownames=1), "sparseMatrix")
rownames(mouse_eye) <- toupper(rownames(mouse_eye))
colnames(mouse_eye) <- paste0("mouse_eye_", colnames(mouse_eye))
mouse_eye.metadata <- read.table("./mouse eye - Raw counts/Metadata.csv", header = T, stringsAsFactors = F, sep=",")
mouse_eye.metadata$Observation <- paste0("mouse_eye_", mouse_eye.metadata$Observation)
mouse_eye.metadata$Cluster_withSample <- paste0(mouse_eye.metadata$Cluster, " - mouse_eye")
mouse_eye.metadata$Cluster_withSample_withExtraInfo <- paste0(mouse_eye.metadata$Cluster, " - ", mouse_eye.metadata$Condition, " - mouse_eye")
rownames(mouse_eye.metadata) <- mouse_eye.metadata$Observation
mouse_eye.metadata <- mouse_eye.metadata[mouse_eye.metadata$Endothelial.cell == "Yes", ]
mouse_eye <- mouse_eye[ , colnames(mouse_eye) %in% rownames(mouse_eye.metadata)]

CarmelietData.singleObj <- CreateSeuratObject(
  counts = mouse_eye, 
  meta.data = mouse_eye.metadata,
  names.delim = NA,
  project = "mouse_eye", 
  assay = "RNA",
  min.cells = 0, 
  min.features = 0
)


# bead assay endothelial cells Seurat object generated previously 
beadAssaySeurat.ECs
# An object of class Seurat 
# 36601 features across 6482 samples within 1 assay 
# Active assay: RNA (36601 features, 2000 variable features)
#  3 layers present: counts, data, scale.data
#  2 dimensional reductions calculated: pca, umap

# keep the common genes between the two objects
commonGenes <- intersect(rownames(CarmelietData.singleObj), rownames(beadAssaySeurat.ECs))
length(commonGenes)
# 13784

# recreate the Seurat object using the common subset of the genes
CarmelietData.singleObj <- CreateSeuratObject(counts = GetAssayData(CarmelietData.singleObj, assay = "RNA", slot="counts")[commonGenes, ], 
                                              meta.data = CarmelietData.singleObj@meta.data,
                                              names.delim = NA,
                                              project = "mouse_eye", 
                                              assay = "RNA", 
                                              min.cells = 0, 
                                              min.features = 0)

# downsample to keep max 1000 cells per ident
Idents(CarmelietData.singleObj) <- "Cluster_withSample"
set.seed(123456)
CarmelietData.singleObj <- subset(CarmelietData.singleObj, downsample=1000)

# Seurat RNA workflow
DefaultAssay(CarmelietData.singleObj) <- "RNA"
CarmelietData.singleObj <- NormalizeData(CarmelietData.singleObj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
CarmelietData.singleObj <- RunUMAP(CarmelietData.singleObj, dims = 1:30, reduction = "pca", return.model = TRUE)
CarmelietData.singleObj <- FindNeighbors(CarmelietData.singleObj, reduction = "pca", dims = 1:30) %>% FindClusters()

# plot original authors' cell annotations
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " eye neovascularization ECs - type annotations.pdf"), width = 8, height = 6)
DimPlot(CarmelietData.singleObj, label=T, repel = T, group.by =  "Cluster",  pt.size = 0.15) +
  ggtitle("Reference dataset EC populations")
dev.off()

# gene expression UMAP
genesToPlot <- c("CXCR4", "TOP2A")
for (myGene in genesToPlot) {
  p <- FeaturePlot(CarmelietData.singleObj, features = myGene, pt.size = 0.75, order=T,cols = c("lightgray", "red"))
  pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " umap in vivo ECs - ", myGene,".pdf"), width = 6.5, height = 5.5)
  print(p)
  dev.off()
}


# predict labels for in vitro cells using the in vivo reference:
beadAssay.predictLabels <- beadAssaySeurat.ECs
predictionAnchors <- FindTransferAnchors(reference = CarmelietData.singleObj,
                                         query = beadAssay.predictLabels,
                                         dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = predictionAnchors, 
                            refdata = CarmelietData.singleObj$Cluster, # "Cluster" has cell types by the original authors
                            dims = 1:30)
beadAssay.predictLabels <- AddMetaData(beadAssay.predictLabels, metadata = predictions)

# map query
beadAssay.predictLabels <- MapQuery(anchorset = predictionAnchors, 
                                    reference = CarmelietData.singleObj, 
                                    query = beadAssay.predictLabels,
                                    refdata = list(celltype = "Cluster"), # "Cluster" has cell types by the original authors
                                    reference.reduction = "pca", reduction.model = "umap")

# how many received a majority (reliable) prediction?
summary(beadAssay.predictLabels$predicted.celltype.score > 0.5)
# Mode   FALSE    TRUE 
# logical     930    5552 
beadAssay.predictLabels$prediction.reliable <- beadAssay.predictLabels$predicted.celltype.score > 0.5

# plot: prediction reliable (majority score) or not
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " ECs - prediction reliable.pdf"), width = 5.5, height = 5)
DimPlot(beadAssay.predictLabels, label=F, group.by = "prediction.reliable",  pt.size = 0.75) + 
  ggtitle("Prediction reliable (majority score)")
dev.off()

# plot prediction scores
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " ECs - prediction score.pdf"), width = 5, height = 4.8)
FeaturePlot(beadAssay.predictLabels, features = "predicted.celltype.score",  pt.size = 0.75) + 
  scale_y_continuous(breaks=seq(-5,5, by = 2.5)) +
  scale_x_continuous(breaks=seq(-3,6, by = 3)) +
  ggtitle("Prediction score")
dev.off()

# remove cells with unreliable prediction scores
beadAssay.predictLabels <- beadAssay.predictLabels[ , beadAssay.predictLabels$prediction.reliable]

# plot: predicted cell types
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " ECs - predicted type.pdf"), width = 5.5, height = 4.5)
DimPlot(beadAssay.predictLabels, label=T, group.by =  "predicted.celltype",  pt.size = 0.75, label.size = 6) + 
  ggtitle("Predicted cell type")
dev.off()

# plot query cells in reference UMAP coordinates

# to match axis limits
temp.for_UMAP_limits <- rbind(Embeddings(CarmelietData.singleObj, reduction = "umap"),
                              Embeddings(beadAssay.predictLabels, reduction = "ref.umap"))
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " ECs - in reference UMAP.pdf"), width = 7.5, height = 6)
DimPlot(beadAssay.predictLabels, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, repel = TRUE,  pt.size = 0.15) + 
  ggtitle("Query cells") + 
  xlim(min(temp.for_UMAP_limits[,1]), max(temp.for_UMAP_limits[,1])) + 
  ylim(min(temp.for_UMAP_limits[,2]), max(temp.for_UMAP_limits[,2])) +
  theme(panel.grid.major=element_line(colour="gray", linetype = "dotted"),  panel.background = element_rect(colour = "black")) + 
  xlab("UMAP_1") + ylab("UMAP_2")
dev.off()






