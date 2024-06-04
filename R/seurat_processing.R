# bead assay scRNA dataset processing using Seurat

library(Seurat) # v4.3
library(ggplot2)
library(dplyr)
set.seed(1234)


current10X_RNA_filtered <- Read10X(data.dir = "filtered_feature_bc_matrix", strip.suffix=T)
dim(current10X_RNA_filtered)
# 36601  9454

beadAssaySeurat <- CreateSeuratObject(counts = current10X_RNA_filtered, 
                                      project = "bead_assay", 
                                      assay = "RNA",
                                      min.cells = 0, 
                                      min.features = 0)

# cell QC filtering

beadAssaySeurat[["percent.mt"]] <- PercentageFeatureSet(beadAssaySeurat, pattern = "^MT-", assay = "RNA")

png(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " before QC filtering.png"), width=1000, height=1000)
print(VlnPlot(beadAssaySeurat, features=c("percent.mt", "nCount_RNA", "nFeature_RNA")))
dev.off()

beadAssaySeurat <- subset(beadAssaySeurat, 
                          subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & 
                            nCount_RNA < 12000 & nCount_RNA > 1000 &
                            percent.mt < 10)

png(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " after QC filtering.png"), width=1000, height=1000)
print(VlnPlot(beadAssaySeurat, features=c("percent.mt", "nCount_RNA", "nFeature_RNA")))
dev.off()

# standard Seurat v4 RNA workflow
DefaultAssay(beadAssaySeurat) <- "RNA"
beadAssaySeurat <- NormalizeData(beadAssaySeurat) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
beadAssaySeurat <- FindNeighbors(beadAssaySeurat)
beadAssaySeurat <- FindClusters(beadAssaySeurat)
beadAssaySeurat <- RunUMAP(beadAssaySeurat, dims = 1:30)

# check QC parameters on umap
png(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " cell QC params.png"), width=500, height=1500)
FeaturePlot(beadAssaySeurat, features = c("percent.mt")) +
  FeaturePlot(beadAssaySeurat, features = c("nCount_RNA")) +
  FeaturePlot(beadAssaySeurat, features = c("nFeature_RNA")) +
  DimPlot(beadAssaySeurat, label = T)
dev.off()

# marker genes to locate endothelial cells and fibroblasts
png(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " cell type marker genes.png"), width=1500, height=1500)
FeaturePlot(beadAssaySeurat, 
            features = c("PECAM1", "CD34", "VWF", "DCN", "COL1A2"),
            cols = c("lightgray", "red"))
dev.off()

beadAssaySeurat
# An object of class Seurat 
# 36601 features across 6946 samples within 1 assay 
# Active assay: RNA (36601 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap





# remove fibroblasts and recluster

# in the clustering above, fibroblasts (DCN-high cells) are cluster 8
beadAssaySeurat.ECs <- subset(beadAssaySeurat, idents = "8", invert=T)
beadAssaySeurat.ECs
# An object of class Seurat 
# 36601 features across 6482 samples within 1 assay 
# Active assay: RNA (36601 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap

DefaultAssay(beadAssaySeurat.ECs) <- "RNA"
beadAssaySeurat.ECs <- NormalizeData(beadAssaySeurat.ECs) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
beadAssaySeurat.ECs <- FindNeighbors(beadAssaySeurat.ECs)
beadAssaySeurat.ECs <- FindClusters(beadAssaySeurat.ECs)
beadAssaySeurat.ECs <- RunUMAP(beadAssaySeurat.ECs, dims = 1:30)
DimPlot(beadAssaySeurat.ECs, label=T)

# for 5 clusters:
beadAssaySeurat.ECs <- FindClusters(beadAssaySeurat.ECs, resolution = 0.4)
DimPlot(beadAssaySeurat.ECs, label=T)

# check QC parameters on umap
png(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " cell QC params - EC clusters.png"), width=500, height=1500)
FeaturePlot(beadAssaySeurat.ECs, features = c("percent.mt")) +
  FeaturePlot(beadAssaySeurat.ECs, features = c("nCount_RNA")) +
  FeaturePlot(beadAssaySeurat.ECs, features = c("nFeature_RNA")) +
  DimPlot(beadAssaySeurat.ECs, label = T)
dev.off()

# add cell cycle scoring 
s.genes <- cc.genes$s.genes # included with Seurat
g2m.genes <- cc.genes$g2m.genes # included with Seurat
beadAssaySeurat.ECs <- CellCycleScoring(beadAssaySeurat.ECs, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
table(beadAssaySeurat.ECs$Phase)
# G1  G2M    S 
# 4072 1253 1157 

# cell cycle plots
png(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " bead assay scRNA - cell cycle phase - EC clusters.png"), width=1500, height=1500)
DimPlot(beadAssaySeurat.ECs, label = F, group.by = "Phase") +
  FeaturePlot(beadAssaySeurat.ECs, features = c("S.Score")) +
  FeaturePlot(beadAssaySeurat.ECs, features = c("G2M.Score")) +
  DimPlot(beadAssaySeurat.ECs, group.by = c("RNA_snn_res.0.4"), label = T)
dev.off()


# find markers for EC clusters (each cluster vs all others)
Idents(beadAssaySeurat.ECs) <- "RNA_snn_res.0.4"
DefaultAssay(beadAssaySeurat.ECs) <- "RNA"
beadAssaySeurat.ECs.markers <- FindAllMarkers(object = beadAssaySeurat.ECs, 
                                              test.use = "wilcox", # wilcox is the default
                                              min.pct = 0.1, # 0.1 is the default
                                              logfc.threshold = 0.25, # 0.25 is the default
                                              only.pos = T) # FALSE is the default
beadAssaySeurat.ECs.markers <- beadAssaySeurat.ECs.markers %>% filter(p_val_adj < 0.05)

top10 <- beadAssaySeurat.ECs.markers %>% group_by(cluster) %>% top_n(n = 10, avg_log2FC)
top15 <- beadAssaySeurat.ECs.markers %>% group_by(cluster) %>% top_n(n = 15, avg_log2FC)

# scale all genes
beadAssaySeurat.ECs.copyforplot <- beadAssaySeurat.ECs
beadAssaySeurat.ECs.copyforplot <- ScaleData(beadAssaySeurat.ECs.copyforplot, features = rownames(beadAssaySeurat.ECs.copyforplot))

# cluster markers heatmap
png(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " beadAssaySeurat.ECs.markers - cluster_markers_heatmap_logFC_top10 - ECs only.png"), width = 1500, height = 1500)
print(DoHeatmap(object = beadAssaySeurat.ECs.copyforplot, features = top10$gene))
dev.off()

# cluster markers heatmap
png(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " beadAssaySeurat.ECs.markers - cluster_markers_heatmap_logFC_top15 - ECs only.png"), width = 1500, height = 1500)
print(DoHeatmap(object = beadAssaySeurat.ECs.copyforplot, features = top15$gene))
dev.off()

rm(beadAssaySeurat.ECs.copyforplot)

# save the full list of cluster markers
write.table(beadAssaySeurat.ECs.markers, file = paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S")," beadAssaySeurat.ECs.markers.txt"))



# gene expression UMAP
myGene <- "MKI67"
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " ECs FBs - ", myGene,".pdf"), width = 6.6, height = 4.9)
FeaturePlot(beadAssaySeurat.ECs, features = myGene,  pt.size = 0.75,  cols = c("lightgray", "red"))
dev.off()

# two genes - gene expression UMAP
myGene1 <- "MATN2"
myGene2 <- "MKI67"
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " ECs - two genes.pdf"), width = 17, height = 4.8)
FeaturePlot(beadAssaySeurat.ECs, features = c(myGene1, myGene2), blend=T, blend.threshold = 0, pt.size = 0.75, 
            cols=c("lightgray","red","blue"))
dev.off()

# Violin plot
genesToPlot <- c("MEOX2", "ANGPT2", "APLN")
for (myGene in genesToPlot) {
  p <- VlnPlot(beadAssaySeurat.ECs, features = myGene, group.by = "RNA_snn_res.0.4", alpha = 1) +
    ylab("Expression level") +
    xlab("Endothelial cell cluster") +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
  p
  pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " violin ECs - ", myGene,".pdf"), width = 4, height = 4.5)
  print(p)
  dev.off()
}

# fraction of MKI67-expressing cells by cluster
myCounts <- GetAssayData(beadAssaySeurat.ECs, assay = "RNA", slot="counts")
counts_MKI67 <- cbind(myCounts["MKI67", ], as.character(beadAssaySeurat.ECs$RNA_snn_res.0.4))
colnames(counts_MKI67) <- c("MKI67", "EC_cluster")
counts_MKI67 <- data.frame(counts_MKI67)
counts_MKI67 <- counts_MKI67 %>% group_by(EC_cluster) %>% summarize(count_expressed = sum(MKI67>0), n=n())
counts_MKI67$frac_expr <- counts_MKI67$count_expressed / counts_MKI67$n
counts_MKI67

# pie chart for cell count in each cluster
cellNumberByCluster <- data.frame(table(beadAssaySeurat.ECs$RNA_snn_res.0.4))
# Compute the position of labels
data <- cellNumberByCluster
data$Var1 <- as.character(data$Var1)
data <- data %>% arrange(desc(Var1)) %>% mutate(prop = Freq / sum(data$Freq) *100) %>% mutate(ypos = cumsum(prop)- 0.5*prop )
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " ECs - pie chart.pdf"), width = 5, height = 4.5)
ggplot(data, aes(x="", y=prop, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  ggtitle("EC cell proportion by cluster") +
  geom_text(aes(y = ypos, label = Var1), color = "black", size=8)
dev.off()





# plot a set of genes as module score

# EC subtype markers listed in Table S2 of Rohlenova et al 2020 (https://doi.org/10.1016/j.cmet.2020.03.009)
# provides 50 markers per annotation

# read EC marker gene lists for identities implicated by label transfer: immature, tip, proliferating
markers.Rohlenova2020.CEC <- read.table("./data/Rohlenova_et_al_2020_EC_subtype_markers.txt", stringsAsFactors = F, sep="\t", header = T)
markers.Rohlenova2020.CEC$gene.capitalized <- toupper(markers.Rohlenova2020.CEC$Gene)

# add module scores
beadAssay.moduleScores <- beadAssaySeurat.ECs
for (currentName in sort(unique(markers.Rohlenova2020.CEC$EC_subset))) {
  currentGenes <- markers.Rohlenova2020.CEC$gene.capitalized[markers.Rohlenova2020.CEC$EC_subset == currentName]
  beadAssay.moduleScores <- AddModuleScore(beadAssay.moduleScores, 
                                           features = list(currentGenes), 
                                           name=paste0("score_",currentName))
}

# plot the signatures of interest individually

# immature
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " ECs - signature expr - immature.pdf"), width = 5, height = 5.3)
FeaturePlot(beadAssay.moduleScores, features = "score_CNV-EC immature1",  pt.size = 0.75, max.cutoff = "q90", min.cutoff = "q10") + 
  ggtitle("Gene signature expression: immature", subtitle = "Rohlenova et al 2020\nEC subset markers\n(50 genes per subset)") +
  scale_y_continuous(breaks=seq(-5,5, by = 2.5)) +
  scale_x_continuous(breaks=seq(-3,6, by = 3))
dev.off()

# proliferating
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " ECs - signature expr - proliferating.pdf"), width = 5, height = 5.3)
FeaturePlot(beadAssay.moduleScores, features = "score_CNV-EC proliferating1",  pt.size = 0.75, max.cutoff = "q90", min.cutoff = "q10") + 
  ggtitle("Gene signature expression: proliferating", subtitle = "Rohlenova et al 2020\nEC subset markers\n(50 genes per subset)") +
  scale_y_continuous(breaks=seq(-5,5, by = 2.5)) +
  scale_x_continuous(breaks=seq(-3,6, by = 3))
dev.off()

# tip cell
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " ECs - signature expr - tip cell.pdf"), width = 5, height = 5.3)
FeaturePlot(beadAssay.moduleScores, features = "score_CNV-EC tip cell1",  pt.size = 0.75, max.cutoff = "q90", min.cutoff = "q10") + 
  ggtitle("Gene signature expression: tip cell", subtitle = "Rohlenova et al 2020\nEC subset markers\n(50 genes per subset)") +
  scale_y_continuous(breaks=seq(-5,5, by = 2.5)) +
  scale_x_continuous(breaks=seq(-3,6, by = 3))
dev.off()






# subset the most tip-like EC cluster from the 5-way EC clustering and run subclustering on that subset

# currently the cells of interest are named cluster 1 as the numbering is 0...4. 
# (later it is named EC_2 when 5-way EC clustering is presented as EC_1...5)
beadAssaySeurat.ECs.subset <- subset(beadAssaySeurat.ECs, idents = "1", invert=F)

DefaultAssay(beadAssaySeurat.ECs.subset) <- "RNA"
beadAssaySeurat.ECs.subset <- NormalizeData(beadAssaySeurat.ECs.subset) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
beadAssaySeurat.ECs.subset <- FindNeighbors(beadAssaySeurat.ECs.subset)
beadAssaySeurat.ECs.subset <- FindClusters(beadAssaySeurat.ECs.subset, resolution = 0.6)
beadAssaySeurat.ECs.subset <- RunUMAP(beadAssaySeurat.ECs.subset, dims = 1:30)
DimPlot(beadAssaySeurat.ECs.subset, label=T)

# run FindAllMarkers
Idents(beadAssaySeurat.ECs.subset) <- "RNA_snn_res.0.6"
DefaultAssay(beadAssaySeurat.ECs.subset) <- "RNA"
beadAssaySeurat.ECs.subset.markers <- FindAllMarkers(object = beadAssaySeurat.ECs.subset, 
                                                     test.use = "wilcox", # wilcox is the default
                                                     min.pct = 0.1, # 0.1 is the default
                                                     logfc.threshold = 0.25, # 0.25 is the default
                                                     only.pos = T) # FALSE is the default
beadAssaySeurat.ECs.subset.markers <- beadAssaySeurat.ECs.subset.markers %>% filter(p_val_adj < 0.05)

top10 <- beadAssaySeurat.ECs.subset.markers %>% group_by(cluster) %>% top_n(n = 10, avg_log2FC)
top15 <- beadAssaySeurat.ECs.subset.markers %>% group_by(cluster) %>% top_n(n = 15, avg_log2FC)

# scale all genes for heatmap
beadAssaySeurat.ECs.subset.copyforplot <- beadAssaySeurat.ECs.subset
beadAssaySeurat.ECs.subset.copyforplot <- ScaleData(beadAssaySeurat.ECs.subset.copyforplot, features = rownames(beadAssaySeurat.ECs.subset.copyforplot))

# cluster markers heatmap
png(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " beadAssaySeurat.ECs.subset.markers - cluster_markers_heatmap_logFC_top10 - ECs only.png"), width = 1500, height = 1500)
print(DoHeatmap(object = beadAssaySeurat.ECs.subset.copyforplot, features = top10$gene))
dev.off()

# cluster markers heatmap
png(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " beadAssaySeurat.ECs.subset.markers - cluster_markers_heatmap_logFC_top15 - ECs only.png"), width = 1500, height = 1500)
print(DoHeatmap(object = beadAssaySeurat.ECs.subset.copyforplot, features = top15$gene))
dev.off()

rm(beadAssaySeurat.ECs.subset.copyforplot)

# save the full list of cluster markers
write.table(beadAssaySeurat.ECs.subset.markers, file = paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S")," beadAssaySeurat.ECs.subset.markers - cluster markers full table.txt"))




