# NicheNet for bead assay cells

library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(circlize)
library(nichenetr) # v2.0.4
library(Seurat)
set.seed(1234)



# Seurat object generated previously; contains endothelial cells and fibroblasts
beadAssaySeurat.FBs_ECs
# An object of class Seurat 
# 36601 features across 6946 samples within 1 assay 
# Active assay: RNA (36601 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap

# cells are annotated as either fibroblast or one of 5 EC clusters (EC 1...5)
Idents(beadAssaySeurat) <- "FB_and_EC_clustering_5way"



# NicheNet analysis based on the tutorials:
# https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_wrapper.md
# https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md

# using the 'NicheNet-v2' networks published as ready-to-use RDS files by the NicheNet authors
lr_network = readRDS("lr_network_human_21122021.rds")
ligand_target_matrix = readRDS("ligand_target_matrix_nsga2r_final.rds")
weighted_networks = readRDS("weighted_networks_nsga2r_final.rds")

lr_network = lr_network %>% distinct(from, to)
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))


# seurat object for nichenet
seuratObj <- beadAssaySeurat.FBs_ECs

# nichenet setup:
# * aim to explain the signature of each EC subpopulation (its markers vs all other EC-s)
# * the ligand sender cells can be any of the five EC populations and also the fibroblasts
celltypes_to_test <- c("EC_1", "EC_2", "EC_3", "EC_4", "EC_5")
nichenet_res.list <- list()
for (currentReceiver in celltypes_to_test) {
  
  print(currentReceiver)
  
  nichenet_output = nichenet_seuratobj_cluster_de(
    seurat_obj = seuratObj, 
    receiver_affected = currentReceiver, # the "affected" cells
    receiver_reference = setdiff(celltypes_to_test, currentReceiver), # the other unaffected cells (= other ECs)
    sender = "all", # can be EC_1...5 or FB
    expression_pct=0.1, # Default: 0.10
    lfc_cutoff = 0.25, # Default: 0.25.
    geneset = "up", # only genes upregulated in "affected" cells
    filter_top_ligands = T, # Default: TRUE.
    top_n_ligands	= 5, # Default: 30.
    top_n_targets = 200, # Default = 200.
    cutoff_visualization = 0.33, # Default = 0.33.
    ligand_target_matrix = ligand_target_matrix, 
    lr_network = lr_network, 
    weighted_networks = weighted_networks
  )
  
  nichenet_res.list[[currentReceiver]] <- nichenet_output
}


# heatmap of ligand expressing populations
nichenet.topLigandsCombined <- sort(unique(unlist(lapply(nichenet_res.list, function(x) {x$top_ligands}))))
nichenet.topLigandsCombined

# TPM table by cluster for all genes
cluster.averages <- AverageExpression(beadAssaySeurat.FBs_ECs,
                                      group.by = "FB_and_EC_clustering_5way",
                                      assays="RNA",
                                      return.seurat = F,
                                      use.scale = F,
                                      use.counts = F)
cluster.averages <- cluster.averages[["RNA"]]
colSums(cluster.averages)
# all are 10000
# calculate per million (=TPM)
cluster.averages <- cluster.averages*100
colSums(cluster.averages)
# all are 1e6

# ligand expression in the possible sender cell types
nichenet.GenesAvgExpr <- cluster.averages[nichenet.topLigandsCombined , ]

temp.TableForExprHeatmap <- nichenet.GenesAvgExpr
colnames(temp.TableForExprHeatmap) <- sub("-", "_", colnames(temp.TableForExprHeatmap))

# scale by max
temp.TableForExprHeatmap.byMax <- t(apply(temp.TableForExprHeatmap, MARGIN=1, function(x) { x / max(x) }))

myPheatmap <- pheatmap::pheatmap(
  mat=temp.TableForExprHeatmap.byMax,
  scale = "none",
  color = colorRampPalette(c("white", RColorBrewer::brewer.pal(n = 9, name = "Reds")[1:5], "darkred"))(50),
  cluster_cols = T,
  cluster_rows = T,
  border_color = "lightgray"
)
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " NicheNet topLigandsCombined - sender expr heatmap.pdf"), width = 3, height = 4.5)
print(myPheatmap)
dev.off()


# collect top ligands and their activities for each EC population

# collect stats
currentActivities.combined <- data.frame()
# collect minitable plots
currentActivities.toptableGrobList <- list()
for (currentReceiver in names(nichenet_res.list)) {
  
  currentActivities <- data.frame(nichenet_res.list[[currentReceiver]]$ligand_activities[1:5, ])
  rownames(currentActivities) <- currentActivities$test_ligand
  currentActivities$receiver <- currentReceiver
  currentActivities$test_ligand <- factor(currentActivities$test_ligand, levels=rev(currentActivities$test_ligand))
  
  # minitable: plot one set of ligands as a table (table Grob in the library gridExtra)
  currentActivities.plotTable <- currentActivities[ , c("rank", "test_ligand", "aupr")]
  rownames(currentActivities.plotTable) <- NULL
  colnames(currentActivities.plotTable) <- c("Rank", "Ligand", "AUPR")
  currentActivities.toptableGrobList[[currentReceiver]] <- gridExtra::arrangeGrob(gridExtra::tableGrob(currentActivities.plotTable, rows=NULL), 
                                                                                   top = currentReceiver)
  
  # collect stats into table
  currentActivities.combined <- rbind(currentActivities.combined, currentActivities)
}

# plot the minitables
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " NicheNet topLigandsCombined - tables.pdf"), width = 8, height = 5)
gridExtra::grid.arrange(grobs=currentActivities.toptableGrobList, top="Ligands predicted to influence EC populations", ncol=3)
dev.off()



# circos plot
# based on the nichenet circos tutorial: https://github.com/saeyslab/nichenetr/blob/master/vignettes/circos.md

# select a nichenet run
names(nichenet_res.list)
# "EC_1" "EC_2" "EC_3" "EC_4" "EC_5"

nichenet_run_name <- "EC_2"


# take the run output
nichenet_output <- nichenet_res.list[[nichenet_run_name]]

# save the auto-generated ligand-receptor heatmap
p <- nichenet_output$ligand_receptor_heatmap +
  ggtitle(paste0("Cell cluster ", nichenet_run_name), subtitle = "Expressed receptors for the prioritized ligands") +
  labs(fill='Prior interaction potential\n(ligand-receptor pair)')
p
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " ", nichenet_run_name, " lig-rec heatmap.pdf"), width =4, height = 3)
p
dev.off()

# extract from the nichenet object
ligand_activities <- nichenet_output$ligand_activities
ligand_target_matrix <- nichenet_output$ligand_target_matrix
ligand_target_df <- nichenet_output$ligand_target_df
ligand_receptor_df <- nichenet_output$ligand_receptor_df

best_upstream_ligands = ligand_activities %>% top_n(5, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
head(best_upstream_ligands)

# ligand grouping (will be used to color ligands)
# here, will just color each ligand separately:
ligand_type_indication_df = tibble(
  ligand_type = c(best_upstream_ligands),
  ligand = c(best_upstream_ligands))

active_ligand_target_links_df = ligand_target_df
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "cluster_marker") %>% inner_join(ligand_type_indication_df) 
# filter links
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.66)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

# ligands and targets with no remaining links will be removed
ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

# give each segment of ligands and targets a specific color and order

# ligand coloring
ligandsToPlot.count <- length(unique(circos_links$ligand))
grid_col_ligand = RColorBrewer::brewer.pal(ligandsToPlot.count, "Dark2")
names(grid_col_ligand) <- unique(circos_links$ligand)

# target coloring: (named vector)
# there was only one target category here:
grid_col_target = c("cluster_marker" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # 2 extra spaces: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

# order ligands and targets
target_order = circos_links$target %>% unique()
ligand_order = unique(circos_links$ligand) # has the 2 extra spaces: make a difference between a gene as ligand and a gene as target!
order = c(ligand_order,target_order)

# define the gaps between the different segments
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

ligandTypes <- unique(circos_links$ligand_type)

gaps = c(
  # since in the current plotting, each ligand is it's own ligand group, generating the gaps between ligands and ligand groups can be simplified:
  rep(width_different_cell, times=length(ligandTypes)-1 ),
  
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "cluster_marker") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)

# Render the circos plot
pdf(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), " ", nichenet_run_name, " ligand-target circos.pdf"), width =10, height = 8)

circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,
             transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,
             annotationTrack = "grid", 
             preAllocateTracks = list(list(track.height = 0.15),
                                      list(track.height = 0.01)))
# customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #

dev.off()
circos.clear()
