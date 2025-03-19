# This scripts plots Fig. 8a; it carries out cell trajectory analysis on CAF compartment

library(ggplot2)
library(viridis)
library(monocle3)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratWrappers)
library(patchwork)
set.seed(42)

# Reading in data ####

load('compartments/CAF.RData')

# Converting Seurat to Monocle ####

DefaultAssay(s_objs) = 'RNA'                                              # automatically adds RNA assay's counts slot
cds = as.cell_data_set(x = s_objs)                                        # automatically adds UMAP and PCA from integrated2 assay as that's the only assay with dim reductions
rowData(cds)$gene_short_name = rowData(cds)$gene_name = rownames(cds)     # NEEDED for gene-based plotting
cds = estimate_size_factors(cds)                                          # NEEDED for plot_genes_in_pseudotime

# Trajectory works better with lower granularity like cell types instead of clusters ####

cds@clusters$UMAP$cluster_result$optim_res$membership = as.character(colData(cds)$cell_type)
cds@clusters$UMAP$clusters = colData(cds)$cell_type
cds@clusters$UMAP$partitions = factor(x = rep(1, ncol(cds)), levels = '1')
names(cds@clusters$UMAP$cluster_result$optim_res$membership) = names(cds@clusters$UMAP$clusters) = names(cds@clusters$UMAP$partitions) = colnames(cds)

# Learning trajectory #####

cds = learn_graph(cds, use_partition = TRUE, verbose = FALSE, close_loop = F,
                  learn_graph_control = list(minimal_branch_len = 11, nn.k = 6))

# trajectory plot

pdf(file = 'a.pdf', width = 8, height = 8)
p_ = plot_cells(cds,
                color_cells_by = "cell_type",
                label_groups_by_cluster = F,
                label_leaves= F,
                label_branch_points = T,
                group_label_size = 3,
                graph_label_size = 3,
                cell_size = .5)+
      theme(text = element_blank(), axis.ticks = element_blank())

s_objs[["integrated"]]@misc[['trajectory']] =  data.frame(x1 = p_$layers[[3]]$data$source_prin_graph_dim_1,  x2 = p_$layers[[3]]$data$target_prin_graph_dim_1,
                                                          y1 =  p_$layers[[3]]$data$source_prin_graph_dim_2, y2 =  p_$layers[[3]]$data$target_prin_graph_dim_2)

breaks_ = c(min(s_objs$nFeature_RNA), mean(range(s_objs$nFeature_RNA)), max(s_objs$nFeature_RNA))
FeaturePlot(s_objs, pt.size = .6, features = 'nFeature_RNA', order = T, coord.fixed = T, label = F, label.size = 2.5, raster = F)+
theme(legend.spacing.y = unit(0.5, 'cm'), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), text = element_text(size = 15, face = 'bold'))+
labs(title = NULL, color = 'Gene count')+
viridis::scale_color_viridis(option = 'H', breaks = breaks_, labels = round(breaks_))+
geom_segment(data =  p_$layers[[3]]$data, aes(x = source_prin_graph_dim_1, xend = target_prin_graph_dim_1,
                                              y = source_prin_graph_dim_2, yend = target_prin_graph_dim_2), color = 'red',linewidth = 1)

DimPlot(s_objs, group.by = 'cell_type', label = T, label.box = T, label.size = 4.5, repel = T, pt.size = 1, raster = F, order = c('ADM','ductal','ductal EMT','acinar','precursor lesion','differentiating ductal'))+
theme(legend.position = 'none', axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
      text = element_text(size = 10), plot.background = element_blank(), panel.background = element_blank())+
geom_segment(data =  p_$layers[[3]]$data, aes(x = source_prin_graph_dim_1, xend = target_prin_graph_dim_1,
                                              y = source_prin_graph_dim_2, yend =  target_prin_graph_dim_2), linewidth = 1)

graphics.off()
