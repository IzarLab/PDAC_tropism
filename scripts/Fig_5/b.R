# This script shows how clustering/annotation is done in nonmalignant exocrine compartment.

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(SeuratWrappers)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(monocle3)

# cell types of integrated data ####

load('compartments/exocrine.RData')

# Clustering ####

DefaultAssay(s_objs) = 'integrated'
cds = as.cell_data_set(s_objs)
cds = cluster_cells(cds, resolution = .0006)
s_objs$seurat_clusters = Idents(s_objs) = clusters(cds)

s_objs$acinar = colMeans(s_objs[["RNA"]]@data[c('CPB1','PRSS1'),])
s_objs$ductal = colMeans(s_objs[["RNA"]]@data[c('CFTR'),, drop = F])
s_objs$low_PanIN = colMeans(s_objs[["RNA"]]@data[c('DCLK1'),,drop = F])
s_objs$lesion = colMeans(s_objs[["RNA"]]@data[c('MUC1', 'MUC3A', 'MUC5AC', 'MUC5B', 'MUC16','TFF1','TFF2','TFF3','KLF5','KRT7','KRT19'),])
FeaturePlot(s_objs, features = c('acinar','ductal','low_PanIN','lesion'),
                 cols = c('grey80','red3'),
                 pt.size = .4, raster = F, order = T,
                 label = T, label.size = 6, repel = F, ncol = 2, max.cutoff = 'q99')

ggsave(filename = 'b.pdf', width = 15, height = 15, units = 'in', device = 'pdf')
