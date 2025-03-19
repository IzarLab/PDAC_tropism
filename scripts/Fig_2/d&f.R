# This script plots Fig. 2d and f (UMAP of clusters and histology of P10). P10 refers to patient PN14 in the Supplementary Table 1.

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)

# Reading in snRNA-seq data of P10 (PN14) ####

load(file = 'Misc/PN14.RData')
DefaultAssay(s_obj) = 'RNA'

# Annotating cancer clusters ####

s_obj$histology = factor('healthy', levels = c('cancer','healthy'))
s_obj$histology[s_obj$seurat_clusters %in% c(2,4,5,8,9,13)] = as.factor('cancer')

# Plotting UMAPs ####

p_ = DimPlot(s_obj, group.by = c('seurat_clusters','histology'), label = T, pt.size = .5, label.size = 4) & NoLegend()
ggsave(plot = p_, filename = 'd&f.pdf', device = 'pdf', width = 14, height = 7)
