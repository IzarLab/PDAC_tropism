# This script plots UMAP of clusters and histology of P10 (PN14 in Supplementary Table 1) patient

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)

load(file = 'Misc/PN14.RData')
DefaultAssay(s_obj) = 'RNA'

s_obj$histology = factor('healthy', levels = c('cancer','healthy'))
s_obj$histology[s_obj$seurat_clusters %in% c(2,4,5,8,9,13)] = as.factor('cancer')

p_ = DimPlot(s_obj, group.by = c('seurat_clusters','histology'), label = T, pt.size = .5, label.size = 4) & NoLegend()
ggsave(plot = p_, filename = 'd_f.pdf', device = 'pdf', width = 14, height = 7)
