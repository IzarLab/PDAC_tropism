# This script plot UMAP projection of cancer compartment

library(Seurat)
library(SeuratObject)
options(Seurat.object.assay.version = "v4")
library(ggplot2)

load(file = 'compartments/cancer.RData')
DefaultAssay(s_objs) = 'RNA'

DimPlot(s_objs, group.by = 'cell_type', label = T, label.size = 5, pt.size = 1, repel = T) & NoLegend()
ggsave(filename = 'cancer_compartment.pdf', device = 'pdf', width = 7, height = 7)
