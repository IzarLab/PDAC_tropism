# This script plots PDAC landscape of CUIMC cohort

library(Seurat)
library(SeuratObject)
options(Seurat.object.assay.version = "v3")
library(SeuratWrappers)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)
library(xlsx)
library(monocle3)
library(parallel)

load(file = 'compartments/healthy_cancer.RData')
DefaultAssay(s_objs) = 'RNA'

p_ = DimPlot(s_objs, group.by = c('subcompartment','condition'), label = T, pt.size = .5, label.size = 4)
ggsave(plot = p_, filename = 'h-i.pdf', device = 'pdf', width = 21, height = 7)
