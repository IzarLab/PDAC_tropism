# This scripts plots Fig. 6a, b and d; generates UMAP and dot plots of immune compartment

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(SeuratWrappers)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)

# Loading data ####

load('compartments/T.RData'); T_ = s_objs
load('compartments/B.RData'); B_ = s_objs
load('compartments/myeloid.RData'); my_ = s_objs; rm(s_objs)
DefaultAssay(T_) = DefaultAssay(my_) = 'RNA'

# UMAP of compartment ####

pdf('a-b&d.pdf', width = 12, height = 12)

B_[["umap"]]@cell.embeddings = B_[["umap"]]@cell.embeddings+5
t_b = merge(x = T_, y = B_, merge.data = T, merge.dr = T)

DimPlot(t_b, pt.size = 2, group.by = 'cell_type', label = T, label.box = T, label.size = 4, repel = T) & NoLegend()

DimPlot(my_, pt.size = 2, group.by = 'cell_type', label = T, label.box = T, label.size = 4, repel = T) & NoLegend()

# Dot plot exhaustion receptors ####

s_ = T_
Idents(s_) = s_$cell_type
DotPlot(object = s_, features = c('HAVCR2','TIGIT','CTLA4','PDCD1'),      # exhaustion markers (inhibitory immune checkpoints)
        scale = F,
        dot.min = 0.05, dot.scale = 15,
        scale.by = 'size',
        cols = c('red','blue'),
        cluster.idents = F, )+
  theme(axis.text.x = element_text(size = 20, angle = -25, hjust = .1, family = 'Helvetica', face = 'bold'),
        axis.text = element_text(size = 20, face = 'bold', family = 'Helvetica'),
        legend.text = element_text(face = 'bold',family = 'Helvetica', size = 20), legend.title = element_text(family = 'Helvetica', face = 'bold', size = 20, vjust = 1),legend.position = 'right')+
  labs(x = NULL, y = NULL)+coord_flip()

graphics.off()

