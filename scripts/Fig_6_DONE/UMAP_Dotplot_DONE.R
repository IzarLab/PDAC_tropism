# This scripts generates UMAP and dot plots of immune compartment

library(Seurat)
library(SeuratObject)
options(Seurat.object.assay.version = "v4")
library(SeuratWrappers)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)

# Loading data ####

load('compartments/T.RData'); T_ = s_objs
load('compartments/myeloid.RData'); my_ = s_objs
load('compartments/healthy_cancer.RData')
DefaultAssay(T_) = DefaultAssay(my_) = 'RNA'

# UMAP of compartment ####

pdf('immune_landscape.pdf', width = 12, height = 12)

DimPlot(T_, pt.size = 2, group.by = 'cell_type', label = T, label.box = T, label.size = 4, repel = T) & NoLegend()

DimPlot(my_, pt.size = 2, group.by = 'cell_type', label = T, label.box = T, label.size = 4, repel = T) & NoLegend()

# Dot plots ####

# exhaustion receptors' ligands

s_ = subset(s_objs, subset = cell_type %in% c('basal-like','Classical','progenitor','differentiating ductal','E2F Targets (S)','G2-M','ciliated','neural-like progenitor',
                                               'CD197+ cDC','CD1c+ cDC','CD141+CD226+ cDC','pDC',
                                              'M1-like','M2-like','M2a-like','M2c-like','M2d-like',
                                              'mast','monocyte'))
s_$cell_type = factor(x = s_$cell_type, levels = c('basal-like','Classical','progenitor','differentiating ductal','E2F Targets (S)','G2-M','ciliated','neural-like progenitor',
                                                   'CD197+ cDC','CD1c+ cDC','CD141+CD226+ cDC','pDC',
                                                   'M1-like','M2-like','M2a-like','M2c-like','M2d-like',
                                                   'mast','monocyte'))
Idents(s_) = s_$cell_type
DotPlot(object = s_, features = c('LGALS9',
                                  'CD80','CD86',
                                  'NECTIN3','NECTIN2','PVR',
                                  'PDCD1LG2','CD274'),
        scale = F,
        dot.min = 0.05, dot.scale = 15,
        scale.by = 'size',
        cols = c('red','blue'),
        cluster.idents = F, )+
  theme(axis.text.x = element_text(size = 20, angle = -25, hjust = .1, family = 'Helvetica', face = 'bold'),
        axis.text = element_text(size = 20, face = 'bold', family = 'Helvetica'),
        legend.text = element_text(face = 'bold',family = 'Helvetica', size = 20), legend.title = element_text(family = 'Helvetica', face = 'bold', size = 20, vjust = 1),legend.position = 'right')+
  labs(x = NULL, y = NULL)

# Dot plot exhaustion receptors ####

s_ = T_
Idents(s_) = s_$cell_type
DotPlot(object = s_, features = c('HAVCR2','TIGIT','CTLA4','PDCD1'),
        scale = F,
        dot.min = 0.05, dot.scale = 15,
        scale.by = 'size',
        cols = c('red','blue'),
        cluster.idents = F, )+
  theme(axis.text.x = element_text(size = 20, angle = -25, hjust = .1, family = 'Helvetica', face = 'bold'),
        axis.text = element_text(size = 20, face = 'bold', family = 'Helvetica'),
        legend.text = element_text(face = 'bold',family = 'Helvetica', size = 20), legend.title = element_text(family = 'Helvetica', face = 'bold', size = 20, vjust = 1),legend.position = 'right')+
  labs(x = NULL, y = NULL)

graphics.off()

