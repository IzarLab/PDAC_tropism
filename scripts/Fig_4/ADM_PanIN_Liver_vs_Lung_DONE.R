# This script plots UMAP projection of the expression of PDAC-recurrence signatures in exocrine compartment

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

load(file = 'compartments/exocrine_differentiating-cancer.RData')

# Metastatic PDAC liver/lung tumor signatures ####

liver_ = c('HPN','SERPING1','UGT2B7','C1S','CYP4F11','C2','C3',
           'CFB','NNMT','SLCO1B3','FAH','XDH','SLC16A2','PECR')
lung_ = c('PGC','KCNJ15','ROS1','CXCL17','SLC26A9','CLDN18','ANKRD29','VSIG2',
          'LAMA3','CEACAM6','ROR1','C4orf19','RTKN2','ATP13A4','CAPN8','TMPRSS2')

genes = liver_
genes = genes[genes %in% rownames(s_objs)]
s_objs@meta.data$avg_liver = colMeans(s_objs[["RNA"]]@data[genes,,drop = F])

genes = lung_
genes = genes[genes %in% rownames(s_objs)]
s_objs@meta.data$avg_lung = colMeans(s_objs[["RNA"]]@data[genes,,drop = F])

s_ = s_objs

# Visualization ####

p_ = FeaturePlot(object = s_, features = c('avg_liver','avg_lung'),
                 order = T, repel = F, blend = T,
                 label = F, label.size = 6, label.color = 'green4',
                 pt.size = 1,
                 coord.fixed = F, raster = F, cols = c('grey90','red3','blue'),
                 max.cutoff = 'q99', min.cutoff = 'q40',
                 ncol = 1)
p_ = p_[[3]]+theme(text = element_text(face = 'bold', size = 20, family = 'Helvetica'),
axis.title = element_blank(),axis.text = element_blank(), axis.ticks = element_blank(),
legend.position = 'none')+
labs(title = 'PDAC-recurrence signature')

ggsave(plot = p_, device = 'pdf', width = 7, height = 7, filename = 'PDAC_recurrence_signature_exo_cancer.pdf')

