# This script plots Fig. 5j; it shows the average expression of PDAC-recurrence signatures in exocrine compartment.

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(SeuratDisk)     # Did you 'brew install hdf5'?
library(MatrixGenerics)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggridges)
library(viridis)
library(hacksig)

# Loading exocrine data ####

load(file = 'compartments/exocrine_differentiating-cancer.RData')
DefaultAssay(s_objs) = 'RNA'

# Reading in PDAC recurrence signatures ####

# lung
pdac_lung_markers = read.delim('Misc/cell_types_markers_main_subtypes_lung.tsv', header = T, sep = '\t')
pdac_lung_markers = pdac_lung_markers[pdac_lung_markers$cluster %in% 'all' & pdac_lung_markers$avg_log2FC > 0,]
pdac_lung_markers$gene = trimws(pdac_lung_markers$gene)
pdac_lung_markers = pdac_lung_markers[pdac_lung_markers$gene %in% rownames(s_objs),]
pdac_lung_markers$rank = -log2(pdac_lung_markers$p_val)*pdac_lung_markers$avg_log2FC      # gene ranks
pdac_lung_markers = pdac_lung_markers[order(pdac_lung_markers$rank, decreasing = T),]     # ordering genes based on rank

# liver
pdac_liver_markers = read.delim('Misc/cell_types_markers_main_subtypes_liver.tsv', header = T, sep = '\t')
pdac_liver_markers = pdac_liver_markers[pdac_liver_markers$cluster %in% 'all' & pdac_liver_markers$avg_log2FC > 0,]
pdac_liver_markers$gene = trimws(pdac_liver_markers$gene)
pdac_liver_markers = pdac_liver_markers[pdac_liver_markers$gene %in% rownames(s_objs),]
pdac_liver_markers$rank = -log2(pdac_liver_markers$p_val)*pdac_liver_markers$avg_log2FC       # gene ranks
pdac_liver_markers = pdac_liver_markers[order(pdac_liver_markers$rank, decreasing = T),]      # ordering genes based on rank

# Visualization ####

s_objs@meta.data$avg_lung = colMeans(s_objs[["RNA"]]@data[pdac_lung_markers$gene[1:100],,drop = F])     # top 100 genes
s_objs@meta.data$avg_liver = colMeans(s_objs[["RNA"]]@data[pdac_liver_markers$gene[1:100],,drop = F])

p_ = FeaturePlot(object = s_objs, features = c('avg_liver','avg_lung'),
                 order = T, blend = T,
                 pt.size = .4,
                 coord.fixed = T, raster = F, cols = c('grey90','red3','blue'),
                 min.cutoff = 'q30', max.cutoff = 'q99')
p_ = p_[[3]]+theme(text = element_text(face = 'bold', size = 20, family = 'Helvetica'),
axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
legend.position = 'none')+
labs(title = 'PDAC recurrence signature')

ggsave(plot = p_, device = 'pdf', width = 7.5, height = 7.5, filename = 'j.pdf')

