# This script generates volcano plot of genes DE between PDAC recurrence sites.

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(viridis)

# Reading in adult PDAC recurrence signatures ####

# pdac-lung signature
lung_markers = read.delim('Misc/cell_types_markers_main_subtypes_lung.tsv', header = T, sep = '\t')
lung_markers = lung_markers[lung_markers$cluster %in% 'all' & lung_markers$avg_log2FC > 0, ]
lung_markers$gene = trimws(lung_markers$gene)
lung_markers$rank = -log2(lung_markers$p_val_adj)

# pdac-liver signature
liver_markers = read.delim('Misc/cell_types_markers_main_subtypes_liver.tsv', header = T, sep = '\t')
liver_markers = liver_markers[liver_markers$cluster %in% 'all' & liver_markers$avg_log2FC > 0, ]
liver_markers$gene = trimws(liver_markers$gene)
liver_markers$rank = -log2(liver_markers$p_val_adj)

dt_ = rbind(lung_markers, liver_markers)
dt_$avg_log2FC[dt_$recurrence %in% 'liver'] = - dt_$avg_log2FC[dt_$recurrence %in% 'liver']
dt_$rank[is.infinite(dt_$rank)] = max(dt_$rank[!is.infinite(dt_$rank)])+ 0.01*max(dt_$rank[!is.infinite(dt_$rank)])

# Volcano plot ####

# adding some markers of parenchymal liver-lung cells
idxs_ = which(dt_$gene %in% c(# some adult lung parenchymal markers
                              'CLDN18', 'PGC', 'MUC5AC', 'MUC1', 'CEACAM6',
                              # some adult liver parenchymal markers
                              'HPN', 'SERPING1', 'UGT2B7', 'C3', 'NNMT'))
labels = dt_[idxs_,]

# adding vertical and horizontal lines corresponding to average L2FC's and largest adjusted p-value (q-value)
vlines_ = c(max(dt_$avg_log2FC[dt_$recurrence == 'liver']), min(dt_$avg_log2FC[dt_$recurrence == 'lung']))
hline_ = labels$rank[which.min(labels$rank)]

# plot
ggplot(data = dt_, aes(x = avg_log2FC, y = rank, color = rank))+
theme(panel.background = element_blank(), axis.line = element_line(color = 'black'), panel.grid.major.y = element_line(color = 'grey70'))+
labs(x = 'Average L2FC', y = '-log2(q)', color = '-log2(q)')+
geom_point(shape = 19)+
geom_vline(xintercept = vlines_, color = c('red','blue'), linetype = 'dashed', size = .8)+
geom_hline(yintercept = hline_ - 50, color = 'grey70', linetype = 'dashed', size = .8)+
annotate(geom = 'text', x = min(dt_$avg_log2FC)+.5, y = hline_ - 30, label = format(2^-hline_, scientific = T), fontface = 'bold')+
scale_color_viridis(option = 'H', direction = -1,  breaks = floor(seq(min(dt_$rank), max(dt_$rank),length.out = 3)))+
geom_label_repel(data = labels, aes(label = gene, fill = recurrence), color = 'black', alpha = .75, fontface = 'bold', show.legend = F)+
scale_fill_manual(values = c('red','blue'))

ggsave(file = 'c.pdf', units = 'in', width = 10, height = 7.5, device = 'pdf')

