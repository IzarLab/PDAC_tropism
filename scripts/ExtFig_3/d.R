# This script takes adult matched primary and metastatic PDACs in liver and lymphnode (PMID: 39294496) and checks the enrichment of PDAC-liver and PDAC-lung signatures in them.

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(SeuratDisk)# Did you 'brew install hdf5'?
library(SeuratWrappers)
library(rstatix)
library(ggplot2)
library(ggridges)
library(patchwork)
library(MatrixGenerics)
library(cowplot)
set.seed(42)

# Reading in data ####
message('Reading in data')

pdac_ = readRDS('Misc/spatial_matched_PDAC_liver_lymphnode.rds')      # there are 91,496 spots; 30 histopathological images and no Visium image!
DefaultAssay(pdac_) = 'SCT'                                           # normalized Spatial is in data slot of SCT assay as they used SCTransform() rather than NormalizeData()

cells_ = colnames(pdac_)[pdac_$Histology %in% c('PDAC','Liver_Mets', 'Lymph node') &
                         pdac_$spot_class != 'reject' &
                         !(pdac_$first_type %in% 'Normal Epithelial cells' | pdac_$second_type %in% 'Normal Epithelial cells') ]

pdac_ = subset(pdac_, cells = cells_)
pdac_$first_type = droplevels(pdac_$first_type)
pdac_$second_type = droplevels(pdac_$second_type)

DimPlot(object = pdac_, reduction = 'umap.ischia14', group.by = 'Histology')+      # there are several umap reductions and this one is the default if default assay is rctd_fullfinal!
DimPlot(object = pdac_, reduction = 'umap.ischia14', group.by = 'first_type')+
DimPlot(object = pdac_, reduction = 'umap.ischia14', group.by = 'second_type')+
DimPlot(object = pdac_, reduction = 'umap.ischia14', group.by = 'spot_class')

# ggsave(filename = 'landscape.pdf', device = 'pdf', width = 15, height = 15)

# Signature score ####
message('signature score')

# reading in liver/lung met PDAC signatures

pdac_liver_markers = read.delim('Misc/cell_types_markers_main_subtypes_liver.tsv', header = T, sep = '\t')
pdac_liver_markers = pdac_liver_markers[pdac_liver_markers$cluster %in% 'all',]
pdac_liver_markers = trimws(pdac_liver_markers$gene)
pdac_liver_markers = pdac_liver_markers[pdac_liver_markers %in% rownames(pdac_)][1:100]      # Att.: rownames(pdac_)? Default assay is SCT

pdac_lung_markers = read.delim('Misc/cell_types_markers_main_subtypes_lung.tsv', header = T, sep = '\t')
pdac_lung_markers = pdac_lung_markers[pdac_lung_markers$cluster %in% 'all',]
pdac_lung_markers = trimws(pdac_lung_markers$gene)
pdac_lung_markers = pdac_lung_markers[pdac_lung_markers %in% rownames(pdac_)][1:100]      # Att.: rownames(pdac_)? Default assay is SCT

# PDAC-recurrence signature score

Idents(pdac_) = pdac_$Histology
pdac_$liv_score = colMeans(pdac_[["SCT"]]@data[pdac_liver_markers,])
pdac_$lng_score = colMeans(pdac_[["SCT"]]@data[pdac_lung_markers,])

dt_ = data.frame(score = c(pdac_$liv_score, pdac_$lng_score), Histology = c(pdac_$Histology, pdac_$Histology), signature = rep(c('Liver','Lung'), times = c(ncol(pdac_),ncol(pdac_))) )

# Cohen's d and Wilcoxon rank-sum test ###

dt_$score = sqrt(dt_$score)
d_ = p_values = numeric(length = 3)
names(d_) = names(p_values) = unique(dt_$Histology)
for(h_ in unique(dt_$Histology))
{
  # Cohen's d
  
  tmp_ = as.data.frame(dt_[dt_$Histology %in% h_,])
  d_[h_] = round(abs(cohens_d(data = tmp_, formula = score ~ signature, var.equal = F)$effsize),2)
  
  # p-value
  
  tmp_liv = tmp_[tmp_$signature %in% 'Liver',]      # tmp_liv and tmp_lng have the same dimensions
  tmp_lng = tmp_[tmp_$signature %in% 'Lung',]
  p_ls = list()       # list of 1000 p-values
  for(i_ in 1:1e3)
  {
    c_sub = sample(x = nrow(tmp_liv), size = 50, replace = F)     # tmp_liv and tmp_lng have the same dimensions
    p_ls[[i_]] = wilcox.test(x = tmp_liv$score[c_sub],
                                 y = tmp_lng$score[c_sub],
                                 alternative = "two.sided")$p.value
  }
  p_values[h_] = format(median(p.adjust(p = p_ls, method = 'BH')), scien = T)
}

# Visualization ####

ggplot(data = dt_, aes(x = Histology, y = score, fill = signature))+
theme(panel.background = element_blank(), axis.line = element_line(color = 'black'),
      text = element_text(face = 'bold', size = 20, color = 'black'),
      plot.title = element_text(face = 'bold', size = 20, hjust = 0.5))+
labs(fill = 'PDAC recurrence\nsignature', title = 'Matched primary and metastatic PDAC samples')+
geom_violin()+
stat_summary(fun = median, geom = "crossbar", 
             width = 0.16, color = "white",
             position = position_dodge(0.9), show.legend = F)+
annotate(geom = 'text', x = c(1,2,3), y = max(dt_$score)+0.05, label = paste0('Cohen\'s d = ',d_), fontface = 'bold', size = 5)+
annotate(geom = 'text', x = c(1,2,3), y = max(dt_$score)+0.01, label = paste0('Wilcoxon q = ',p_values), fontface = 'bold', size = 5)+
scale_fill_manual(values = c('red','blue'))

ggsave(filename = 'd.pdf', device = 'pdf', width = 14, height = 8, )

