# This script plots Extended Data Fig. 3b-c; it finds top n informative genes from PDAC-recurrence signatures that separate out lung and liver single cells
# in healthy adult lung-liver single-cell data. To find top n, it runs PCA on each signature and extracts genes with positive-direction rotations.
# This scripts is used with Fig_4/a.R and Fig_4/g-h.R scripts to check the informativeness of top n genes for MET500 and CCLE-MetMap500 datasets.

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(SeuratDisk)# Did you 'brew install hdf5'?
library(MatrixGenerics)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggridges)
library(viridis)
library(hacksig)

# Gene score ####
message('Gene score')

# reading in healthy adult liver-lung single cell data

load(file = 'Misc/merged_human_adult_lung_liver_processed.RData')

# reading in liver/lung met PDAC signatures

pdac_lung_markers = read.delim('Misc/cell_types_markers_main_subtypes_lung.tsv', header = T, sep = '\t')
pdac_lung_markers = pdac_lung_markers[pdac_lung_markers$cluster %in% 'all',]
pdac_lung_markers$gene = trimws(pdac_lung_markers$gene)
pdac_lung_markers = pdac_lung_markers[pdac_lung_markers$gene %in% rownames(lng_liv),]

pdac_liver_markers = read.delim('Misc/cell_types_markers_main_subtypes_liver.tsv', header = T, sep = '\t')
pdac_liver_markers = pdac_liver_markers[pdac_liver_markers$cluster %in% 'all',]
pdac_liver_markers$gene = trimws(pdac_liver_markers$gene)
pdac_liver_markers = pdac_liver_markers[pdac_liver_markers$gene %in% rownames(lng_liv),]

# gene scores (loadings)

lng_sig = pdac_lung_markers$gene[pdac_lung_markers$avg_log2FC > 0]
liv_sig = pdac_liver_markers$gene[pdac_liver_markers$avg_log2FC > 0]
mat_ = t(as.matrix(lng_liv[["RNA"]]@data[c(liv_sig, lng_sig),]))

cl_ = makeCluster(getOption("cl.cores", 2))
pcs_ = parLapply(X = list(liver = liv_sig, lung = lng_sig), fun = function(sig_, m_)
                                                                  {
                                                                    library(parallel)
                                                                    library(doParallel)
                                                                    return(prcomp(x = m_[,sig_], rank. = 2))
                                                                  }, m_ = mat_, cl = cl_)
stopCluster(cl_)

loadings_liv = sort(pcs_$liver$rotation[,"PC1"], decreasing = T)
loadings_liv = names(loadings_liv[loadings_liv > 0])[1:200]     # top 200 genes

loadings_lng = sort(pcs_$lung$rotation[,"PC1"], decreasing = T)
loadings_lng = names(loadings_lng[loadings_lng > 0])[1:200]     # top 200 genes

# visualization

pdf(file = 'b-c.pdf', width = 7.5, height = 15, fonts = 'Helvetica')

lng_liv$liv_loadings = colMeans(lng_liv[["RNA"]]@data[loadings_liv,])
lng_liv$lng_loadings = colMeans(lng_liv[["RNA"]]@data[loadings_lng,])
p_ = FeaturePlot(lng_liv, features = c('liv_loadings','lng_loadings'),
                 pt.size = .1, order = T,ncol = 1,
                 repel = T, label = T, label.size = 3,
                 cols = c('grey90','red3'),
                 min.cutoff = c(NA,NA), max.cutoff = c(NA,NA),
                 raster = F)
p_[[1]] = p_[[1]] + labs(title = "PDAC-liver signature's\naverage laoding")
p_[[2]] = p_[[2]] + labs(title = "PDAC-lung signature's\naverage laoding")
plot(p_)

graphics.off()
