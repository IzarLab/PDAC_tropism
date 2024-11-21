# This script computes single-sample enrichment scores (SSES) of PDAC recurrence signatures in metastases provided by MET500's paper.
# MET500 cohort is bulk RNA-seq gene expression data of 20 metastatic (not primary!) cancer types, including PDAC.

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
library(parallel)

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

loadings_lng = sort(pcs_$lung$rotation[,"PC1"], decreasing = T)
loadings_lng = names(loadings_lng[loadings_lng > 0])[1:200]      # top 200 genes

loadings_liv = sort(pcs_$liver$rotation[,"PC1"], decreasing = T)
loadings_liv = names(loadings_liv[loadings_liv > 0])[1:200]      # top 200 genes

# Reading in met500 data ####

load('Misc/MET500.RData')
rownames(sampleMeta) = sampleMeta$Sample_id

# converting PDAC-liver/lung-progression gene symbols to Ensembl IDs

loadings_lng = unique(rownames(geneMeta)[geneMeta$gene %in% loadings_lng])
loadings_lng = loadings_lng[loadings_lng %in% rownames(exprsMat)]     # only retains those genes measured in met500

loadings_liv = unique(rownames(geneMeta)[geneMeta$gene %in% loadings_liv])
loadings_liv = loadings_liv[loadings_liv %in% rownames(exprsMat)]     # only retains those genes measured in met500

# sub-matrices of liver and lung ####

lung_smpls = sampleMeta$Sample_id[sampleMeta$biopsy_tissue %in% 'lung']
liver_smpls = sampleMeta$Sample_id[sampleMeta$biopsy_tissue %in% 'liver']
exprsMat = exprsMat[c(loadings_lng, loadings_liv), c(lung_smpls, liver_smpls)]

# Conversion back to HGNC symbols ####

loadings_lng = geneMeta[loadings_lng,'gene']
loadings_liv = geneMeta[loadings_liv,'gene']
rownames(exprsMat) = geneMeta[rownames(exprsMat),'gene']

# Signature single sample scores ####

rslt_ = as.data.frame(hack_sig(expr_data = exprsMat, signatures = list(pdac_lung = loadings_lng, pdac_liver = loadings_liv), method = "singscore", direction = 'up'))
rslt_$metastases = sampleMeta[rslt_$sample_id, "biopsy_tissue"]

rslt_ = as.data.frame(pivot_longer(data = rslt_,
                                   cols = c(pdac_lung, pdac_liver),     # Columns to pivot
                                   names_to = "signature",              # Name of the new column that will contain 'pdac_lung' and 'pdac_liver'
                                   values_to = "score"))                # Name of the new column that will contain the values

# Visaulization ####

# effect size by Cohen's d

d_ = numeric(length = 2)
names(d_) = unique(rslt_$signature)
for(s_ in names(d_))
{
  tmp_ = as.data.frame(rslt_[rslt_$signature %in% s_,])
  d_[s_] = round(abs(cohens_d(data = tmp_, formula = score ~ metastases, var.equal = F)$effsize),2)
}

# plotting and testing by Wilcoxon rank-sum test

pdf(file = 'a.pdf', width = 7.5, height = 7.5, fonts = 'Helvetica')

# scaling between 0-1

rslt_liv = rslt_[rslt_$signature %in% 'pdac_liver',]
rslt_liv$score = (rslt_liv$score - min(rslt_liv$score))/diff(range(rslt_liv$score))

rslt_lng = rslt_[rslt_$signature %in% 'pdac_lung',]
rslt_lng$score = (rslt_lng$score - min(rslt_lng$score))/diff(range(rslt_lng$score))

rslt_ = rbind(rslt_liv, rslt_lng)

rslt_$signature = factor(x = rslt_$signature, levels = c('pdac_lung','pdac_liver'))
rslt_$metastases = factor(x = rslt_$metastases, levels = c('lung','liver'))

# plotting

p_ =  ggplot(data = rslt_, aes(x = signature, y = score, fill = metastases))+
      theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_line(color = 'black'), axis.line = element_line(color = 'black'), panel.background = element_blank())+
      geom_boxplot()+
      geom_jitter(position = position_dodge(.75), color = "black", alpha = .7, size = 1)+
      scale_fill_manual(values = c('blue','red'))+
      annotate(geom = 'text', x = c(1,2), y = c(1.02,1.02), label = paste0('Cohen\'s d = ',d_))+
      stat_compare_means(method = 'wilcox.test', label.y = max(rslt_$score)+0.05, paired = F, method.args = list(alternative = 'two.sided'))      # = wilcox.test(x = a_, y = b_,  paired = FALSE, alternative = "two.sided")$p.value
plot(p_)
graphics.off()

