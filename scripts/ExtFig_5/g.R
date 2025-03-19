# This script plots Extended Data Fig. 5g; it computes single-sample gene set enrichment (ssGSE) scores of
# PDAC-recurrencesignatures in CCLE-MetMap500 cohort.
# CCLE-MetMap500's cohort is bulk RNAseq data of cell lines from various cancer types including PDAC.

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

# genes to remove
genes_remove = read.table(file = 'Misc/gdx/GDX_gene.tsv', header = T, sep = '\t', quote = "", as.is = T, check.names = F)
rownames(genes_remove) = genes_remove$gene

pdac_liver_markers = read.delim('Misc/cell_types_markers_main_subtypes_liver.tsv', header = T, sep = '\t')
pdac_liver_markers = pdac_liver_markers[pdac_liver_markers$cluster %in% 'all',]
pdac_liver_markers$gene = trimws(pdac_liver_markers$gene)
pdac_liver_markers = pdac_liver_markers[pdac_liver_markers$gene %in% rownames(lng_liv),]
pdac_liver_markers = pdac_liver_markers[!pdac_liver_markers$gene %in% rownames(genes_remove),]

pdac_lung_markers = read.delim('Misc/cell_types_markers_main_subtypes_lung.tsv', header = T, sep = '\t')
pdac_lung_markers = pdac_lung_markers[pdac_lung_markers$cluster %in% 'all',]
pdac_lung_markers$gene = trimws(pdac_lung_markers$gene)
pdac_lung_markers = pdac_lung_markers[pdac_lung_markers$gene %in% rownames(lng_liv),]
pdac_lung_markers = pdac_lung_markers[!pdac_lung_markers$gene %in% rownames(genes_remove),]

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
loadings_lng = names(loadings_lng[loadings_lng > 0])[1:200]     # top 200 lineage-oriented genes; for pdac cell line you can use: loadings_lng = pdac_lung_markers$gene[1:100]

loadings_liv = sort(pcs_$liver$rotation[,"PC1"], decreasing = T)
loadings_liv = names(loadings_liv[loadings_liv > 0])[1:200]     # top 200 lineage-oriented genes; for pdac cell line you can use: loadings_liv = pdac_liver_markers$gene[1:100]

# Reading in MetMap500 data ####

# tropism data
lung_tropism = read.delim(file = 'Misc/lung_metastasis_potential.txt', sep = '\t', header = T)
rownames(lung_tropism) = lung_tropism$CCLE_name

liver_tropism = read.delim(file = 'Misc/liver_metastasis_potential.txt', sep = '\t', header = T)
rownames(liver_tropism) = liver_tropism$CCLE_name

# adding preference column
tropism_ = data.frame(lung_mean = lung_tropism$mean, liver_mean = liver_tropism$mean, preference = 'Liver', row.names = rownames(lung_tropism))     # lung_tropism and liver_tropism have identical row names (CCLE_name)
tropism_$preference[tropism_$liver_mean < tropism_$lung_mean] = 'Lung'

# meta data file
sampleMeta = read.delim(file = 'Misc/MetMap_cell_line_annotation.txt', sep = '\t', header = T)
rownames(sampleMeta) = sampleMeta$depmap_id
sampleMeta = sampleMeta[sampleMeta$corrected_site_of_origin %in% 'primary',]      # only cell lines taken form primary tumor

# expression matrix
exprsMat = t(read.delim(file = 'Misc/OmicsExpressionProteinCodingGenesTPMLogp1.csv', sep = ',', check.names = F, as.is = T, header = T, row.names = 1))
rownames(exprsMat) = sub(x = rownames(exprsMat), pattern = ' \\(.+', replacement = '')

# Cancer expression matrices ####

loadings_lng = loadings_lng[loadings_lng %in% rownames(exprsMat)]
loadings_liv = loadings_liv[loadings_liv %in% rownames(exprsMat)]

# filtering cancer cell lines with non-metastatic, weakly metastatic, predominantely metastatic to one site or few cell clines of a cancer type
lines_ = sampleMeta[sampleMeta$cancer_type %in% c('esophageal cancer','head and neck cancer','colorectal cancer',
                                                  'pancreatic cancer',
                                                  'endometrial cancer','breast cancer','thyroid cancer'),
                    c("CCLE_name","depmap_id")]

lines_$preference = tropism_[lines_$CCLE_name, "preference"]
lines_ = lines_[! is.na(lines_$preference),]
mat_ = exprsMat[c(loadings_lng, loadings_liv),lines_$depmap_id]

# Signature single sample scores ####

rslt_ = as.data.frame(hack_sig(expr_data = mat_, signatures = list(pdac_lung = loadings_lng, pdac_liver = loadings_liv), method = "ssgsea", direction = 'up'))
rslt_$tropism = lines_[rslt_$sample_id, "preference"]

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
  d_[s_] = round(abs(cohens_d(data = tmp_, formula = score ~ tropism, var.equal = F)$effsize),2)
}

# plotting and testing by Wilcoxon rank-sum test

pdf(file = 'g.pdf', width = 7.5, height = 7.5, fonts = 'Helvetica')

# scaling between 0-1

rslt_liv = rslt_[rslt_$signature %in% 'pdac_liver',]
rslt_liv$score = (rslt_liv$score - min(rslt_liv$score))/diff(range(rslt_liv$score))

rslt_lng = rslt_[rslt_$signature %in% 'pdac_lung',]
rslt_lng$score = (rslt_lng$score - min(rslt_lng$score))/diff(range(rslt_lng$score))

rslt_ = rbind(rslt_liv, rslt_lng)

rslt_$signature = factor(x = rslt_$signature, levels = c('pdac_lung','pdac_liver'))
rslt_$tropism = factor(x = rslt_$tropism, levels = c('Lung','Liver'))

# plotting

p_ =  ggplot(data = rslt_, aes(x = signature, y = score, fill = tropism))+
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_line(color = 'black'), axis.line = element_line(color = 'black'), panel.background = element_blank())+
  geom_boxplot()+
  stat_summary(fun = median, position = position_dodge(.75), geom ='point', size = 15, colour = "white", shape = 95)+
  # geom_jitter(position = position_dodge(.75), color = "black", alpha = 1)+
  scale_fill_manual(values = c('blue','red'))+
  stat_compare_means(method = 'wilcox.test', label.y = max(rslt_$score)+0.05, paired = F, method.args = list(alternative = 'two.sided'))+      # = wilcox.test(x = a_, y = b_,  paired = FALSE, alternative = "two.sided")$p.value
  annotate(geom = 'text', x = c(1,2), y = c(1.02,1.02), label = paste0('Cohen\'s d = ',d_))
plot(p_)

graphics.off()
