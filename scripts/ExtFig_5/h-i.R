# This script plots Extended Data Fig. 5h-i; it takes primary tumors of adult lung and liver and checks the
# enrichment of PDAC recurrence signatures in them.

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(SeuratDisk)     # Did you 'brew install hdf5'?
library(SeuratWrappers)
library(rstatix)
library(ggplot2)
library(ggridges)
library(patchwork)
library(MatrixGenerics)
library(cowplot)
library(doParallel)
library(parallel)

# Reading in data ####
message('Reading in data')

load(file = 'Misc/merged_human_adult_cancer_lung_liver_processed.RData')
DefaultAssay(lng_liv) = 'RNA'

# Signature score ####
message('signature score')

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

# signature score
# Att.: there should be no scaling at any time before PCA, i.e neither use scale slot nor prcomp(scale = T)

liv_sig = pdac_liver_markers$gene[pdac_liver_markers$avg_log2FC > 0]
lng_sig = pdac_lung_markers$gene[pdac_lung_markers$avg_log2FC > 0]
mat_ = t(as.matrix(lng_liv[["RNA"]]@data[c(liv_sig, lng_sig),]))

cl_ = makeCluster(getOption("cl.cores", 2))
pcs_ = parLapply(X = list(liver = liv_sig, lung = lng_sig), fun = function(sig_, m_)
{
  library(parallel)
  library(doParallel)
  return(prcomp(x = m_[,sig_], rank. = 2))      # scaling is not a good idea
}, m_ = mat_, cl = cl_)
stopCluster(cl_)

liv_score = -pcs_$liver$x[,"PC1"]     # negative projection for PDAC-liver signature; -1 is add to ensure positive values
lng_liv$liv_score = liv_score

lng_score = -pcs_$lung$x[,"PC1"]      # negative projection for PDAC-lung signature; -1 is add to ensure positive values
lng_liv$lng_score = lng_score

# Visualization ####
message('Visualization')

pdf(file = 'h-i.pdf', width = 7.5, height = 7.5)

# scatter plot

Idents(lng_liv) = lng_liv$histology
p_ = FeaturePlot(object = lng_liv, features = c('liv_score','lng_score'),
                 order = T, repel = F,
                 label = T, label.size = 4, label.color = 'black',
                 pt.size = .4,
                 coord.fixed = F, raster = F, cols = c('grey90','red3'),
                 max.cutoff = c('q99','q99'), min.cutoff = c(0,0),
                 ncol = 1)
p_[[1]] = p_[[1]]+scale_color_gradient(low = 'grey90', high = 'red3', breaks = range(p_[[1]]$data$liv_score), labels = c('low','high'))
p_[[2]] = p_[[2]]+scale_color_gradient(low = 'grey90', high = 'red3', breaks = range(p_[[2]]$data$lng_score), labels = c('low','high'))
p_[[1]] = p_[[1]] + theme(legend.position = 'right', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
  labs(title = 'PDAC-liver signature score')
p_[[2]] = p_[[2]] + theme(legend.position = 'right', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
  labs(title = 'PDAC-lung signature score')
plot(p_)

# Violin plot

p_ =  VlnPlot(lng_liv, features = 'liv_score', group.by = "histology", pt.size = 0, combine = FALSE, log = F, sort = F, cols = c('red','blue'))[[1]]+
      theme(legend.position = 'none', axis.text.x = element_text(angle = 90), text = element_text(face = 'bold'))+
      labs(title = 'PDAC-liver signature', y = 'Signature score', x = NULL)+
      stat_summary(fun = median, geom ='point', size = 12, colour = "white", shape = 95)

q_ =  VlnPlot(lng_liv, features = 'lng_score', group.by = "histology", pt.size = 0, combine = FALSE, log = F, sort = F, cols = c('red','blue'))[[1]]+
      theme(legend.position = 'none',  axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90),
            text = element_text(face = 'bold'))+
      labs(title = 'PDAC-lung signature', y = 'Signature score', x = NULL)+
      stat_summary(fun = median, geom='point', size = 12, colour = "white", shape = 95)

r_ = range(range(p_$data$liv_score), range(q_$data$lng_score))
f_ = p_ + scale_y_continuous(limits = r_) + q_ + scale_y_continuous( limits = r_ )
plot(f_)

graphics.off()

# Testing ####
message('Testing')

set.seed(42)
liv_cells = names(lng_liv$cell_type[lng_liv$tissue %in% 'liver'])     # cancer liver cells
lng_cells = names(lng_liv$cell_type[lng_liv$tissue %in% 'lung'])      # cacner lung cells
for(c_ in list(liv_cells, lng_cells))
{
  # Cohen's d (effect size)
  
  dt_ = data.frame(val = c(lng_liv$liv_score[c_],      # cancer liver/lung cells with pdac-liver signature
                           lng_liv$lng_score[c_]),     # cancer liver/lung cells with pdac-lng signature
                   type = rep(x = c('pdac-liver sig','pdac-lung sig'), times = c(length(c_),length(c_))) )
  cat('PDAC recurrence signatures in cells of ', lng_liv$tissue[c_][1],': d = ',
      abs(cohens_d(formula = val ~ type, var.equal = F, data = dt_)$effsize), '\n')
  
  # Wilcox rank-sum test
  
  p_values = list()
  
  for(i_ in 1:1e3)
  {
    c_sub = sample(c_, size = 50, replace = F)
    p_values[[i_]] = wilcox.test(x = lng_liv$liv_score[c_sub],
                                 y = lng_liv$lng_score[c_sub],
                                 alternative = "two.sided")$p.value
  }
  message('PDAC recurrence signatures in cells of ', lng_liv$tissue[c_][1],': q =',
          format(median(p.adjust(p = p_values, method = 'BH')), scien = T))
}

