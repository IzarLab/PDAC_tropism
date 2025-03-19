# This script performs regression analysis to regress out the putative effect of copy number variation from log-normalized gene expressions before DGEA.
# It also performs effect size analysis between the CNV matrices of LIV-MR and LUN-MR groups.
# This script does not correspond to any figures within the manuscript but the rebuttal for reviewers.

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratWrappers)
library(rstatix)
library(ggplot2)
library(ggrepel)
library(viridis)

message("\n-------------------------------------------------
        Running differential gene expression analysis corrected for copy number variation &
        Cohen's d effect analysis")

# Reading in liver/lung met PDAC-recurrence signatures ####

pdac_liver_markers = read.delim('Misc/cell_types_markers_main_subtypes_liver.tsv', header = T, sep = '\t')
pdac_liver_markers = pdac_liver_markers[pdac_liver_markers$cluster %in% 'all',]
pdac_liver_markers = pdac_liver_markers[1:500,]                     # top 500 genes; inferCNV for low-quality genes is quite noisy; top 500 high-quality DE genes; try 1000, or all the genes
pdac_liver_markers$gene = trimws(pdac_liver_markers$gene)           # because of Excel issue, genes were padded with one white-space character
pdac_liver_markers$avg_log2FC = -pdac_liver_markers$avg_log2FC      # the original DGEA liver was taken as reference so negative L2FC means genes upregulated in liver singnature

pdac_lung_markers = read.delim('Misc/cell_types_markers_main_subtypes_lung.tsv', header = T, sep = '\t')
pdac_lung_markers = pdac_lung_markers[pdac_lung_markers$cluster %in% 'all' ,]
pdac_lung_markers = pdac_lung_markers[1:500,]
pdac_lung_markers$gene = trimws(pdac_lung_markers$gene)

markers_org = rbind(pdac_lung_markers, pdac_liver_markers)
markers_org$rank = -log2(markers_org$p_val_adj)                     # the rank of genes
markers_org$rank[is.infinite(markers_org$rank)] = max(markers_org$rank[!is.infinite(markers_org$rank)])+ 0.01*max(markers_org$rank[!is.infinite(markers_org$rank)])
mrks_ = c(pdac_liver_markers$gene, pdac_lung_markers$gene)          # combined genes in PDAC-recurrence signatures (LIV-MR and LUN-MR signatures)

# Reading in CNVs ####

load(file = 'Misc/CNV_cancer_compartment.RData')      # loads infercnv.observations.txt (non-reference cells) from all patients; inferCNV generates values in the range of [0,2] cut in 15 intervals; the 8th one with 1 at center corresponds to neutral value
                                                      # not all genes are measured for every patient sample, as many genes could be filtered out by inferCNV during denoising
                                                      # the higher the quality of cell, the more accurate inferCNV
# combined matrix of all patient CNV matrices

mat_ = matrix(data = -1,                                        # CNV matrix of all cancer cells and PDAC-recurrence signature genes; some genes were filtered out by inferCNV during smoothing and thus not measured
              nrow = length(cells_), ncol = length(mrks_),      # rows are all cancer cells, columns genes
              dimnames = list(cells_, mrks_))
for(m_ in obs_mat)                                              # for each sample's matrix of cancer cells (aka non-reference cells)
{
  gs_ = colnames(m_)[colnames(m_) %in% mrks_]     # genes of m_ that are present in the PDAC-recurrence signatures
  m_ = m_[, gs_]
  mat_[ rownames(m_), gs_] = m_
}
rm(obs_mat)

# Denoising ####
# this denoising is different from the inferCNV's internal denoising step
# here, all the values within the neutral interval are replaced with the neutral value (center of the interval)

mat_t = mat_[mat_ != -1]                                                                      # -1 means not measured CNVs
den_ = density(mat_t)                                                                         # density of measured CNVs
neutral_peak_ = den_$x[which.max(den_$y)]                                                     # the value corresponding to CNV neutral is the crest of the density curve
breaks_ = seq(min(mat_t), max(mat_t), length.out = 16)                                        # inferCNV generates values in the range of [0,2] cut in 15 intervals
lower_ = max(which(breaks_ < neutral_peak_))                                                  # the lower bound of the neutral interval
mat_[ (breaks_[lower_] <= mat_ & mat_ <= breaks_[lower_+1]) | mat_ == -1] = neutral_peak_     # all values (including not measured genes) within the neutral range are replaced with the neutral value
mat_[mat_ < neutral_peak_] = -mat_[mat_ < neutral_peak_ ]                                     # negative means deletion/loss

# Loading malignant compartment ####

load('compartments/cancer_subtypes.RData')                          # malignant exocrine compartment of main subtypes

# creating Seurat object
s_ = CreateSeuratObject(counts = s_objs[['RNA']]@data[mrks_,],      # Seurat's FindMarkers uses data slot
                        meta.data = s_objs@meta.data,
                        min.cells = 1, min.features = 1)

# subsetting mat_ on cancer cells of main subtypes
cells_ = colnames(s_objs)                                           # cancer cells
cells_ = sub(x = cells_, pattern = '_[0-9]+', replacement = '')     # _[0-9][0-9] is added to barcodes by Seurat
cells_ = paste0(s_objs$orig.ident,'_', cells_)                      # inferCNV precedes barcodes with patient IDs like PN10_GCACTCTTCCCACTTG-1
mat_ = t(mat_[cells_, ]); rm(s_objs); gc()

# Regressing out copy number variation --------------------------------------------------------------------------------

expr_mat = as.matrix(s_[['RNA']]$data)                              # converting sparse expression matrix of Seurat object to dense matrix
zeros_ = which(expr_mat == 0)                                       # zero-valued genes means they were not measured by snRNA-seq; they MUST remain zero
adjusted_expr = list()
r_ = NULL                                                           # Pearson correlation coefficient between observed and predicted gene expressions
for(i_ in 1:nrow(expr_mat))
{
  lm_ = lm(expr_mat[i_,] ~ mat_[i_,])
  r_ = c(r_, sqrt(summary(lm_)$r.squared))
  adjusted_expr[[i_]] = residuals(lm_)      # regressing out
}
adjusted_expr = do.call(adjusted_expr, what = rbind)
rownames(adjusted_expr) = rownames(mat_)
adjusted_expr[zeros_] = 0                                           # zero-valued genes 
s_[['RNA']]$data = as(adjusted_expr, "dgCMatrix")                   # adding adjusted expression matrix to Seurat object
message('Average correlation between observed and predicted gene expression is ', round(mean(r_),2))

## DGEA ####

markers_ = FindMarkers(object = s_,
                       group.by = 'condition', ident.1 = 'lung', ident.2 = 'liver',
                       min.pct = 0.09, logfc.threshold = .40)
markers_$cluster = 'all'
markers_$gene = rownames(markers_)
pos_inds = 0 < markers_$avg_log2FC
m_p = markers_[pos_inds,]                             # contains genes upregulated in pdac-lung patients
m_p = m_p[order(m_p$avg_log2FC, decreasing = T),]     # ordering by L2FC
m_p$recurrence = 'lung'

m_n = markers_[!pos_inds,]                            # contains genes upregulated in pdac-liver patients
m_n = m_n[order(m_n$avg_log2FC, decreasing = F),]
m_n$recurrence = 'liver'

markers_ = rbind(m_p, m_n)
markers_$rank = -log2(markers_$p_val_adj)
markers_$rank[is.infinite(markers_$rank)] = max(markers_$rank[!is.infinite(markers_$rank)])+ 0.01*max(markers_$rank[!is.infinite(markers_$rank)])

## Statics ####

t1_ = length(pdac_liver_markers$gene[pdac_liver_markers$gene %in% m_n$gene])/length(pdac_liver_markers$gene)*100
t2_ = length(pdac_liver_markers$gene[pdac_liver_markers$gene %in% m_p$gene])
message('After adjustment, ',t1_,'% of LIN-MR signature was recovered, and ', t2_, ' LIV-MR genes were found upregualted in the lung group')

t1_ = length(pdac_lung_markers$gene[pdac_lung_markers$gene %in% m_p$gene])/length(pdac_lung_markers$gene)*100
t2_ = length(pdac_lung_markers$gene[pdac_lung_markers$gene %in% m_n$gene])
message('After adjustment, ',t1_,'% of LUN-MR signature was recovered, and ', t2_, ' LUN-MR genes were found upregualted in the liver group')

## Volcano plot ####

p_ = list(); i_ = 1# list of plots
for(m_ in list(markers_org, markers_))
{
  # adding some markers of parenchymal liver-lung cells
  idxs_ = which(m_$gene %in% c(# some adult lung parenchymal markers
                               'CLDN18', 'PGC', 'MUC5AC', 'MUC1', 'CEACAM6',
                               # some adult liver parenchymal markers
                               'HPN', 'SERPING1', 'UGT2B7', 'C3', 'NNMT'))
  labels = m_[idxs_,]
  
  # adding vertical and horizontal lines corresponding to average L2FC's and largest adjusted p-value (q-value)
  vlines_ = c(max(m_$avg_log2FC[m_$recurrence == 'liver']), min(m_$avg_log2FC[m_$recurrence == 'lung']))
  hline_ = labels$rank[which.min(labels$rank)]
  
  # plot
  t_ = if(i_ == 1){ 'DGEA with CNV' }else{ 'DGEA without CNV' }
  p_[[i_]] =  ggplot(data = m_, aes(x = avg_log2FC, y = rank, color = rank, shape = recurrence))+
              theme(panel.background = element_blank(), axis.line = element_line(color = 'black'), panel.grid.major.y = element_line(color = 'grey70'),
                    text = element_text(face = 'bold', size = 20), plot.title = element_text(hjust = 0.5, vjust = 0.5), legend.title = element_text(vjust = 3))+
              labs(x = 'Average L2FC', y = '-log2(q)', color = '-log2(q)', shape = 'Recurrence', title = t_)+
              geom_point(color = 'black', size = 2)+geom_point(size = 1)+
              geom_vline(xintercept = vlines_, color = c('red','blue'), linetype = 'dashed', linewidth = .8)+
              geom_hline(yintercept = hline_ - 10, color = 'grey70', linetype = 'dashed', linewidth = .8)+
              annotate(geom = 'text', x = min(m_$avg_log2FC)+2, y = hline_ +10, size = 5, label = paste0('q = ',format(signif(2^-hline_,digits = 2), scientific = T)), fontface = 'bold')+
              scale_color_viridis(option = 'inferno', direction = -1,  breaks = floor(seq(min(m_$rank), max(m_$rank),length.out = 3)))+
              geom_label_repel(data = labels, aes(label = gene, fill = recurrence), color = 'black', box.padding = .5, label.padding = .2, alpha = .85, fontface = 'bold', show.legend = F)+
              scale_fill_manual(values = c('red','skyblue2'))
  i_ = i_ + 1
}

p_ = p_[[1]]+p_[[2]]
ggsave(plot = p_, file = 'DGEA_CNV_regression.pdf', units = 'in', width = 15, height = 7.5, device = 'pdf')

# Differential gene copy number analysis #--------------------------------------------------------------------------------

liv_mat = mat_[,cells_[s_$condition %in% 'liver']]      # CNV matrix of pdac-liver cells
lng_mat = mat_[,cells_[s_$condition %in% 'lung']]       # CNV matrix of pdac-lung cells

effsiz_ = list()                                        # a list to store Cohen's d effect size
for(i_ in 1:nrow(mat_))
{
  dt_ = data.frame(val = c(liv_mat[i_,],
                           lng_mat[i_,]),
                   type = rep(x = c('pdac-liver cells',
                                    'pdac-lung cells'),
                              times = c(ncol(liv_mat),
                                        ncol(lng_mat))) )
  effsiz_[[mrks_[i_]]] = as.data.frame(cohens_d(formula = val ~ type, var.equal = F, data = dt_)[,c('effsize','magnitude')])      # Cohen's d
}
effsiz_ = do.call(effsiz_, what = rbind)

constant_vals = is.nan(effsiz_$effsize)                 # cohens_d() generates NaN for effect size of constant values (neutral CNV); these are set to 0 (negligible)
effsiz_$effsize[constant_vals] = 0
effsiz_$magnitude[constant_vals] = 'negligible'

effsiz_$effsize = abs(effsiz_$effsize)                  # effect size could be negative, while only its magnitude matters
stats_ = sort(table(effsiz_$magnitude))

message("\nDifferential gene copy number analysis: Cohen's d (effect size):")
print(stats_)

## Volcano plot ####

effsiz_$gene = rownames(effsiz_)
mean_lng = rowMeans(lng_mat)
mean_liv = rowMeans(liv_mat)
effsiz_$mean_diff = mean_lng - mean_liv     # rownames of effsiz_, lng_mat and liv_mat are the same
effsiz_[pdac_lung_markers$gene,'recurrence'] = 'lung'
effsiz_[pdac_liver_markers$gene,'recurrence'] = 'liver'

# adding some markers of parenchymal liver/lung cells
idxs_ = which(effsiz_$gene %in% c(# some adult lung parenchymal markers
                                  'CLDN18', 'PGC', 'MUC5AC', 'MUC1', 'CEACAM6',
                                  # some adult liver parenchymal markers
                                  'HPN', 'SERPING1', 'UGT2B7', 'C3', 'NNMT'))
labels = effsiz_[idxs_,]

# coloration
cutoff_ = 0.8                                                                     # cohens_d(): based on (Cohen 1992), |d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
col_bins = seq(min(effsiz_$effsize), max(effsiz_$effsize), length.out = 100)      # color bins: effect sizes are binned by 100
cols_ = inferno(n = 100)                                                          # 100 colors from inferno color map, one assigned to each color bin
names(cols_) = col_bins
cols_[col_bins < cutoff_] = 'grey70'                                              # effect sizes with magnitudes of small to moderate are set to gray (insignificant)
effsiz_$color = NA                                                                # final color of each gene on volcano plot
for(i_ in 1:nrow(effsiz_))
{
  indx_ = min(which(effsiz_$effsize[i_] <= col_bins))                             # into which bin does the current effect size fall?
  effsiz_$color[i_] = cols_[indx_]                                                # color of the bin is assigned to the gene
}

# plot
p_ =  ggplot(data = effsiz_, aes(x = mean_diff, y = effsize, shape = recurrence))+
      theme(panel.background = element_blank(), axis.line = element_line(color = 'black'), panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
            text = element_text(face = 'bold', size = 20), plot.title = element_text(vjust = 0.5, hjust = 0.5))+
      labs(x = 'Mean difference (lung vs liver)', y = 'Effect size (d)', title = 'Gene copy number', shape = 'Recurrence')+
      geom_point(color = 'black', size = 2)+geom_point(color = effsiz_$color, size = 1)+
      coord_trans(y = 'log1p', x = 'log1p')+     # does not transform value but only visualization
      geom_hline(yintercept = cutoff_, color = 'grey70', linetype = 'dashed', linewidth = .8)+
      annotate(geom = 'text', x = min(effsiz_$mean_diff)+0.12, y = cutoff_ + 0.02, size = 5, label = 'min large effect size', color = 'red3', fontface = 'bold')+
      geom_label_repel(data = labels, aes(label = gene, fill = recurrence), color = 'black', box.padding = .5, label.padding = .2, alpha = .75, fontface = 'bold', show.legend = F)+
      scale_fill_manual(values = c('red','skyblue2'))

ggsave(file = 'CNV_DA.pdf', units = 'in', width = 10, height = 7.5, device = 'pdf')
