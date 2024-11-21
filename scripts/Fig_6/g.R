# This script generates a dot plot of cll-cell interactions called by CellPhoneDB for immune compartment

library(ggplot2)
library(gridExtra)
library(cowplot)
library(patchwork)
library(viridis)
library(tidyr)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(Matrix)

# Loading cellphonedb results ####

load(file = 'Misc/CellPhoneDB_cancer_T-NK_myeloid.RData')

sig_means_ = rslt_$significant_means                                                # mean expression of significant interactions
rownames(sig_means_) = sig_means_$interacting_pair
pvals_ = rslt_$pvalues
rownames(pvals_) = pvals_$interacting_pair
pvals_ = pvals_[sig_means_$interacting_pair,]                                       # p-value of significant interactions

sigMeans_info_cols = which(!grepl(pattern = ' vs ', x = colnames(sig_means_)))      # 'vs' is the delimiter of cell types
pvals_info_cols = which(!grepl(pattern = ' vs ', x = colnames(pvals_)))

# Filtering interaction cells ####

pairs_ = colnames(sig_means_)[-sigMeans_info_cols]             # first sigMeans_info_cols columns of sig_means_ are metadata
mate1_ = sub(pattern = ' (vs) .+', replacement = '', x = pairs_, perl = T)
mate2_ = sub(pattern = '.+ (vs) ', replacement = '', x = pairs_, perl = T)

sub1_ = c("NK","CD8+ Tex","CD8+ Tm",   "Treg","CD4+ Tex","effector CD4+ T","CD4+ Tm")
sub2_ = c('cancer',   'M1-like','M2-like','M2a-like','M2c-like','M2d-like','monocyte',    'CD141+CD226+ cDC','CD197+ cDC','CD1c+ cDC','pDC')

cells_ = colnames(sig_means_)[which(mate1_ %in% sub1_ & mate2_ %in% sub2_ |
                                    mate2_ %in% sub1_ & mate1_ %in% sub2_) + max(sigMeans_info_cols)]

# MHCII possibilities!

cells_ = cells_[!cells_ %in% c('cancer vs CD4+ Tex', 'cancer vs Treg',
                               'CD141+CD226+ cDC vs CD8+ Tex','CD141+CD226+ cDC vs NK',
                               'CD197+ cDC vs CD8+ Tex'      ,'CD197+ cDC vs NK',
                               'CD1c+ cDC vs CD8+ Tex'       ,'CD1c+ cDC vs NK',
                               'M1-like vs CD8+ Tex'         ,'M1-like vs NK',
                               'M2-like vs CD8+ Tex'         ,'M2-like vs NK',
                               'M2a-like vs CD8+ Tex'        ,'M2a-like vs NK',
                               'M2c-like vs CD8+ Tex'        ,'M2c-like vs NK',
                               'M2d-like vs CD8+ Tex'        ,'M2d-like vs NK',
                               'monocyte vs CD8+ Tex'        ,'monocyte vs NK')]

# filtering ligand-receptor pairs

intacts_ = rownames(sig_means_)
intacts_ = as.data.frame(do.call(args = lapply(X = intacts_, FUN = function(e_){ e_ = strsplit(x = e_, split = '_', perl = T)[[1]]; return(c(ligand = e_[1], receptor = e_[2]))}), what = rbind))
l_ = c('CD274','PDCD1LG2','PVR','NECTIN2','NECTIN3','CD80','CD86','LGALS9')
r_ = c('PDCD1','TIGIT','CTLA4','HAVCR2')
intacts_ = which(intacts_$ligand %in% l_ & intacts_$receptor %in% r_)

sig_means_ = cbind(sig_means_[intacts_,sigMeans_info_cols], sig_means_[intacts_,cells_])
pvals_ = cbind(pvals_[intacts_,pvals_info_cols], pvals_[intacts_, cells_])      # first pvals_info_cols columns of sig_means_ are metadata

# Correcting p-values ####

qvals_ = as.numeric(as.matrix(pvals_[,-pvals_info_cols]))     # first 11 columns of p-values are meta data; as.numeric() acts column-wise; check this out: as.data.frame(matrix(data = qvals_, byrow = F, nrow = nrow(pvals_), dimnames = list(rownames(pvals_), colnames(pvals_)[-c(1:11)])))
qvals_ = p.adjust(p = qvals_, method = 'BH')
pvals_[,-pvals_info_cols] = as.data.frame(matrix(data = qvals_, nrow = nrow(pvals_), ncol = ncol(pvals_)-max(pvals_info_cols), byrow = F))
sig_means_[,-sigMeans_info_cols][ pvals_[,-pvals_info_cols] > 0.05] = NA

# Removing all-NA rows/columns ####

null_intacts = which(rowSums(sig_means_[,-sigMeans_info_cols], na.rm = T) == 0)     # rowSums with na.rm = T, will return 0 for an all-NA row
if(length(null_intacts) == 0){ null_intacts = -c(1:nrow(sig_means_))}

null_cells = which(colSums(sig_means_[,-sigMeans_info_cols], na.rm = T) == 0)       # any all-NA columns?
if(length(null_cells) == 0){ null_cells = -c(1:ncol(sig_means_[,-sigMeans_info_cols]))}

sig_means_ = sig_means_[-null_intacts, -(null_cells+max(sigMeans_info_cols))]
pvals_ = pvals_[-null_intacts, -(null_cells+max(pvals_info_cols))]

# Plotting ####

sig_means_tmp = sig_means_[, -sigMeans_info_cols]
sig_means_tmp = sig_means_tmp[,sort(colnames(sig_means_tmp))]
sig_means_tmp = cbind(interaction = sig_means_$interacting_pair, sig_means_tmp)
sig_means_tmp = gather(data = sig_means_tmp,
                       2:ncol(sig_means_tmp),     # break these column names
                       key = cells,               # combine them column names into this column
                       value = means)             # store their corresponding values in this column
sig_means_tmp$cells = gsub(x = sig_means_tmp$cells, pattern = ' vs ', replacement = ' <-> ')
sig_means_tmp$cells = factor(x = sig_means_tmp$cells, levels = sort(unique(sig_means_tmp$cells)))

breaks_ = range(sig_means_tmp$means, na.rm = T)
ggplot(data = sig_means_tmp, aes(x = cells, y = interaction))+
  theme(axis.text.x = element_text(face = 'bold', family = 'Helvetica', size = 10, angle = -45, vjust = 2, hjust = -.05),
        axis.text.y = element_text(face = 'bold', family = 'Helvetica', size = 8),
        legend.title = element_text(family = 'Helvetica', face = 'bold', size = 10, vjust = 4), legend.text = element_text(family = 'Helvetica', face = 'bold'),
        panel.background = element_rect(fill = NA),
        panel.grid = element_line(color= 'grey70', linewidth = 0.05))+
  labs(size = 'Mean interaction\nexpression', color = 'Mean interaction\nexpression', x = NULL, y = NULL)+
  geom_point(aes(color = means, size = means), alpha = .8)+
  geom_point(aes(size = means), color = 'black',pch = 21)+
  scale_size_continuous(breaks = sort(c(breaks_, mean(breaks_))), labels = sort(c(breaks_, mean(breaks_))))+
  scale_color_viridis(option = 'H', breaks = c(breaks_, mean(breaks_)))

ggsave(filename = 'g.pdf',width = 12, height = 5.5, device = 'pdf')

