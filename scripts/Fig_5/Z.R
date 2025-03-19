# This scripts generates a dot plot of PanIN and ductal signatures' score par patient. These two are not comparable.
# This script does not correspond to any figure within the paper but the rebuttal for reviewers.
# rows: signatures
# column: patients
# dot size: fraction of cells significantly expressing a signature
# dot color: signature expression score
# As single-cell suffer from sparsity and imprecise gene expression measurements:
# 1) ALRA wrapped in Seurat package was used to decrease sparsity from 89% to 56%
# 2) an algorithm was used to deal with inconsistent gene expressions
# The output of ALRA:
# The matrix went from 10.58% nonzero to 43.84% nonzero

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratWrappers)
library(ggplot2)
library(viridis)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)
library(parallel)

# Reading in data ####
message('Reading in data')

# signatures

panin_ = c(# PanIN-III
           'MUC5B','MUC3A','EFNA5','LYPD6B','NRG3','MIR31HG','IL1RAPL2','CAPN8','WDR72','NRG1',
           'CRISP3','DUOX2','MOGAT1','ANKRD36C','DMBT1','KRT19','GDA','SHROOM3','KLF5',
           'PCSK5','AHR','RELN','ONECUT2','DLGAP1','KIAA1211',
           # ductal EMT (low-grade PanIN)
           'CREB5','TNFAIP2','ACSM3','CASC15','TGFB2','KCND2','MMP7','NCAM1','SLCO3A1','TSHZ2',
           'FRMD5','NCEH1','NRG3','CSMD2','MIR4435-2HG','SNCAIP','ANXA2','DCLK1','DCDC2','PRAG1',
           'FAM155A','BIRC3','KRT7','UGCG','TNC')

ductal_ = c('SLC4A4','KCNJ15','CADPS','ADAMTS9-AS2','BICC1','PDE3A','PCDH9','CFTR','DZIP1','PKHD1',
            'FRAS1','PBX3','GUCY1A2','AUTS2','ESRRG','KCNJ16','SPP1','PDGFD','ATP1A1','MUC6','CHST9',
            'SLC3A1','VAV3','NRP1','DCDC2','AKAP7','GLIS3','PLEKHS1','POU6F2','LINC00472','C6','CASR',
            'COL18A1','CTTNBP2','ERBB4','SOX9-AS1','TRABD2B','SCTR','SGCZ','MTSS1','FLRT2',
            'LRP1B','PDE7A','MAML2','FAM189A1','SLC17A4','VEPH1','ANXA4','CFAP221','WWTR1')

# cancer data

if(!file.exists('compartments/exocrine_imputed.RData'))
{
  load(file = 'compartments/exocrine.RData')
  s_objs = subset(s_objs, subset = cell_type %in% c('PanIN','ductal/ductal-like','low-grade PanIN'))      # only ductal/ductal-like and PanIN subpopulations are required
  s_objs$cell_type = as.character(s_objs$cell_type)
  s_objs$cell_type[s_objs$cell_type %in% c('PanIN','low-grade PanIN')] = 'PanIN'
  annot_ = s_objs@meta.data
  
  # Imputation ####
  
  s_objs = RunALRA(s_objs, assay = 'RNA', k.only = TRUE)      # only computes optimal k WITHOUT performing ALRA; by default it uses normalized data in slot data
  s_objs = RunALRA(s_objs, assay = 'RNA')                     # runs ALRA with the chosen k
  s_objs = as.matrix(s_objs[['alra']]@data)
  save(s_objs, annot_, file = 'compartments/exocrine_imputed.RData')
}else
{
  load(file = 'compartments/exocrine_imputed.RData')
}

# Correcting gene expressions ####
message('Correcting expression values')

sig_markers = list('PanIN' = panin_, 'ductal/ductal-like' = ductal_)
for(s_ in names(sig_markers))
{
  cells_ = rownames(annot_)[annot_$cell_type %in% s_]
  for(g_ in rownames(s_objs))
  {
    exprs_ = s_objs[g_, cells_]
    
    non_zero_exprs = exprs_ > 0
    if(.10*length(exprs_) <= length(which(non_zero_exprs)))     # first, makes sure the gene is expressed in at least 10% of cells; Default min.pct in FindMarkers() is 1% and we used 9% for DGEA; see the bottom of https://github.com/KlugerLab/ALRA about NCAM1 and CR2
    {
      new_val = quantile(x = exprs_[non_zero_exprs], probs = .5)      # new value (median of non-zero values) to replace inaccurate expression; 
      exprs_[exprs_ < new_val] = new_val
    }
    s_objs[g_,cells_] = exprs_                                  # otherwise use the current measurements
  }
}

# Signature expression score ####
message('Pathway expression score')

# for each patient

annot_[,names(sig_markers)] = 0     # adding pathways as columns to meta data tables
for(p_ in unique(annot_$orig.ident))
{
  # extracting current patient's cancer cells
  cat('processing ', p_, '\n')
  
  # computing pathway (signature) score
  
  cells_ = rownames(annot_)[annot_$orig.ident %in% p_]      # cells of current patient
  
  for(i_ in 1:length(sig_markers))
  {
    pathnm = names(sig_markers)[i_]
    cat('\t> ', pathnm, '\n')
    genes_ = sig_markers[[i_]]      # genes in the current signature
    
    # signal signature expression
    
    pathway_exprs = colMeans(s_objs[genes_,cells_])     # current signal signature mean expression for each cell
    
    # background signature expression
    
    bkgrnd_exprs = list()
    bkgrnd_genes = rownames(s_objs)[!rownames(s_objs) %in% genes_]      # all other genes not present in signal signature comprise background genes
    mat_bgrn = s_objs[bkgrnd_genes, cells_]
    for(e_ in 1:1e3)
    {
      bkgrnd_genes_e = sample(x = bkgrnd_genes, size = length(genes_))      # sample the same number of genes as in signal signature from background genes
      bkgrnd_exprs[[e_]] = colMeans(mat_bgrn[bkgrnd_genes_e,])
    }
    bkgrnd_exprs = do.call(bkgrnd_exprs, what = rbind)                  # 1e3 background signature mean expressions for each cell in current cell type data
    
    # identification of significant cells
    
    probs_ = numeric(length = length(cells_))     # each cell in current cell type data is assigned a p-value
    names(probs_) = cells_
    for(cell_ in cells_)
    {
      probs_[cell_] = round(1- ecdf(bkgrnd_exprs[,cell_])(pathway_exprs[cell_]), 5)     # 1- (fraction/frequency of background mean expressions less than or equal to signal mean expression for current cell)
    }
    sig_cells = names(probs_[probs_ <= 0.01])     # retaining only significant cells (significance level 1%)
    
    # significant signature expression
    
    pathway_exprs = pathway_exprs[sig_cells] - colMeans(bkgrnd_exprs[, sig_cells, drop = F])      # pathway (signal signature) score (avoiding batch effect between cells)
    annot_[sig_cells, pathnm] = pathway_exprs
  }
}

# Creating dot plot ####

# formatting data frame table for ggplot2

dt_ = list(); i_ = 1
for(p_ in unique(annot_$orig.ident))
{
  for(s_ in names(sig_markers))
  {
    exprs_ = annot_[annot_$orig.ident %in% p_, c("condition", s_)]
    dt_[[i_]] = data.frame(patient = p_,
                           pathway = s_,
                           condition = exprs_$condition[1],
                           exprs = mean(exprs_[,s_]),
                           fraction = round(length(which(0 < exprs_[,s_]))/nrow(exprs_)*100, 2))
    i_ = i_ + 1
  }
}
dt_ = do.call(dt_, what = rbind)
dt_$exprs[which.min(dt_$exprs)] = round(min(dt_$exprs),2)
dt_$exprs[which.max(dt_$exprs)] = round(max(dt_$exprs),2)

# plotting

dt_$patient = factor(dt_$patient, levels = c('PN9','PN18',                'PN1','PN10','PN23','PN27','PN3','PN8',
                                             'PN4','PN13','PN15','PN17',  'PN14','PN16','PN19','PN2','PN24','PN28','PN5','PN6','PN7'))
dt_$patient = droplevels(dt_$patient)

ggplot(data = dt_, aes(y = pathway, x = patient, size = fraction, color = exprs, shape = condition))+
theme(panel.background = element_blank(), panel.grid = element_line(color = 'black',linewidth = .05),
      legend.position = 'right', legend.text = element_text(size = 15), legend.key = element_blank(),
      legend.title = element_text(size = 15), legend.spacing.y = unit(.6, 'cm'),
      axis.text.x = element_text(angle = 90),
      text =  element_text(face = 'bold', size = 15, family = 'Helvetica'))+
labs(y = NULL, x = NULL, color = 'Signature score', size = 'Fraction', shape = 'Initial metastatsis')+
geom_point()+
scale_size_continuous(range = c(5,6))+
scale_color_viridis(discrete = F, option = 'H', breaks = c(min(dt_$exprs), round(mean(range(dt_$exprs)),2), max(dt_$exprs)) )+
guides(color = guide_colorbar(ticks.colour = NA, label.position = 'right'),
       size = guide_legend(override.aes = list(size = c(5,5.5,6,6.5))),
       shape = guide_legend(override.aes = list(size = 4:5)))+
coord_fixed(ratio = 3)

ggsave(filename = 'Z.pdf',
       width = 9, height = 5, units = 'in')
