# This script plots CNV heatmap for all patients in cancer compartment

library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)
library(pheatmap)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratWrappers)
library(patchwork)
library(cowplot)

# Reading in CNVs ####

load(file = 'Misc/CNV_cancer_compartment.RData')

# Reading in ordered genes file ####

genes_chr = read.delim(file = 'Misc/gene_ordering_file.txt',                                                # ordered gene-chromosome file; MUST BE THE ONE USED BY inferCNV
                       header = F, sep = '\t', quote = "", as.is = T, check.names = F, col.names = c('gene','chr','str','end'))
rownames(genes_chr) = genes_chr$gene
genes_ = genes_chr$gene[genes_chr$gene %in% genes_]                                                         # only genes measured by every sample in the cancer compartment; see above
names(genes_) = genes_chr[genes_,"chr"]                                                                     # adding chromosome names to each gene

# Forming universal obs mat for pheatmap() ####

mat_ = matrix(data = 1, nrow = length(cells_), ncol = length(genes_), dimnames = list(cells_, genes_))      # CNV matrix of all cancer cells and genes measured 
for(m_ in obs_mat)                                                                                          # for each sample's matrix of non-reference cancer cells
{
  m_ = m_[,genes_]      # genes_ are present in every sample; one heatmap for all patients!
  mat_[ rownames(m_), genes_] = m_
}
rm(obs_mat); gc()

# Preparing annotation table for pheatmap() ####

annot_$cell_type = as.character(annot_$cell_type)
annot_$cell_type[annot_$cell_type %in% c('E2F Targets (S)','G2-M')] = 'mitotic'
annot_$cell_type = factor(annot_$cell_type, levels = c("neural-like progenitor",
                                                       "ciliated",
                                                       "mitotic",
                                                       "Classical",
                                                       "basal-like",
                                                       "progenitor",
                                                       "differentiating ductal"))

annoCol = list(cell_type = c("neural-like progenitor" = 'purple',
                             "ciliated"               = 'red',
                             "mitotic"                = 'orange',
                             "Classical"              = 'green',
                             "basal-like"             = 'pink4',
                             "progenitor"             = 'brown',
                             "differentiating ductal" = 'yellow3'),
               condition = c('liver' = 'red', 'lung' = 'blue'))

# Ordering cancer cell types for every sample in the merged meta data table (annot_) ####

a_ = list()
cells_ = NULL
for(s_ in unique(annot_$orig.ident))
{
  a_[[s_]] = annot_[annot_$orig.ident %in% s_,]
  a_[[s_]] = a_[[s_]][order(a_[[s_]]$cell_type),]
  cells_ = c(cells_, rownames(a_[[s_]]) )
}
annot_ = do.call(a_, what = rbind); rm(a_); gc()
rownames(annot_) = cells_
annot_ = annot_[, -1]                                                                                       # first column is orig.ident
mat_ = mat_[rownames(annot_),]

# Breaks ####

den_ = density(mat_)
neutral_peak_ = den_$x[which.max(den_$y)]                                                                   # neutral CNV value

breaks_ = seq(min(mat_), max(mat_), length.out = 16)                                                        # inferCNV has 15 intervals
lower_ = max(which(breaks_ < neutral_peak_))                                                                # neutral interval is here
dels_ = colorRampPalette(colors =  c('blue4','white'))(8)[max(7-(lower_-2),1):7]                            # based on which interval in neutral, some of the 8 blue shades are selected
neutral_ = colorRampPalette(colors =  c('white'))(1)
amps_ = colorRampPalette(colors =  c('white','red4'))(8)[2:(1+min(7,15-lower_))]                            # based on which interval in neutral, some of the 8 red shades are selected
cols_ = c(dels_, neutral_, amps_)

# Row labels and gaps ####

row_labs = rep('',nrow(mat_))
smpls_ = unlist(strsplit(x = rownames(mat_), split = '[_]'))
smpls_ = smpls_[seq(1,length(smpls_), by = 2)]
smpls_idx = sapply(X = unique(smpls_), FUN = function(s_){ max(which(smpls_ %in% s_))})                     # end position of each sample on the heatmap
row_labs[smpls_idx] = names(smpls_idx)
row_gaps = smpls_idx

# Col labels and gaps ####

gaps_ = sapply(X = unique(names(genes_)), FUN = function(chr_){ max(which(names(genes_) %in% chr_)) })      # end position of chromosomes on the heatmap
chr_index = sapply(X = unique(names(genes_)), FUN = function(chr_){ mean(which(names(genes_) %in% chr_)) })
collabs = rep('', ncol(mat_))
collabs[chr_index] = names(chr_index)

# Heatmap ####

p_ = pheatmap(mat = mat_,
              annotation_row = annot_, annotation_colors = annoCol,
              gaps_col = gaps_, labels_col = collabs,labels_row = row_labs, gaps_row = row_gaps,
              legend_labels = c('loss','neutral','gain'), legend_breaks = c(min(mat_),1,max(mat_)),
              angle_col = 90,
              color = cols_ , border_color = NA,
              cluster_rows = F, cluster_cols = F, treeheight_row = 0, treeheight_col = 0, clustering_method = 'complete',
              legend = T,
              show_rownames = T, show_colnames = T, silent = T)
p_$gtable$grobs[[2]]$gp = gpar(col = "white", fontsize = 10)      # col labels: p_$gtable
p_$gtable$grobs[[3]]$gp = gpar(col = "white", fontsize = 10)      # row lables
p_$gtable$grobs[[5]]$gp = gpar(col = 'white', fontsize = 10)      # row_annotation_names
p_$gtable$grobs[[6]]$gp = gpar(col = 'white', fontsize = 10)      # annotation_legend
p_$gtable$grobs[[7]]$gp = gpar(col = "white", fontsize = 10)      # legend labels

png(file = 'g.png', width = 10, height = 10, units = 'in', res = 300, bg = 'black')
  grid::grid.newpage()
  grid::grid.draw(p_$gtable)
graphics.off()

