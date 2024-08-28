# This script generates a heatmap of canonical markers of macrophage subtypes

library(Seurat)
library(SeuratObject)
options(Seurat.object.assay.version = "v4")
library(pheatmap)
library(grid)
library(viridis)

load(file = 'compartments/myeloid.RData')
s_objs = subset(s_objs, subset = cell_type %in% c('M1-like','M2a-like','M2c-like','M2d-like'))
DefaultAssay(s_objs) = 'RNA'

dt_ = matrix(data = 0, ncol = 21, nrow = length(unique(s_objs$cell_type)), dimnames = list(types = c('M1-like',
                                                                                                     'M2a-like',
                                                                                                     'M2c-like',
                                                                                                     'M2d-like'),
                                                                                           Gene = c('CD80','CD86','FCGR1A','FCGR3A',
                                                                                                    'CD163','MRC1','CD200R1','NRP1',
                                                                                                    'SCARB1','TLR8','CCL18',
                                                                                                    'VEGFA','ITGB3','CCL18_rep','FLT1','SEMA3C','CXCL8','CCL2','ADM','CSF1','MARCO')) )
for(r_ in rownames(dt_))
{
  s_ = subset(s_objs, subset = cell_type == r_)
  for(g_ in colnames(dt_))
  {
    if(g_ == 'CCL18_rep'){ dt_[r_,g_] = mean(colMeans(s_[["RNA"]]@data['CCL18',,drop = F])); next() }
    dt_[r_,g_] = mean(colMeans(s_[["RNA"]]@data[g_,,drop = F]))
  }
}

dt_ = apply(X = dt_, MARGIN = 2, FUN = function(c_){ (c_ - min(c_))/diff(range(c_)) })


p_ =
  pheatmap( mat = dt_,
            cluster_cols = F, cluster_rows = F, gaps_col = c(4,8,11),
            color = viridis(n = 100, option = 'B', end = .85), angle_col = 315, breaks = seq(from = 0, to = 1, length.out = 101),
            legend_labels = c('0','1'), legend_breaks = c(0,1),
            fontsize = 25, cellheight = 35, cellwidth = 40, silent = T)

p_$gtable$grobs[[2]]$gp = gpar(fontface = "bold", fontfamily = 'Helvetica')     # column labels
p_$gtable$grobs[[3]]$gp = gpar(fontface = "bold", fontfamily = 'Helvetica')     # row labels
p_$gtable$grobs[[4]]$gp = gpar(fontface = "bold", fontfamily = 'Helvetica')     # legend
p_$gtable$grobs[[4]]$children[[1]]$hjust = -4
p_$gtable$grobs[[4]]$children[[2]]$hjust = -3

pdf(file = 'macropahges_canonical_markers_heatmap.pdf', width = 15, height = 5.5)
grid::grid.newpage()
grid::grid.draw(p_$gtable)
graphics.off()
