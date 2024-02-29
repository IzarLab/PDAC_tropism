# This scripts generates a heatmap of canonical markers of CAF subpopulations

library(Seurat)
library(SeuratObject)
options(Seurat.object.assay.version = "v4")
library(pheatmap)
library(grid)
library(viridis)

load(file = 'compartments/CAF.RData')
s_objs = subset(s_objs, subset = cell_type %in% c('CXCL14 fibroblast','proaxonogenic fibroblast','inflammatory fibroblast','myofibroblast','activated stellate'))
DefaultAssay(s_objs) = 'RNA'

dt_ = matrix(data = 0, ncol = 23, nrow = length(unique(s_objs$cell_type)), dimnames = list(types = c('activated stellate',
                                                                                                     'proaxonogenic fibroblast',
                                                                                                     'myofibroblast',
                                                                                                     'inflammatory fibroblast',
                                                                                                     'CXCL14 fibroblast'),
                                                                                            Gene = c('COL1A1','SFRP2','LUM','BGN','MYH11',
                                                                                                     'NRP1','NRXN3','SEMA6D','NFASC','CXCL12',
                                                                                                     'POSTN','JPH2','SPEG','MYH11_rep','ACKR3',
                                                                                                     'CXCL2','ICAM1','NFATC1','NFATC2',
                                                                                                     'CXCL14','MEOX1','ANO2','SEZ6L')) )
for(r_ in rownames(dt_))
{
  s_ = subset(s_objs, subset = cell_type == r_)
  for(g_ in colnames(dt_))
  {
    if(g_ == 'MYH11_rep')
    {
      dt_[r_,'MYH11_rep'] = mean(colMeans(s_[["RNA"]]@data['MYH11',,drop = F]))
    }else
    {
      dt_[r_,g_] = mean(colMeans(s_[["RNA"]]@data[g_,,drop = F]))
    }
  }
}

dt_ = apply(X = dt_, MARGIN = 2, FUN = function(c_){ (c_ - min(c_))/diff(range(c_)) })


p_ =
  pheatmap( mat = dt_,
               cluster_cols = F, cluster_rows = F, gaps_col = c(5,10,15,19),
               color = viridis(n = 100, option = 'B', end = .85), angle_col = 315, breaks = seq(from = 0, to = 1, length.out = 101),
               legend_labels = c('0','1'), legend_breaks = c(0,1),
               fontsize = 25, cellheight = 35, cellwidth = 40, silent = T)

p_$gtable$grobs[[2]]$gp = gpar(fontface = "bold", fontfamily = 'Helvetica')     # column labels
p_$gtable$grobs[[3]]$gp = gpar(fontface = "bold", fontfamily = 'Helvetica')     # row labels
p_$gtable$grobs[[4]]$gp = gpar(fontface = "bold", fontfamily = 'Helvetica')     # legend
p_$gtable$grobs[[4]]$children[[1]]$hjust = -4
p_$gtable$grobs[[4]]$children[[2]]$hjust = -3

pdf(file = 'CAF_canonical_markers_heatmap.pdf', width = 20, height = 5.5)
grid::grid.newpage()
grid::grid.draw(p_$gtable)
graphics.off()
