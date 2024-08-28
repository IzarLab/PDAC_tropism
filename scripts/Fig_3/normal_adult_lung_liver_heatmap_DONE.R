# This script generates correlation matrix of liver/lung metastatic PDAC signature genes across healthy human adult liver and lung cells

library(Seurat)
library(SeuratObject)
library(pheatmap)
library(grid)
library(gridExtra)
library(viridis)

# Reading in human adult lung-liver merged data ####

load('Misc/merged_human_adult_lung_liver_processed.RData')

# Reading in PDAC-liver/lung markers ####

pdac_liver_markers = c('HPN','SERPING1','UGT2B7','C1S','CYP4F11','C2','C3',
                       'CFB','NNMT','SLCO1B3','FAH','XDH','SLC16A2','PECR')

pdac_lung_markers  = c('PGC','KCNJ15','ROS1','CXCL17','SLC26A9','CLDN18','ANKRD29','VSIG2',
                       'LAMA3','CEACAM6','ROR1','C4orf19','RTKN2','ATP13A4','CAPN8','TMPRSS2')

# Plotting ####

# heatmap

mat_exp = lng_liv[['RNA']]$data[c(pdac_lung_markers, pdac_liver_markers),]
annot_ = data.frame('PDAC recurrence signature' = rep(c('lung','liver'), c(length(pdac_lung_markers), length(pdac_liver_markers))), stringsAsFactors = F, check.names = F, row.names = c(pdac_lung_markers, pdac_liver_markers))
annot_$`PDAC recurrence signature` = factor(annot_$`PDAC recurrence signature`, levels = c('lung','liver'))
annot_col = list('PDAC recurrence signature' = c(lung = 'blue', liver = 'red'))

corr_ = cor(t(as.matrix(mat_exp)))      # correlation

# heatmap
p_ = pheatmap(mat = corr_,
              annotation_row  = annot_, annotation_col  = annot_, annotation_colors = annot_col,
              color = inferno(n = 100),
              border_color = NA,
              show_colnames = T, show_rownames = T,
              legend_labels = c('Low','Medium','High'), legend_breaks = seq(from = min(corr_), to = max(corr_), length.out = 3),
              cluster_rows = T, cluster_cols = T, clustering_method = 'ward.D2',treeheight_col = 0, treeheight_row = 0, cutree_rows = 2, cutree_cols = 2,
              main = paste0('Correlatin matrix'),
              silent = T)

#  print(p_$gtable): list items
p_$gtable$grobs[[3]]$gp = gpar(fontface = "bold", fontfamily = 'Helvetica')       # column labels
p_$gtable$grobs[[4]]$gp = gpar(fontface = "bold", fontfamily = 'Helvetica')       # row labels
p_$gtable$grobs[[9]]$gp = gpar(fontface = "bold", fontfamily = 'Helvetica')       # annotation legend labels
p_$gtable$grobs[[10]]$gp = gpar(fontface = "bold", fontfamily = 'Helvetica')      # legend labels

p_$gtable$grobs[[6]]$gp = gpar(fontsize = 0, fontfamily = 'Helvetica')      # column annotation
p_$gtable$grobs[[8]]$gp = gpar(fontsize = 0, fontfamily = 'Helvetica')      # row annotation

pdf(file = 'human_liver_lung_atlas_heatmap.pdf', width = 11, height = 9, fonts = 'Helvetica')
  grid::grid.newpage()
  grid::grid.draw(p_$gtable)
graphics.off()

