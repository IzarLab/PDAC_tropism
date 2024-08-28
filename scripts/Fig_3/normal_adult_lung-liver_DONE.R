# This script takes healthy adult lung and liver data and checks liver/lung metastatic PDAC signatures in them

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggridges)
library(patchwork)
library(cowplot)
set.seed(42)


load(file = 'Misc/merged_human_adult_lung_liver_processed.RData')

pdf(file = 'Human_Adult_Lung_Liver_Atlas_lung.pdf', width = 15, height = 10)

DimPlot(lng_liv, pt.size = 1, cells = colnames(lng_liv), raster = F, group.by = 'cell_type', label = T, label.size = 7, repel = T, label.box = T)+
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())+
  labs(title = NULL) & NoLegend()

# liver/lung metastatic PDAC signatures

pdac_liver_markers = c('HPN','SERPING1','UGT2B7','C1S','CYP4F11','C2','C3',
                       'CFB','NNMT','SLCO1B3','FAH','XDH','SLC16A2','PECR')

pdac_lung_markers  = c('PGC','KCNJ15','ROS1','CXCL17','SLC26A9','CLDN18','ANKRD29','VSIG2',
                       'LAMA3','CEACAM6','ROR1','C4orf19','RTKN2','ATP13A4','CAPN8','TMPRSS2')

# Visualization ####

# Density plots

for(c_ in list(c('Alveolar Epithelial Type 1','Alveolar Epithelial Type 2'), 'Hepatocytes'))
{
  s_ = subset(x = lng_liv, subset = cell_type %in% c_)
  
  liver_mat = s_[['RNA']]$data[pdac_liver_markers,,drop = F]
  lung_mat = s_[['RNA']]$data[pdac_lung_markers,,drop = F]
  liv_den = density(colMeans(liver_mat), adjust = 3)
  lng_den = density(colMeans(lung_mat), adjust = 3)
  dt_ = data.frame(exprs = c(liv_den$x, lng_den$x), density = c(liv_den$y, lng_den$y), condition = rep(x = c('Liver', 'Lung'), times = c(length(liv_den$x), length(lng_den$x))), stringsAsFactors = T)
  
  t_ = if('Hepatocytes' %in% c_){'Adult human liver parenchymal cells'}else{'Adult human lung parenchymal cells'}
  p_ =  ggplot(data = dt_, aes(x = exprs, y = density, fill = condition, color = condition))+
    theme(plot.background = element_blank(), panel.background = element_rect(fill = NA), axis.line.y = element_line(color = 'black'), axis.line.x = element_line(color = 'black'),
          legend.position = 'bottom', text = element_text(face = 'bold',size = 20), legend.key=element_rect(fill="white"), plot.title = element_text(hjust = 0.5))+
    labs(x = 'Signature expression', y = 'Density', color = 'PDAC recurrence signature', title = t_)+
    geom_line(linewidth = 1.5)+
    geom_area(position = "identity", alpha = .8, show.legend = F)+
    geom_vline(xintercept = c(liv_den$x[which.max(liv_den$y)], lng_den$x[which.max(lng_den$y)]), linetype = 'dashed' )+
    scale_color_manual(values = c('red','blue'))+
    scale_fill_manual(values = c('red','blue'))+
    scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0,0.05)))+
    facet_grid(condition ~ ., scales='free')
  plot(p_)
}

# Violin plot

lng_liv@meta.data$avg_liv = colMeans(lng_liv[['RNA']]$data[pdac_liver_markers,])
lng_liv@meta.data$avg_lng = colMeans(lng_liv[['RNA']]$data[pdac_lung_markers,])

lng_liv$cell_type = factor(x = lng_liv$cell_type, levels = c('Hepatocytes','Cholangiocytes',
                                                   'Alveolar Epithelial Type 1','Alveolar Epithelial Type 2','Club','Mucous','Ciliated','Basal'))
p_ =  VlnPlot(lng_liv, features = 'avg_liv', group.by = "cell_type", pt.size = 0, combine = FALSE, log = T, sort = F)[[1]]+
      theme(legend.position = 'none', axis.text.x = element_text(angle = 90), text = element_text(face = 'bold'))+
      labs(title = 'PDAC-liver signature', y = 'Signature expression', x = NULL)+
      stat_summary(fun = median, geom ='point', size = 19, colour = "black", shape = 95)
q_ =  VlnPlot(lng_liv, features = 'avg_lng', group.by = "cell_type", pt.size = 0, combine = FALSE, log = T, sort = F)[[1]]+
      theme(legend.position = 'none',  axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90),
            text = element_text(face = 'bold'))+
      labs(title = 'PDAC-lung signature', y = 'Signature expression', x = NULL)+
      stat_summary(fun = median, geom='point', size = 19, colour = "black", shape = 95)

r_ = range(range(p_$data$avg_liv), range(q_$data$avg_lng))
f_ = p_ + scale_y_continuous(limits = r_) + q_ + scale_y_continuous( limits = r_ )
plot(f_)

graphics.off()

# Testing ####

set.seed(42)
hepatocytes_ = names(lng_liv$cell_type[lng_liv$cell_type %in% 'Hepatocytes'])
alveoli_ = names(lng_liv$cell_type[lng_liv$cell_type %in% c('Alveolar Epithelial Type 1','Alveolar Epithelial Type 2')])
for(c_ in list(hepatocytes_, alveoli_))
{
  p_values = list()
  
  for(i_ in 1:1000)
  {
    c_sub = sample(c_, size = 50, replace = F)
    p_values[[i_]] = wilcox.test(x = lng_liv$avg_liv[c_sub],
                                 y = lng_liv$avg_lng[c_sub],
                                 alternative = "two.sided")$p.value
  }
  message(format(median(p.adjust(p = p_values, method = 'BH')), scien = T))
}
