# This script checks the correlation of liver/lung metastatic PDAC signature genes across liver and lung metastases provided by Moffitt's paper
# Moffitt's cohort is bulk microarray gene expression data of liver and lung metastatic PDACs

library(ggplot2)
library(Seurat)
library(pheatmap)
library(gridExtra)
library(grid)

# Reading in liver/lung met PDAC signatures ####

pdac_liv_prog_mrks = c('HPN','SERPING1','UGT2B7','C1S','CYP4F11','C2','C3',
                       'CFB','NNMT','SLCO1B3','FAH','XDH','SLC16A2','PECR')

pdac_lng_prog_mrks  = c('PGC','KCNJ15','ROS1','CXCL17','SLC26A9','CLDN18','ANKRD29','VSIG2',
                        'LAMA3','CEACAM6','ROR1','C4orf19','RTKN2','ATP13A4','CAPN8','TMPRSS2')

# Reading in Moffitt data ####

load('Misc/Moffitt.RData')

# Adding PDAC-liver/lung-progression outdated gene symbols to the updated ones missing from Moffitt data ####

# gene ID conversion table
geneMeta = read.delim(file = 'Misc/HGNC.tsv', sep = '\t', check.names = T, as.is = T, na.strings = '', quote = '', header = T)[,c('Approved.symbol','Previous.symbols')]
rownames(geneMeta) = geneMeta$Approved.symbol

updated_gene_symbols = pdac_lng_prog_mrks[!pdac_lng_prog_mrks %in% rownames(exprsMat)]
outdated_gene_symbols = geneMeta[updated_gene_symbols, "Previous.symbols"]
names(outdated_gene_symbols) = updated_gene_symbols
outdated_gene_symbols = outdated_gene_symbols[!is.na(outdated_gene_symbols)]
pdac_lng_prog_mrks[pdac_lng_prog_mrks %in% names(outdated_gene_symbols)] = outdated_gene_symbols
pdac_lng_prog_mrks = pdac_lng_prog_mrks[pdac_lng_prog_mrks %in% rownames(exprsMat)]     # only retains those genes measured in met500

updated_gene_symbols = pdac_liv_prog_mrks[!pdac_liv_prog_mrks %in% rownames(exprsMat)]
outdated_gene_symbols = geneMeta[updated_gene_symbols, "Previous.symbols"]
names(outdated_gene_symbols) = updated_gene_symbols
outdated_gene_symbols = outdated_gene_symbols[!is.na(outdated_gene_symbols)]
pdac_liv_prog_mrks[pdac_liv_prog_mrks %in% names(outdated_gene_symbols)] = outdated_gene_symbols
pdac_liv_prog_mrks = pdac_liv_prog_mrks[pdac_liv_prog_mrks %in% rownames(exprsMat)]     # only retains those genes measured in met500

# sub-matrices of liver and lung

lung_smpls = names(sampleMeta)[sampleMeta %in% 'Lung_Metastasis']
liver_smpls = names(sampleMeta)[sampleMeta %in% 'Liver_Metastasis']
exprsMat = exprsMat[c(pdac_lng_prog_mrks, pdac_liv_prog_mrks), c(lung_smpls, liver_smpls)]

# Gene correlation heatmap plotting ####

# Across Moffitt cohort ####

annot_ = data.frame('PDAC metastasis signature' = rep(c('Lung','Liver'), c(length(pdac_lng_prog_mrks), length(pdac_liv_prog_mrks))),
                   stringsAsFactors = T, check.names = F,
                   row.names = c(pdac_lng_prog_mrks, pdac_liv_prog_mrks))
corr_ = cor(t(log1p(exprsMat)))

# plotting

p_ = pheatmap(mat = corr_,
              annotation_row = annot_, annotation_col = annot_, annotation_colors = list('PDAC metastasis signature' = c(Lung = 'blue', Liver = 'red')), annotation_names_row = F, annotation_names_col = F,
              color = viridis::inferno(n = 100),
              border_color = NA,
              cluster_rows = T, cluster_cols = T, clustering_method = 'ward.D', cutree_cols = 2, cutree_rows = 2, treeheight_col = 0, treeheight_row = 0,
              silent = T,
              main = 'Moffitt correlation matrix')
p_$gtable$grobs[[3]]$gp = gpar(fontfamily = 'Helvetica', fontface= 'bold')      # column names
p_$gtable$grobs[[4]]$gp = gpar(fontfamily = 'Helvetica', fontface= 'bold')      # row names
p_$gtable$grobs[[7]]$gp = gpar(fontfamily = 'Helvetica', fontface= 'bold')      # annotation legend keys
p_$gtable$grobs[[8]]$gp = gpar(fontfamily = 'Helvetica', fontface= 'bold')     # legend keys

pdf(file='Moffitt_gene_correlation_matrix_all_samples.pdf', width = 11, height = 8)
  grid::grid.newpage()
  grid::grid.draw(p_$gtable)
graphics.off()

# Per met site ####

pdf(file='Moffitt_gene_corelations_violin_met_samples.pdf', width = 9, height = 7)
for(m_ in list('lung', 'liver'))
{
  if(m_ == 'lung'){ m_ = pdac_lng_prog_mrks; t_ = 'PDAC-lung signature' }else{ m_ = pdac_liv_prog_mrks; t_ = 'PDAC-liver signature'}
  
  cor_lung = cor(t(log1p(exprsMat[m_, lung_smpls])))
  cor_lung = 0.5 * log((1 + cor_lung)/(1 - cor_lung))      # standardizing correlations by Fisher's z transformation
  
  cor_liver = cor(t(log1p(exprsMat[m_, liver_smpls])))
  cor_liver = 0.5 * log((1 + cor_liver)/(1 - cor_liver))
  
  # testing
  
  corrs_lung  = cor_lung [upper.tri(cor_lung, diag = F)]
  corrs_liver = cor_liver[upper.tri(cor_liver, diag = F)]
  set.seed(42)
  p_values = list()
  for(i_ in 1:1000)
  {
    corrs_liver_sub = sample(corrs_liver, size = 50, replace = F)
    corrs_lung_sub  = sample(corrs_lung,  size = 50, replace = F)
    
    p_values[[i_]] = wilcox.test(x = corrs_liver_sub,
                                 y = corrs_lung_sub,
                                 alternative = "two.sided")$p.value
  }
  p_values = format(median(p.adjust(p = p_values, method = 'BH')), scien = T)
  
  # plotting
  
  # violin
  dt_ = data.frame(correlation = c(lung = corrs_lung, liver = corrs_liver), signature = rep( c('in lung smaples','in liver samples'), c(length(corrs_lung),length(corrs_liver))), stringsAsFactors = T)
  q_ =  ggplot(data = dt_, aes(x = signature, y = correlation, fill = signature))+
        theme(text = element_text(face = 'bold', size = 15, family = 'Helvetica'), legend.position = 'none', plot.title = element_text(hjust = 0.5),
              plot.background = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = 'black'), axis.line.y = element_line(color = 'black'),
              panel.grid.minor.y = element_line(color = 'grey60'))+
        labs(y = 'Standardized correlation', title =  paste0(t_,' (p = ',p_values,')'))+
        geom_violin(trim=F)+
        geom_boxplot(width=0.1, size = 1)+
        geom_jitter(shape=16, position=position_jitter(0.2), cex = 1)+
        scale_fill_manual(values = c('red', 'blue'))
        grid::grid.draw(q_)
  
  # heatmap
  
  corr_ = cbind(cor_lung, cor_liver)
  corr_[is.infinite(corr_)] = 1
  corr_[corr_ < 0 ] = 0
  
  p_ = pheatmap( mat = corr_,
                 legend_breaks = c(0,max(corr_)/2,max(corr_)), legend_labels = c('low','mid','high'),
                 color = viridis::inferno(n = 100),
                 border_color = NA,
                 gaps_row = length(m_), gaps_col = length(m_),
                 cluster_rows = F, cluster_cols = F,
                 silent = T,
                 main = paste0(t_,' (p = ',p_values,')'))
  p_$gtable$grobs[[1]]$gp = gpar(fontfamily = 'Helvetica', fontface= 'bold')      # column names
  p_$gtable$grobs[[3]]$gp = gpar(fontfamily = 'Helvetica', fontface= 'bold')      # row names
  p_$gtable$grobs[[4]]$gp = gpar(fontfamily = 'Helvetica', fontface= 'bold')      # annotation legend keys
  p_$gtable$grobs[[5]]$gp = gpar(fontfamily = 'Helvetica', fontface= 'bold')     # legend keys
  
  grid::grid.newpage()
  grid::grid.draw(p_$gtable)
}
graphics.off()

