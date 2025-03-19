# This scripts plots Extended Data Fig. 4a-c; it generates a dot plot of PDAC-recurrence signatures' score per patient.
# rows: PDAC-recurrence signatures
# column: patients
# dot size: fraction of cells significantly expressing that signature
# dot color: signature expression score
# As single-cell suffer from sparsity and imprecise gene expression measurements:
# 1) ALRA wrapped in Seurat package was used to decrease sparsity from 90% to 60%
# 2) an algorithm was used to deal with inconsistent gene expressions
# The output of ALRA:
# The matrix went from 10.00% nonzero to 39.47% nonzero

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

# Reading in PDAC-recurrence signatures ####
message('Reading in data')

# lung
pdac_lung_markers = read.delim('Misc/cell_types_markers_main_subtypes_lung.tsv', header = T, sep = '\t')
pdac_lung_markers = pdac_lung_markers[pdac_lung_markers$cluster %in% 'all' & pdac_lung_markers$avg_log2FC >= 1,]
pdac_lung_markers = trimws(pdac_lung_markers$gene)

# liver
pdac_liver_markers = read.delim('Misc/cell_types_markers_main_subtypes_liver.tsv', header = T, sep = '\t')
pdac_liver_markers = pdac_liver_markers[pdac_liver_markers$cluster %in% 'all' & pdac_liver_markers$avg_log2FC >= 1,]
pdac_liver_markers = trimws(pdac_liver_markers$gene)

# Imputation ####
message('Imputation')

if(!file.exists('compartments/cacner_imputed.RData'))
{
  load(file = 'compartments/cancer.RData')
  s_objs = subset(s_objs, cells = colnames(s_objs)[!s_objs$cell_type %in% c('G2-M','E2F Targets (S)','neural-like progenitor')] )     # all cancer cells except mitotic and NLP
  annot_ = s_objs@meta.data
  
  s_objs = RunALRA(s_objs, assay = 'RNA', k.only = TRUE)      # only computes optimal k WITHOUT performing ALRA; by default it uses normalized data in slot data
  s_objs = RunALRA(s_objs, assay = 'RNA')                   # runs ALRA with the chosen k
  save(s_objs, file = 'compartments/cacner_imputed.RData.RData')
  s_objs = as.matrix(s_objs[['alra']]@data)
}else
{
  load(file = 'compartments/cacner_imputed.RData')
  annot_ = s_objs@meta.data
  s_objs = as.matrix(s_objs[['alra']]@data)
}

# Correcting gene expressions ####
message('Correcting expression values')

pdac_sig_markers = list(lung = pdac_lung_markers, liver = pdac_liver_markers)
for(s_ in names(pdac_sig_markers))
{
  cells_ = rownames(annot_)[annot_$condition %in% s_]
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

# pathways of interest

pathws_ = list(pdac_lung = pdac_lung_markers, pdac_liver = pdac_liver_markers)
names(pathws_) = c('PDAC-lung recurrence signature', 'PDAC-liver recurrence signature')

# for each patient

annot_[,names(pathws_)] = 0     # adding pathways as columns to meta data tables
for(p_ in unique(annot_$orig.ident))
{
  # extracting current patient's cancer cells
  cat('processing ', p_, '\n')
  
  # computing pathway (signature) score
  
  cells_ = rownames(annot_)[annot_$orig.ident %in% p_]      # cells of current patient
  
  for(i_ in 1:length(pathws_))
  {
    pathnm = names(pathws_)[i_]
    cat('\t> ', pathnm, '\n')
    genes_ = pathws_[[i_]]      # genes in the current signature
    
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

# Plotting ####

pdf(file = 'a-c.pdf', width = 9, height = 5)

## dot plot ####

# formatting data frame table for ggplot2
dt_ = list(); i_ = 1
for(p_ in unique(annot_$orig.ident))      # for each patient
{
  for(pathw_ in names(pathws_))     # for each signature (LIV-MR and LUN-MR)
  {
    exprs_ = annot_[annot_$orig.ident %in% p_, c("condition", pathw_)]                                    # signature score per cell
    dt_[[i_]] = data.frame(patient = p_,
                           pathway = pathw_,
                           condition = exprs_$condition[1],
                           exprs = mean(exprs_[,pathw_]),                                                 # signature score per patient
                           fraction = round(length(which(0 < exprs_[,pathw_]))/nrow(exprs_)*100, 2))      # fraction of cells expressing the signature
    i_ = i_ + 1
  }
}
dt_ = do.call(dt_, what = rbind)
dt_$exprs[dt_$exprs == min(dt_$exprs)] = round(min(dt_$exprs),2)      # 0 scores are generated here
dt_$exprs[dt_$exprs == max(dt_$exprs)] = round(max(dt_$exprs),2)
dt_$fraction[dt_$exprs == 0] = 0


dt_$patient = factor(dt_$patient, levels = c('PN9','PN18',                'PN1','PN10','PN23','PN27','PN3','PN8',
                                             'PN4','PN13','PN15','PN17',  'PN14','PN16','PN19','PN2','PN24','PN28','PN5','PN6','PN7'))
p_ =  ggplot(data = dt_, aes(y = pathway, x = patient, size = fraction, color = exprs, shape = condition))+
      theme(panel.background = element_blank(), panel.grid = element_line(color = 'black',linewidth = .05),
            legend.position = 'right', legend.text = element_text(size = 15), legend.key = element_blank(),
            legend.title = element_text(size = 15), legend.spacing.y = unit(.6, 'cm'),
            axis.text.x = element_text(angle = 90),
            text =  element_text(face = 'bold', size = 10, family = 'Helvetica'))+
      labs(y = NULL, x = NULL, color = 'Signature score', size = 'Fraction', shape = 'Initial metastatsis')+
      geom_point()+
      scale_color_viridis(discrete = F, option = 'H', breaks = c(min(dt_$exprs), round(diff(range(dt_$exprs))/2,2), max(dt_$exprs)) )+
      guides(color = guide_colorbar(ticks.colour = NA, label.position = 'right'),
             size = guide_legend(override.aes = list(size = 1:5)),
             shape = guide_legend(override.aes = list(size = 4:5)))+
      coord_fixed(ratio = 3)

plot(p_)

## density plots ####

# formatting data frame table for ggplot2
liv_1_lung_2 = annot_[annot_$condition %in% 'liver',]                           # initial liver recurrence, secondary lung recurrence
liv_1_lung_2$initial_score = liv_1_lung_2$`PDAC-liver recurrence signature`     # initial recurrence score
liv_1_lung_2$seond_score = liv_1_lung_2$`PDAC-lung recurrence signature`        # secondary recurrence score
liv_1_lung_2 = liv_1_lung_2[0 < liv_1_lung_2$seond_score,]                      # only non-zero secondary cell scores; this is because many of these cells
                                                                                # are zero for secondary signature that will predominate the density curve

lng_1_liv_2  = annot_[annot_$condition %in% 'lung', ]                           # initial lung recurrence, secondary liver recurrence
lng_1_liv_2$initial_score = lng_1_liv_2$`PDAC-lung recurrence signature`
lng_1_liv_2$seond_score = lng_1_liv_2$`PDAC-liver recurrence signature`
lng_1_liv_2 = lng_1_liv_2[0 < lng_1_liv_2$seond_score,]

cols_ = c("#FF0000",  # red
          "#A52A2A",  # brown
          "#0000FF",  # blue
          "#FFA500",  # orange
          "#800080",  # purple
          "#FFFF00",  # yellow
          "#FFC0CB",  # pink
          "#FC8D62",  # Salmon
          "#00FFFF",  # cyan
          "#00FF00",  # green
          "#FF00FF",  # magenta
          "#006400",  # darkgreen
          "#00008B")  # darkblue

j_ = 1
for(dt_ in list(liv_1_lung_2, lng_1_liv_2))     # for each signature
{
  den_1 = density(dt_$initial_score)      # density of non-zero initial recurrence scores
  max_x1 = den_1$x[which.max(den_1$y)]
  
  den_2 = density(dt_$seond_score)        # density of non-zero secondary recurrence scores
  max_x2 = den_2$x[which.max(den_2$y)]
  
  # for each cell
  for(i_ in 1:nrow(dt_))
  {
    # initial recurrence
    initial_score = dt_$initial_score[i_]
    ub_ = which(initial_score <= den_1$x)[1]      # closest density bin to the score[i_]
    dt_[i_,'den_x1'] = den_1$x[ub_]
    dt_[i_,'den_y1'] = den_1$y[ub_]
    
    # secondary recurrence
    seond_score = dt_$seond_score[i_]
    ub_ = which(seond_score <= den_2$x)[1]        # closest density bin to the score[i_]
    dt_[i_,'den_x2'] = den_2$x[ub_]
    dt_[i_,'den_y2'] = den_2$y[ub_]
  }
  dt_ = dt_[, c("orig.ident",
                'den_x1', 'den_y1',
                'den_x2', 'den_y2')]
  dt_$secondary = F
  dt_$secondary[dt_$orig.ident %in% c('PN9','PN18','PN4','PN13','PN15','PN17')] = T
  
  # density plot
  t_ = NULL
  if(j_ == 1){t_ = 'LUN-MR signatrue in LIV-MR patients'; s_ = 0.5}else{t_ = 'LIV-MR signatrue in LUN-MR patients'; s_ = 0.1}; j_ = j_ + 1
  p_ =  ggplot(dt_, aes(x = den_x2, y = den_y2)) +
        theme(panel.background = element_blank(),plot.title = element_text(hjust = 0.5),
              axis.line = element_line(color = 'black'),
              text = element_text(size = 20, face = 'bold', family = 'Helvetica'))+
        labs(title = t_,x = "Gene Signature Score (x)", y = "Density (y)", color = 'Patient')+
        
        # track 1: Curve for x and y
        geom_line(color = "black", linewidth = 1.2)+
        geom_line(aes(x = den_x1, y = den_y1), color = "grey70", linewidth = 1.2, linetype = 'dashed')+
        scale_x_continuous(breaks =        c(max_x2, mean(c(max_x2,max_x1)), max_x1, max(dt_$den_x1, dt_$den_x2)),
                           labels = round(c(max_x2, mean(c(max_x2,max_x1)), max_x1, max(dt_$den_x1, dt_$den_x2)),2) )+
        
        # track 2: ticks for orig.ident
        geom_rug(aes(x = den_x2, color = orig.ident), sides = "t", linewidth = s_, position = "jitter", show.legend = T)+
        
        # track 3: ticks for secondary metastasis, slightly below the x-axis
        geom_rug(data = dt_[dt_$secondary, ], aes(x = den_x2, color = secondary), sides = "b", position = "jitter", linewidth = s_)+
        scale_color_manual(values = c(cols_[1:length(unique(dt_$orig.ident))],"#8B0000"))+
        guides(color = guide_legend(override.aes = list(linewidth = 2)))

  plot(p_)
}

graphics.off()
