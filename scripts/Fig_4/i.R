# This script plots Fig. 4i; it performs Wilcox rank-sum test between LIV-MR and LUN-MR genes to find genes
#  with significantly different gene dosage effect (GDX) produced by Echidna tool (https://github.com/azizilab/echidna).

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratWrappers)
library(rstatix)
library(ggplot2)
library(ggrepel)
library(viridis)

# Reading in liver/lung met PDAC-recurrence signatures ####

# reading gene ordering file
genes_annot = read.delim(file = 'Misc/gene_ordering_file.txt', header = F, sep = '\t', as.is = T, check.names = F)
colnames(genes_annot) = c('gene', 'chr', 'str', 'end')
rownames(genes_annot) = genes_annot$gene

# PDAC-liver metastasis (LIV-MR) signature
pdac_liver_markers = read.delim('Misc/cell_types_markers_main_subtypes_liver.tsv', header = T, sep = '\t')
pdac_liver_markers = pdac_liver_markers[pdac_liver_markers$cluster %in% 'all',]
pdac_liver_markers$avg_log2FC = -pdac_liver_markers$avg_log2FC      # the original DGEA liver was taken as reference so negative L2FC means genes upregulated in liver singnature
pdac_liver_markers$gene = trimws(pdac_liver_markers$gene)           # because of Excel issue, genes were padded with one white-space character
pdac_liver_markers = pdac_liver_markers[pdac_liver_markers$gene %in% genes_annot$gene,]
pdac_liver_markers = pdac_liver_markers[1:500,]                     # top 500 genes; inferCNV for low-quality genes is quite noisy; top 500 high-quality DE genes; try 1000, or all the genes

# PDAC-lung metastasis (LIV-MR) signature
pdac_lung_markers = read.delim('Misc/cell_types_markers_main_subtypes_lung.tsv', header = T, sep = '\t')
pdac_lung_markers = pdac_lung_markers[pdac_lung_markers$cluster %in% 'all' ,]
pdac_lung_markers$gene = trimws(pdac_lung_markers$gene)
pdac_lung_markers = pdac_lung_markers[pdac_lung_markers$gene %in% genes_annot$gene,]
pdac_lung_markers = pdac_lung_markers[1:500,]

mrks_ = c(pdac_liver_markers$gene, pdac_lung_markers$gene)      # combined genes in LIV-MR and LUN-MR signatures

# Extracting gene dosage effect (GDX) and gene dosage (CNV) for LIV-MR and LUN-MR genes ####

liv_samples = c("PN1","PN10","PN18","PN23","PN27","PN8","PN9",'PN3')
lng_samples = c("PN13","PN14","PN15","PN16","PN17","PN19","PN2","PN24","PN28","PN4","PN5","PN6","PN7")

gdx_liv = gdx_lng = cnv_liv = cnv_lng = vector(mode = 'list', length = length(mrks_))
names(gdx_liv) = names(gdx_lng) = names(cnv_liv) = names(cnv_lng) = mrks_

for(fs_ in list.files(path = 'Misc/gdx/', pattern = 'delta_eta'))
{
  message(fs_)

  fs_ = strsplit(fs_, split = '_')[[1]][1]

  gdx_ = read.csv(file = paste0('Misc/gdx/',fs_,'_gde.csv'), header = T, as.is = T, check.names = F, row.names = 1)           # gene dosage effect
  cnv_ = read.csv(file = paste0('Misc/gdx/',fs_,'_delta_eta.csv'), header = T, as.is = T, check.names = F, row.names = 1)     # inferred CNV

  common_genes = intersect(intersect(mrks_, rownames(gdx_)),rownames(cnv_))     # cnv and gdx are matched
  gdx_ = gdx_[common_genes,]
  cnv_ = cnv_[common_genes,]

  if(fs_ %in% liv_samples)
  {
    for(i_ in 1:nrow(gdx_))
    {
      dt_ = gdx_[i_,]
      gene_ = rownames(dt_)
      if(!gene_ %in% mrks_){ next() }
      gdx_liv[[gene_]] = append(gdx_liv[[gene_]], as.numeric(dt_))
      cnv_liv[[gene_]] = append(cnv_liv[[gene_]], as.numeric(dt_))
    }
  }else
  {
    for(i_ in 1:nrow(gdx_))
    {
      dt_ = gdx_[i_,]
      gene_ = rownames(dt_)
      if(!gene_ %in% mrks_){ next() }
      gdx_lng[[gene_]] = append(gdx_lng[[gene_]], as.numeric(dt_))
      cnv_lng[[gene_]] = append(cnv_lng[[gene_]], as.numeric(dt_))
    }
  }
}

## separating GDX and CNV for LIV-MR and LUN-MR ####

# loading cancer data to extract gene expressions
load('compartments/cancer_subtypes.RData')
cells_ = s_objs$condition
cancer_ = as.matrix(s_objs[['RNA']]@data); rm(s_objs)

# for ALL genes in LIV-MR and LUN-Mr signatures
  
dt_liv = dt_lng = list()      # a table to store GDX and gene dosage (CNV) for LIV-MR or LUN-MR genes
for(i_ in 1:length(mrks_))
{
  if(!is.null(gdx_liv[[i_]]) & !is.null(gdx_lng[[i_]]))     # making sure a gene has been measured both in lIV-MR and LUN-MR cells
                                                            # this is necessary for the CNV heatmap underneath the lollipops plot
  {
    # LIV-MR cells
    dt_liv_t = data.frame(gene               = mrks_[i_],
                          gene_dosage_effect = abs(mean(gdx_liv[[i_]])),      # gene names of gdx_liv/lng and cnv_liv/lng match mrks_,
                          gene_dosage        = mean(cnv_liv[[i_]]),           # CNV was matched to gdx so it cannot be NA,
                          expression         = mean(cancer_[mrks_[i_], cells_ %in% 'liver']))
    # LUN-MR cells
    dt_lng_t = data.frame(gene               = mrks_[i_],
                          gene_dosage_effect = abs(mean(gdx_lng[[i_]])),
                          gene_dosage        = mean(cnv_lng[[i_]]),
                          expression         = mean(cancer_[mrks_[i_], cells_ %in% 'lung']))
    
    dt_liv_t$signature = dt_lng_t$signature = if(mrks_[i_] %in% pdac_liver_markers$gene){ 'LIV-MR' }else{ 'LUN-MR' }
      
    # Wilcox test on GDX
    dt_liv_t$p_val = dt_lng_t$p_val = wilcox.test(x = gdx_liv[[i_]], y = gdx_lng[[i_]], alternative = 'two.sided')$p.value
    
    dt_liv[[mrks_[i_]]] = dt_liv_t
    dt_lng[[mrks_[i_]]] = dt_lng_t
  }
}

# dt_liv and dt_lng have the same genes in the same order as well as the same p_values for wilcox test on GDX
dt_liv = do.call(dt_liv, what = rbind)
dt_lng = do.call(dt_lng, what = rbind)
dt_liv$p_val = dt_lng$p_val = p.adjust(p = dt_liv$p_val, method = 'BH')

# Plotting ####

# for LIV-MR and LUN-MR

p_ = list()     # list of plots
for(k_ in 1:2)
{
  dt_ = if(k_ == 1){ dt_liv }else{ dt_lng }
  
  ## ordering genes in dt_ (LIV-MR/LUN-MR table) according to their genomic loci ####
  
  genes_annot_tmp = genes_annot
  ordered_genes = rownames(genes_annot_tmp)[rownames(genes_annot_tmp) %in% rownames(dt_)]
  genes_annot_tmp = genes_annot_tmp[ordered_genes,]
  dt_ = dt_[ordered_genes,]
  
  ## shifting (offsetting) starts of each gene's locus by 10e5 according to its chromosome number ####
  # this step is necessary for plotting density of genomic loci so that genes at chr12 follow those within chr11 on x-axis

  chrs_ = unique(genes_annot_tmp$chr)     # chromosomes are already sorted; these are chromosomes covered in dt_

  # shifting
  for(i_ in 2:length(chrs_))
  {
    ind_c = which(genes_annot_tmp$chr %in% chrs_[i_])
    ind_p = ind_c[1]-1
    str_ =  genes_annot_tmp$str[ind_p] + 10e5                                                         # the end of previous chromosome + 10e5
    genes_annot_tmp$str[ind_c] = genes_annot_tmp$str[ind_c]-(genes_annot_tmp$str[ind_c[1]]-str_)      # shifting all str loci of current chromosome
  }
  dt_$str = round(genes_annot_tmp$str/1e5)      # 10e5/1e5 = 10 this guarantees that first str of each chromosome is 10 points greater than last chromosome's str
  dt_$str = dt_$str - min(dt_$str)+1
  dt_$chr = genes_annot_tmp$chr                 # just for control

  # making sure all genes have nonoverlapping unique loci after rounding
  while(T)
  {
    ind_ = which(diff(dt_$str) %in% 0) + 1
    if(length(ind_) == 0){ break }
    dt_$str[ind_] = dt_$str[ind_] + 1
  }

  ## preparing for plotting: labels, coloration, normalization, etc ####
  
  # smoothing CNV by Guassian kernel
  dt_n = data.frame(gene = NA,
                    gene_dosage_effect = NA,
                    gene_dosage = rep(0,max(dt_$str)),
                    expression = NA,
                    signature = NA,
                    p_val = NA,
                    str = 1:max(dt_$str),
                    chr = NA)
  dt_n[dt_$str,] = dt_
  dt_n$gene_dosage_smooth = ksmooth(y = dt_n$gene_dosage, x = 1:nrow(dt_n), bandwidth = 500)$y
  dt_ = dt_n; rm(dt_n); gc()
  
  # normalizing gene expression between [0,1]
  dt_$expression = (dt_$expression - min(dt_$expression, na.rm = T))/diff(range(dt_$expression, na.rm = T))
  
  # determining the number of significant digits for GDX and gene dosage (CNV)
  dt_[dt_$gene_dosage_effect %in% min(dt_$gene_dosage_effect,na.rm = T),"gene_dosage_effect"] = round(min(dt_$gene_dosage_effect,na.rm = T),2)
  dt_[dt_$gene_dosage_effect %in% max(dt_$gene_dosage_effect,na.rm = T),"gene_dosage_effect"] = round(max(dt_$gene_dosage_effect,na.rm = T),2)
  dt_[dt_$gene_dosage_smooth %in% min(dt_$gene_dosage_smooth,na.rm = T),"gene_dosage_smooth"] = round(min(dt_$gene_dosage_smooth,na.rm = T),2)
  dt_[dt_$gene_dosage_smooth %in% max(dt_$gene_dosage_smooth,na.rm = T),"gene_dosage_smooth"] = round(max(dt_$gene_dosage_smooth,na.rm = T),2)
  
  # labels
  labs_x = sapply(X = na.omit(unique(dt_$chr)), FUN = function(c_){ max(dt_$str[dt_$chr %in% c_])})     # chromosome labels
  labs_pnt = rep(NA, nrow(dt_))
  labs_pnt[which(dt_$p_val <= 0.01)] = dt_$gene[which(dt_$p_val <= 0.01)]
  
  # subseeting table
  col_ = t_ = NA
  if(k_ == 1)
  {
    dt_[dt_$signature %in% 'LUN-MR', c("expression","gene_dosage_effect")] = NA     # not plotting genes from LUN-MR signature
    col_ = 'red'            # coloration for lollipops' sticks
    t_ = 'LIV-MR cells'     # subtitle
  }else
  {
    dt_[dt_$signature %in% 'LIV-MR', c("expression","gene_dosage_effect")] = NA     # not plotting genes from LIV-MR signature
    col_ = 'blue'
    t_ = 'LUN-MR cells'
  }
  
  ## plotting ####
  
  p_[[k_]] =  ggplot(data = dt_, aes(x = str))+
              theme(panel.background = element_blank(), panel.grid = element_blank(),
                    axis.line = element_line(color = 'black'),
                    axis.text.x = element_text(angle = 90),
                    text = element_text(face = 'bold', size = 20), plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
                    legend.title = element_text(vjust = 2.5))+
              labs(title = "Gene dosage effect (GDX)", subtitle = t_, x = "Gene", y = "Expression", size = 'GDX', color = 'CNV', fill= 'Signature')+
              # coord_trans(y = 'log1p')+                                                                                                                             # does not apply log1p on statics before plotting
              
              # track 1: density of genomic loci
              geom_density(data = dt_[!is.na(dt_$gene_dosage_effect),], aes(y = after_stat(scaled)), adjust = 0.1, fill = "grey", alpha = 0.35, color = NA)+

              # track 2: lollipops of gene expression and GDX
              geom_segment(aes(x = str, xend = str, y = 0, yend = expression), color = col_, linewidth = 0.1)+                                                        # sticks

              geom_point(data = dt_[!is.na(dt_$gene_dosage_effect),], aes(y = expression, size = gene_dosage_effect, fill = signature), shape = 21, alpha = 0.7)+     # points
              scale_fill_manual(values = col_)+                                                                                                                       # manual color scale for signature
              scale_size_continuous(breaks = round(c(min(dt_$gene_dosage_effect,na.rm = T),mean(dt_$gene_dosage_effect,na.rm = T),max(dt_$gene_dosage_effect,na.rm = T)),2))+
              guides(fill = guide_legend(override.aes = list(size = 6, alpha = 1)))+
              
              geom_label_repel(aes(y = expression, label = labs_pnt), box.padding = .75, max.overlaps = 20,                                                           # significant genes
                               fill = 'yellow', color = 'black', fontface = 'bold', segment.color = "black", alpha = 0.7)+
              
              # track 3: CNV heatmap
              geom_rect(xmin = min(dt_$str), xmax = max(dt_$str), ymin = -0.5, ymax = -0.0115, fill = "grey96", alpha = 0.3, inherit.aes = FALSE)+                    # background of CNV heatmap
              geom_rug(aes(y = .5, color = gene_dosage_smooth), sides = "b", position = "identity", linewidth = .1, alpha = 1)+                                       # bars for CNV
              scale_color_gradient2(low = "blue", mid = "grey95", high = "red", midpoint = 0, breaks = round(c(min(dt_$gene_dosage_smooth),0,max(dt_$gene_dosage_smooth)),2))+
              scale_x_continuous(breaks = labs_x, labels = names(labs_x), expand = c(0.01,0))                                                                         # chromosome labels
}

## saving plot ####

p_ = p_[[1]]/p_[[2]]
jpeg(filename = 'i.jpg', width = 30, height = 16, units = 'in', res = 200)
plot(p_)
graphics.off()
# ggsave(plot = p_, filename = 'GDX_exclusive.pdf', device = 'pdf', width = 30, height = 16)

write.table(dt_liv[dt_liv$p_val <= 0.01,c("gene","signature","gene_dosage_effect")], file = 'GDX_gene.tsv', quote = F, sep = '\t', row.names = F, col.names = T)      # pvalues and genes are the smae in dt_liv and dt_lng

