# This script plots Fig. 2b and e (CNV heatmap and UMAP for P10). P10 refers to patient PN14 in the Supplementary Table 1.

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(infercnv)
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(viridis)

# STEP 1: Reading in data ####

genes_chr = read.delim(file = 'Misc/gene_ordering_file.txt',
                       header = F, sep = '\t', quote = "", as.is = T, check.names = F, col.names = c('gene','chr','str','end'))     # gene-chromosome file; MUST BE THE ONE USED BY inferCNV
rownames(genes_chr) = genes_chr$gene

load(file = 'Misc/PN14.RData')                                                                                                      # snRNAseq data, reference and observation matrices returned by inferCNV

ref_mat = ref_mat[grepl(x = rownames(ref_mat), pattern = 'PN14_'),]                                                                 # refernce cells in this sample
rownames(ref_mat) = gsub(pattern = 'PN14_', replacement = '', x = rownames(ref_mat))                                                # removing PNXX_ from the beginning of reference cell names

not_processed_cells = which(!colnames(s_obj) %in% c(rownames(ref_mat), rownames(obs_mat)))                                          # only a subset of T cells from each patient comprises reference cells, the rest is absent in ref_mat too 
abs_mat = matrix(data = 1, nrow = length(not_processed_cells), ncol = ncol(obs_mat), dimnames = list(colnames(s_obj)[not_processed_cells], colnames(obs_mat)))                                                                                                                                                                                                              

cnv_mat = do.call(list(ref_mat, obs_mat, abs_mat), what = rbind); rm(ref_mat, obs_mat, abs_mat); gc()
cnv_mat = cnv_mat[colnames(s_obj),]                                                                                                                                         # reordering cells in cnv_mat to match original count matrix
genes_ = colnames(cnv_mat)                                                                                                                                                  # all genes in this samples
chrs_ = genes_chr[genes_,"chr"]                                                                                                                                             # chromosome of each gene
names(genes_) = chrs_

# STEP 2: CNV heatmap ####

clusters_ = s_obj$seurat_clusters                                                       # Att.: clusters to exclude from CNV mat: s_obj$seurat_clusters[!s_obj$seurat_clusters %in% c(15:24)]
clusters_ = droplevels(clusters_)
clusters_= sort(clusters_)                                                              # sorting clusters in ascending order

mat_ = cnv_mat[names(clusters_),]                                                       # rearranging CNV matrix by cluster order
breaks_ = seq(min(mat_), max(mat_), length.out = 16)                                    # break points for heatmap; in inferCNV there are 15 intervals/color shades. In pheatmap breaks should be one element larger than color vector

gaps_cols = sapply(X = unique(chrs_), FUN = function(c_){ max(which(chrs_ == c_))})     # vertical gaps in heatmap, one for each chromosome
plot_list = list()                                                                      # list to store plots
for(clst_ in levels(clusters_))
{
  # main heatmap
  message(clst_)
  tmp_ = mat_[clusters_ %in% clst_,]
  p_ = pheatmap(mat = tmp_,
                color = colorRampPalette(colors =  c('blue4','white','red4'))(15),
                border_color = NA,
                gaps_col = gaps_cols,
                cluster_cols = F, cluster_rows = T, clustering_method = 'ward.D', treeheight_row = 0,
                show_rownames = F, show_colnames = F,
                cellheight = 32/nrow(tmp_),
                breaks = breaks_,
                legend = F,
                silent = T)
  plot_list[[clst_]] = p_[[4]]

  # right heatmap
  labs_row = rep('',nrow(tmp_))
  labs_row[round(length(labs_row)/2)] = clst_
  p_ = pheatmap(mat = matrix(0, ncol = 1, nrow = nrow(tmp_)),
                color = 'white',
                border_color = NA,
                cluster_rows = F, cluster_cols = F,
                show_rownames = T,labels_row = labs_row,
                show_colnames = F,
                fontsize_row = 10,
                cellheight = 32/nrow(tmp_), cellwidth = .75,
                breaks = 0,
                legend = F,
                silent = T)
  p_[[4]]$grobs[[2]]$gp$col = 'white'
  plot_list[[length(plot_list)+1]] = p_[[4]]
}

# adding chromosome names to the bottom of heatmap
labels_col = chrs_
labels_col[-gaps_cols] = ''                                                             # column names should contain only one chromosome name for each column
p_ = pheatmap(mat = matrix(0, nrow = 1, ncol = ncol(mat_)),
              color = 'white',
              border_color = NA,
              gaps_col = gaps_cols,
              cluster_rows = F, cluster_cols = F,
              show_rownames = F, show_colnames = T, labels_col = labels_col,
              cellheight = 0.75,
              angle_col = 90, fontsize_col = 10,
              breaks = 0,
              legend = F,
              silent = T)
p_[[4]]$grobs[[2]]$gp$col = 'white'
plot_list[[length(plot_list)+1]] = p_[[4]]

# plotting
h_ = .6*length(levels(clusters_))
jpeg(file = paste0('b.png'), width = 1.7*h_, height = h_, units = 'in', res = 600, bg = 'black')
  grid.arrange(grobs= plot_list, nrow = length(levels(clusters_))+1, ncol = 2, as.table = T, widths = c(1, 0.03))
graphics.off()


# STEP 3: Calculating average CNV per chromosome ####

neutral_range = levels(cut(x = as.matrix(cnv_mat), breaks = 15, include.lowest = T, right = T))[8]      # the 8th interval in inferCNV is considered wildtype CNV
neutral_range = c(lower = as.numeric( sub("\\((.+),.*", "\\1", neutral_range) ),
                  upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", neutral_range) ))                  # extracting lower and upper bounds of 8th interval
cnv_mat[ neutral_range['lower'] <= cnv_mat & cnv_mat <= neutral_range['upper']] = 1                     # all values falling in the 8th interval are set to 1 (the centroid)
cnv_mat = abs(cnv_mat-1)                                                                                # to make sure neutral values are assigned black color; amplification or deletion contribute synergistically

s_obj$avg_cnv = 0                                                                                       # average CNV across all chromosomes
i_ = 0                                                                                                  # chromosome counter
for(chr_ in unique(chrs_))
{
  s_obj[[paste0('avg_cnv_',chr_)]] = rowMeans(cnv_mat[, names(genes_) %in% chr_])
  s_obj$avg_cnv = s_obj$avg_cnv + s_obj[[paste0('avg_cnv_',chr_)]]
  i_ = i_+1
}
s_obj$avg_cnv = s_obj$avg_cnv/i_

# STEP 4: Plotting average CNV for each cell per chromosome ####

pdf(file = paste0('e.pdf'), width = 40, height = 40)

  # average CNV per chromosome
  p_ = list()     # plot list
  for(chr_ in unique(chrs_))
  {
    col_ = paste0('avg_cnv_',chr_)      # column of meta.data containing average cnv for chr
    breaks_ = c(min(s_obj[[col_]]), mean(range(s_obj[[col_]])), max(s_obj[[col_]]))
    p_[[chr_]] =  FeaturePlot(s_obj, features = col_, pt.size = .4, order = T, coord.fixed = F)+
                  theme(legend.spacing.y = unit(0.5, 'cm'), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), text = element_text(size = 22, face = 'bold'), plot.title = element_text(size = 25, face = 'bold'))+
                  labs(title = chr_, color = 'CNV level')+
                  scale_color_viridis(option = 'A', breaks = breaks_, labels = round(breaks_,2))
  }

  # average CNV across entire genome
  breaks_ = c(min(s_obj[['avg_cnv']]), mean(range(s_obj[['avg_cnv']])), max(s_obj[['avg_cnv']]))
  p_[[length(p_)+1]] =  FeaturePlot(s_obj, features = 'avg_cnv', pt.size = .4, order = T, coord.fixed = F)+
                        theme(legend.spacing.y = unit(0.5, 'cm'), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), text = element_text(size = 22, face = 'bold'), plot.title = element_text(size = 25, face = 'bold'))+
                        labs(title = 'Genome', color = 'CNV level')+
                        scale_color_viridis(option = 'A', breaks = breaks_, labels = round(breaks_,2))
  do.call(what = grid.arrange, args = c(grobs = p_, nrow = 5))
graphics.off()

