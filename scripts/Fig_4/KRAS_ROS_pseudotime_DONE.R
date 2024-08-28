# This script carries out pseudo-time KRAS/ROS signature expression analysis in exocrine compartment

library(ggplot2)
library(viridis)
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
set.seed(42)

# Reading in data ####

load('compartments/exocrine_differentiating-cancer.RData')

# Converting Seurat to Monocle ####

# breaking the loop between cancer and ADM

cells_ = read.delim('Misc/loop_cells.txt', header = F, quote = '')$V1
s_objs = subset(s_objs, cells = colnames(s_objs)[!colnames(s_objs) %in% cells_])

# conversion

DefaultAssay(s_objs) = 'RNA'                                              # as.cell_data_set automatically adds RNA assay's counts slot
cds = as.cell_data_set(x = s_objs)                                        # automatically adds UMAP and PCA from integrated2 assay as that's the only assay with dim reductions
rowData(cds)$gene_short_name = rowData(cds)$gene_name = rownames(cds)     # NEEDED for gene-based plotting
cds = estimate_size_factors(cds)                                          # NEEDED for plot_genes_in_pseudotime

# Learning trajectory #####

# trajectory works better with lower granularity like cell types instead of clusters

cds@clusters$UMAP$cluster_result$optim_res$membership = as.character(colData(cds)$cell_type)
cds@clusters$UMAP$clusters = colData(cds)$cell_type
cds@clusters$UMAP$partitions = factor(x = rep(1, ncol(cds)), levels = '1')
names(cds@clusters$UMAP$cluster_result$optim_res$membership) = names(cds@clusters$UMAP$clusters) = names(cds@clusters$UMAP$partitions) = colnames(cds)

# learning cell trajectory

cds = learn_graph(cds, use_partition = TRUE, verbose = FALSE, close_loop = F,
                  learn_graph_control = list(minimal_branch_len = 18))

# Pseudo-time analysis ####

# ordering cells in pseudo-time time

cds = order_cells(cds, root_cells = colnames(s_objs)[s_objs$cell_type %in% 'acinar'][1])
s_objs$pseudotime = pseudotime(cds)

# MSigDB markers

ros_ = c('ABCC1','ATOX1','CAT','CDKN2D','EGLN2','ERCC2','FES','FTL','G6PD','GCLC','GCLM','GLRX','GLRX2',
         'GPX3','GPX4','GSR','LAMTOR5','HHEX','HMOX2','IPCEF1','JUNB','LSP1','MBP','MGST1','MPO','MSRA',
         'NDUFA6','NDUFB4','NDUFS2','NQO1','OXSR1','PDLIM1','PFKP','PTPA','PRDX1','PRDX2','PRDX4','PRDX6',
         'PRNP','SBNO2','SCAF4','SELENOS','SOD1','SOD2','SRXN1','STK25','TXN','TXNRD1','TXNRD2')

kras_up = c('ABCB1','ACE','ADAM17','ADAM8','ADAMDEC1','AKAP12','AKT2','ALDH1A2','ALDH1A3','AMMECR1','ANGPTL4',
            'ANKH','ANO1','ANXA10','APOD','ARG1','ATG10','AVL9','BIRC3','BMP2','BPGM','BTBD3','BTC','C3AR1',
            'CA2','CAB39L','CBL','CBR4','CBX8','CCL20','CCND2','CD37','CDADC1','CFB','CFH','CFHR2','CIDEA',
            'CLEC4A','CMKLR1','CPE','CROT','CSF2','CSF2RA','CTSS','CXCL10','CXCR4','DCBLD2','DNMBP','DOCK2',
            'DUSP6','ADGRL4','EMP1','ENG','EPB41L3','EPHB2','EREG','ERO1A','ETS1','ETV1','ETV4','ETV5','EVI5',
            'F13A1','F2RL1','CCSER2','FBXO4','FCER1G','FGF9','FLT4','FUCA1','G0S2','GABRA3','GADD45G','GALNT3',
            'GFPT2','GLRX','GNG11','GPNMB','ADGRA2','GPRC5B','GUCY1A1','GYPC','HBEGF','HDAC9','H2BC3','HKDC1',
            'HOXD11','HSD11B1','ID2','IGF2','IGFBP3','IKZF1','IL10RA','IL1B','IL1RL2','IL2RG','IL33','IL7R','INHBA',
            'IRF8','ITGA2','ITGB2','ITGBL1','JUP','KCNN4','KIF5C','KLF4','LAPTM5','LAT2','LCP1','LIF','LY96','MAFB',
            'MALL','MAP3K1','MAP4K1','MAP7','MMD','MMP10','MMP11','MMP9','MPZL2','MTMR10','MYCN','NAP1L2','NGF','NIN',
            'NR0B2','NR1H4','NRP1','PCP4','PCSK1N','PDCD1LG2','PECAM1','PEG3','PIGR','PLAT','PLAU','PLAUR','PLEK2',
            'PLVAP','PPBP','PPP1R15A','PRDM1','PRKG2','PRRX1','PSMB8','PTBP2','PTCD2','PTGS2','PTPRR','RABGAP1L','RBM4',
            'RBP4','RELN','RETN','RGS16','SATB1','SCG3','SCG5','SCN1B','SDCCAG8','SEMA3B','SERPINA3','PRELID3B','SLPI',
            'SNAP25','SNAP91','SOX9','SPARCL1','SPON1','SPP1','SPRY2','ST6GAL1','STRN','TFPI','TLR8','TMEM100','TMEM158',
            'TMEM176A','TMEM176B','TNFAIP3','TNFRSF1B','TNNT2','TOR1AIP2','TPH1','TRAF1','TRIB1','TRIB2','TSPAN1','TSPAN13',
            'TSPAN7','USH1C','USP12','VWA5A','WDR33','WNT7A','YRDC','ZNF277','ZNF639')

DefaultAssay(s_objs) = 'RNA'

genes = ros_
genes = genes[genes %in% rownames(s_objs)]
s_objs@meta.data$avg_ros = colMeans(s_objs[['RNA']]@data[genes,,drop = F])

genes = kras_up
genes = genes[genes %in% rownames(s_objs)]
s_objs@meta.data$avg_kras = colMeans(s_objs[['RNA']]@data[genes,,drop = F])

# plotting

pdf(width = 17, height = 15, file = 'pseudotime_KRAS_ROS_.pdf')

s_ = s_objs
s_@meta.data = s_@meta.data[order(s_$pseudotime),]              # removing gaps on x-axis

# ros

dt_ = data.frame(y = s_$avg_ros, x = order(s_$pseudotime))      # removing gaps on x-axis
dt_ = dt_[dt_$y > 0,]
ggplot(data = dt_, aes(x = x, y = y))+
theme(plot.background = element_blank(), panel.background = element_blank(), axis.line = element_line(color = 'black'),
      text = element_text(face = 'bold', family = 'Helvetica', size = 25),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = 'bold', size = 15, hjust = .5),
      legend.position = 'right')+
labs(x = 'Pseudotime', y = 'Expression', title = 'ROS', color = 'Pseudotime')+
geom_point(aes(color = x), size = 2)+
geom_smooth(method = 'loess', formula = y ~ x, se = F, linewidth = 4, span = 5, color = 'cyan')+
scale_color_viridis(option = 'B', breaks = c(min(dt_$x), mean(dt_$x), max(dt_$x)),
                   labels = c('acinar','ductal/ductal-like','PanIN/cancer'), name = 'Pseudotime')

# kras

dt_ = data.frame(y = s_$avg_kras, x = order(s_$pseudotime))      # removing gaps on x-axis
dt_ = dt_[dt_$y > 0,]
ggplot(data = dt_, aes(x = x, y = y))+
theme(plot.background = element_blank(), panel.background = element_blank(), axis.line = element_line(color = 'black'),
      text = element_text(face = 'bold', family = 'Helvetica', size = 25),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = 'bold', size = 15, hjust = .5),
      legend.position = 'right')+
labs(x = 'Pseudotime', y = 'Expression', title = 'KRAS', color = 'Pseudotime')+
geom_point(aes(color = x), size = 2)+
geom_smooth(method = 'loess', formula = y ~ x, se = F, linewidth = 4, span = 5, color = 'cyan')+
scale_color_viridis(option = 'B', breaks = c(min(dt_$x), mean(dt_$x), max(dt_$x)),
                    labels = c('acinar','ductal/ductal-like','PanIN/cancer'), name = 'Pseudotime')

graphics.off()
