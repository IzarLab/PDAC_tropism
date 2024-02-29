# This scripts carries out pseudo-time pathway expression analysis on CAF compartment

library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)
library(dbscan)
library(xlsx)
library(monocle3)
library(parallel)

# Loading data ####

load('compartments/CAF.RData')

# subsetting

s_objs = subset(s_objs, subset = cell_type %in% c('hypoactive myofibroblast','myofibroblast','degrading myofibroblast','nascent proaxonogenic fibroblast', 'proaxonogenic fibroblast'))

# Converting Seurat to Monocle ####

DefaultAssay(s_objs) = 'RNA'                                              # automatically adds RNA assay's counts slot
cds = as.cell_data_set(x = s_objs)                                        # automatically adds UMAP and PCA from integrated assay as that's the only assay with dim reductions
rowData(cds)$gene_short_name = rowData(cds)$gene_name = rownames(cds)     # NEEDED for gene-based plotting
cds = estimate_size_factors(cds)                                          # NEEDED for plot_genes_in_pseudotime

# Trajectory works better with lower granularity like cell types instead of clusters ####

cds@clusters$UMAP$cluster_result$optim_res$membership = as.character(colData(cds)$cell_type)
cds@clusters$UMAP$clusters = colData(cds)$cell_type
cds@clusters$UMAP$partitions = factor(x = rep(1, ncol(cds)), levels = '1')
names(cds@clusters$UMAP$cluster_result$optim_res$membership) = names(cds@clusters$UMAP$clusters) = names(cds@clusters$UMAP$partitions) = colnames(cds)

# Learning trajectory #####

cds = learn_graph(cds, use_partition = TRUE, verbose = FALSE, close_loop = F,
                  learn_graph_control = list(minimal_branch_len = 15, nn.k = 15))

# ordering cells in pseudo-time time

cds = order_cells(cds, root_cells = colnames(cds)[cds@colData$cell_type %in% 'hypoactive myofibroblast'][100])      # choose the starting cells (roots)
s_objs$pseudotime = pseudotime(cds)

# plotting

pdf(file = 'pseudotime_CAF.pdf', width = 8, height = 8)

p_ =  plot_cells(cds,cell_size = 1,
                 color_cells_by = "pseudotime",
                 label_cell_groups = FALSE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE,
                 graph_label_size = 4, trajectory_graph_color = "cyan", trajectory_graph_segment_size = 1.5)+
      theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
            legend.text = element_text(face = 'bold', family = 'Helvetica', size = 10), legend.title = element_text(face = 'bold', family = 'Helvetica', size = 15),
            text = element_text(size = 10), plot.background = element_blank(), panel.background = element_blank())

Idents(s_objs) = s_objs$cell_type
p_ =  FeaturePlot(s_objs, features = 'pseudotime', label = T, label.size = 5, pt.size = 1.5, repel = T, label.color = 'green3')+
      theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
            legend.text = element_text(face = 'bold', family = 'Helvetica', size = 10), legend.title = element_text(face = 'bold', family = 'Helvetica', size = 15),
            text = element_text(size = 10), plot.background = element_blank(), panel.background = element_blank())+
      labs(color = 'Pseudotime', title = NULL)+
      geom_segment(data =  p_$layers[[3]]$data, aes(x = source_prin_graph_dim_1, xend = target_prin_graph_dim_1,
                                                    y = source_prin_graph_dim_2, yend = target_prin_graph_dim_2), linewidth = 1.2, color = 'cyan')+
      scale_color_viridis(option = 'C', breaks = c(min(s_objs$pseudotime),mean(s_objs$pseudotime),max(s_objs$pseudotime)), labels = c('h. myCAF', 'd. myCAF', 'pCAF'))

p_$layers[[2]]$aes_params$fontface = 'bold'
p_$layers[[2]]$aes_params$family = 'Helvetica'
plot(p_)

# Pseudo-time pathway expression analysis ####

# MSigDB/KEGG/GO markers; check with enrichr's pathways

EMT_ = c('ABI3BP','ACTA2','ADAM12','ANPEP','APLP1','AREG','BASP1','BDNF','BGN','BMP1','CADM1','CALD1','CALU','CAP2','CAPG',
         'CD44','CD59','CDH11','CDH2','CDH6','COL11A1','COL12A1','COL16A1','COL1A1','COL1A2','COL3A1','COL4A1','COL4A2',
         'COL5A1','COL5A2','COL5A3','COL6A2','COL6A3','COL7A1','COL8A2','COMP','COPA','CRLF1','CCN2','CTHRC1','CXCL1','CXCL12',
         'CXCL6','CCN1','DAB2','DCN','DKK1','DPYSL3','DST','ECM1','ECM2','EDIL3','EFEMP2','ELN','EMP3','ENO2','FAP','FAS',
         'FBLN1','FBLN2','FBLN5','FBN1','FBN2','FERMT2','FGF2','FLNA','FMOD','FN1','FOXC2','FSTL1','FSTL3','FUCA1','FZD8','GADD45A',
         'GADD45B','GAS1','GEM','GJA1','GLIPR1','COLGALT1','GPC1','GPX7','GREM1','HTRA1','ID2','IGFBP2','IGFBP3','IGFBP4','IL15',
         'IL32','IL6','CXCL8','INHBA','ITGA2','ITGA5','ITGAV','ITGB1','ITGB3','ITGB5','JUN','LAMA1','LAMA2','LAMA3','LAMC1','LAMC2',
         'P3H1','LGALS1','LOX','LOXL1','LOXL2','LRP1','LRRC15','LUM','MAGEE1','MATN2','MATN3','MCM7','MEST','MFAP5','MGP','MMP1',
         'MMP14','MMP2','MMP3','MSX1','MXRA5','MYL9','MYLK','NID2','NNMT','NOTCH2','NT5E','NTM','OXTR','PCOLCE','PCOLCE2','PDGFRB',
         'PDLIM4','PFN2','PLAUR','PLOD1','PLOD2','PLOD3','PMEPA1','PMP22','POSTN','PPIB','PRRX1','PRSS2','PTHLH','PTX3','PVR','QSOX1',
         'RGS4','RHOB','SAT1','SCG2','SDC1','SDC4','SERPINE1','SERPINE2','SERPINH1','SFRP1','SFRP4','SGCB','SGCD','SGCG','SLC6A8',
         'SLIT2','SLIT3','SNAI2','SNTB1','SPARC','SPOCK1','SPP1','TAGLN','TFPI2','TGFB1','TGFBI','TGFBR3','TGM2','THBS1','THBS2','THY1',
         'TIMP1','TIMP3','TNC','TNFAIP3','TNFRSF11B','TNFRSF12A','TPM1','TPM2','TPM4','VCAM1','VCAN','VEGFA','VEGFC','VIM','WIPF1','WNT5A')

pCAF = c('BOC', 'CACNB2', 'COL4A4', 'EPHA3', 'EPHA6', 'NFASC', 'NFIB', 'NRP1', 'NRXN3', 'PAK3', 'PTCH1', 'SCN7A', 'SEMA6D', 'SLIT2',
         'SPTBN1', 'TRPC4')

myCAF = c("CACNA1G","VCL","GUCY1B1","MYLK","ITGB5","PXN","MYL6","SORBS1","CACNA1I","MYL9","MYL12A","MYL10","MYL7","ACTA2","ALDH2",
          "MYL12B","SORBS3","CALD1","MYH11","ANXA1","DYSF","TLN1","PDE5A","TPM1","TPM3","PAK1","GUCY1A2","ACTG2","LMOD1","GUCY1A1",
          "TPM4","TRIM72","MYL11","PAK2","CAV3","ANXA2","MYL6B","CACNA1H","ANXA6","TPM2","CALM1","ITGA1","MYL5","MYH11")

#___________________ plotting EMT ___________________

genes = EMT_
genes = genes[genes %in% rownames(s_objs)]
s_objs@meta.data$avg_emt = colMeans(s_objs[['RNA']]@data[genes,,drop = F])

# feature plot

p_ = FeaturePlot(s_objs, features = 'avg_emt', order = T,
            label = T, label.size = 5, repel = T,
            min.cutoff = 'q50', max.cutoff = 'q99',
            pt.size = 1, cols = c('grey90','red3'))+
theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
labs(title = 'EMT')
p_[[1]]$layers[[2]]$aes_params$fontface = 'bold'
plot(p_)

# pseudo-time

dt_ = data.frame(y = s_objs$avg_emt, x = s_objs$pseudotime)
p_ =  ggplot(data = dt_, aes(x = x, y = y))+
      theme(plot.background = element_blank(), panel.background = element_blank(), axis.line = element_line(color = 'black'),
            axis.title = element_text(face = 'bold', size = 15),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(face = 'bold', size = 15, hjust = .5),
            legend.position = 'none')+
      labs(x = 'Pseudotime', y = 'Expression', title = 'EMT')+
      geom_point(aes(color = x), size = .8)+
      geom_smooth(method = 'loess', formula = y ~ x, se = F, linewidth = 2, span = .9)+
      scale_color_viridis(option = 'B', aesthetics = 'color')
plot(p_)

#___________________ plotting pCAF ___________________

genes = pCAF
genes = genes[genes %in% rownames(s_objs)]
s_objs@meta.data$avg_pcaf = colMeans(s_objs[['RNA']]@data[genes,,drop = F])

# feature plot

p_ = FeaturePlot(s_objs, features = 'avg_pcaf', order = T,
                 label = T, label.size = 5, repel = T,
                 min.cutoff = 'q50', max.cutoff = 'q99',
                 pt.size = 1, cols = c('grey90','red3'))+
  theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
  labs(title = 'Axon guidance/axonogenesis')
p_[[1]]$layers[[2]]$aes_params$fontface = 'bold'
plot(p_)

# pseudo-time

dt_ = data.frame(y = s_objs$avg_pcaf, x = s_objs$pseudotime)

p_ =  ggplot(data = dt_, aes(x = x, y = y))+
      theme(plot.background = element_blank(), panel.background = element_blank(), axis.line = element_line(color = 'black'),
            axis.title = element_text(face = 'bold', size = 15),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(face = 'bold', size = 15, hjust = .5),
            legend.position = 'none')+
      labs(x = 'Pseudotime', y = 'Expression', title = 'Axon guidance/axonogenesis program')+
      geom_point(aes(color = x), size = .8)+
      geom_smooth(method = 'loess', formula = y ~ x, se = F, linewidth = 2, span = .9)+
      scale_color_viridis(option = 'B', aesthetics = 'color')
plot(p_)

#___________________ plotting myCAF ___________________

genes = myCAF
genes = genes[genes %in% rownames(s_objs)]
s_objs@meta.data$avg_mycaf = colMeans(s_objs[['RNA']]@data[genes,,drop = F])

# feature plot

p_ = FeaturePlot(s_objs, features = 'avg_mycaf', order = T,
                 label = T, label.size = 5, repel = T,
                 min.cutoff = 'q40', max.cutoff = 'q99',
                 pt.size = 1, cols = c('grey90','red3'))+
      theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
      labs(title = 'Muscle contraction')
p_[[1]]$layers[[2]]$aes_params$fontface = 'bold'
plot(p_)

# pseudo-time

dt_ = data.frame(y = s_objs$avg_mycaf, x = s_objs$pseudotime)

p_ =  ggplot(data = dt_, aes(x = x, y = y))+
      theme(plot.background = element_blank(), panel.background = element_blank(), axis.line = element_line(color = 'black'),
            axis.title = element_text(face = 'bold', size = 15),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(face = 'bold', size = 15, hjust = .5),
            legend.position = 'none')+
      labs(x = 'Pseudotime', y = 'Expression', title = 'Muscle contraction')+
      geom_point(aes(color = x), size = .8)+
      geom_smooth(method = 'loess', formula = y ~ x, se = F, linewidth = 2, span = .9)+
      scale_color_viridis(option = 'B', aesthetics = 'color')
plot(p_)

graphics.off()
