# This script plots Fig. 5a,c,d,f and g (expression of some important markers in exocrine compartment).

library(Seurat)
library(SeuratObject)
options(Seurat.object.assay.version = "v3")
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)

# Loading data ####

load(file = 'compartments/cancer.RData'); cancer_ = s_objs; rm(s_objs)
DefaultAssay(cancer_) = 'RNA'

load(file = 'compartments/exocrine.RData'); exo_ = s_objs; rm(s_objs)
DefaultAssay(exo_) = 'RNA'

load(file = 'compartments/exocrine_differentiating-cancer.RData'); exo_can = s_objs; rm(s_objs)
DefaultAssay(exo_can) = 'RNA'
exo_can$cell_type[exo_can$cell_type %in% c('PanIN','differentiating cancer')] = 'PanIN/cancer'

# Markers ####

diff_ = c('TFF1','TFF2','TFF3','KLF4','SOX9')     # TFF1-3 are also markers of classical subtype

acinar_ = c('CPB1','PRSS1')

ductal_ = c('CFTR')

DCLK1_ = c('DCLK1')

lesion_ = c('MUC1', 'MUC3A', 'MUC5AC', 'MUC5B', 'MUC16','TFF1','TFF2','TFF3','KLF5','KRT7','KRT19')

# MSigDB HALL marks markers; check with enrichr's pathways
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

# Plotting ####

pdf('a&c-d&f-g.pdf', width = 7, height = 7)

# Differentiating cancer

genes = diff_
genes = unique(genes[genes %in% rownames(cancer_)])
cancer_@meta.data$avg = colMeans(cancer_[["RNA"]]@data[genes,,drop = F])
s_ = cancer_
Idents(s_) = s_$cell_type

FeaturePlot(object = s_, features = 'avg',
            order = T, repel = F,
            label = T, label.size = 6, label.color = 'blue',
            pt.size = .7,
            coord.fixed = F, raster = F, cols = c('grey90','red3'),
            max.cutoff = 'q99', min.cutoff = 'q60',
            ncol = 1)+
theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
labs(title = 'Differentiation')

# Acinar

genes = acinar_
genes = unique(genes[genes %in% rownames(exo_)])
exo_@meta.data$avg = colMeans(exo_[["RNA"]]@data[genes,,drop = F])
s_ = exo_
Idents(s_) = s_$cell_type

FeaturePlot(object = s_, features = 'avg',
            order = T, repel = F,
            label = T, label.size = 6, label.color = 'blue',
            pt.size = .7,
            coord.fixed = F, raster = F, cols = c('grey90','red3'),
            max.cutoff = 'q99', min.cutoff = 'q10',
            ncol = 1)+
theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
labs(title = 'Acinar')

# Ductal

genes = ductal_
genes = unique(genes[genes %in% rownames(exo_)])
exo_@meta.data$avg = colMeans(exo_[["RNA"]]@data[genes,,drop = F])
s_ = exo_
Idents(s_) = s_$cell_type

FeaturePlot(object = s_, features = 'avg',
            order = T, repel = F,
            label = T, label.size = 6, label.color = 'blue',
            pt.size = .7,
            coord.fixed = F, raster = F, cols = c('grey90','red3'),
            max.cutoff = 'q99', min.cutoff = 'q10',
            ncol = 1)+
theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
labs(title = 'Ductal')

# DCLK1

genes = DCLK1_
genes = unique(genes[genes %in% rownames(exo_)])
exo_@meta.data$avg = colMeans(exo_[["RNA"]]@data[genes,,drop = F])
s_ = exo_
Idents(s_) = s_$cell_type

FeaturePlot(object = s_, features = 'avg',
            order = T, repel = F,
            label = T, label.size = 6, label.color = 'blue',
            pt.size = .7,
            coord.fixed = F, raster = F, cols = c('grey90','red3'),
            max.cutoff = 'q90', min.cutoff = NA,
            ncol = 1)+
theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
labs(title = 'DCLK1+')

# Lesion

genes = lesion_
genes = unique(genes[genes %in% rownames(exo_can)])
exo_can@meta.data$avg = colMeans(exo_can[["RNA"]]@data[genes,,drop = F])
s_ = exo_can
Idents(s_) = s_$cell_type

FeaturePlot(object = s_, features = 'avg',
            order = T, repel = F,
            label = T, label.size = 6, label.color = 'blue',
            pt.size = .7,
            coord.fixed = F, raster = F, cols = c('grey90','red3'),
            max.cutoff = 'q95', min.cutoff = NA,
            ncol = 1)+
theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
labs(title = 'Lesion')

# EMT

genes = EMT_
genes = unique(genes[genes %in% rownames(exo_can)])
exo_can@meta.data$avg = colMeans(exo_can[["RNA"]]@data[genes,,drop = F])
s_ = exo_can
Idents(s_) = s_$cell_type

FeaturePlot(object = s_, features = 'avg',
            order = T, repel = F,
            label = T, label.size = 6, label.color = 'blue',
            pt.size = .7,
            coord.fixed = F, raster = F, cols = c('grey90','red3'),
            max.cutoff = '.45', min.cutoff = 'q60',
            ncol = 1)+
theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
labs(title = 'EMT')

graphics.off()

