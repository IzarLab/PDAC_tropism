# This scripts plots UMAP projection of pathways active in major CAF subtypes

library(Seurat)
library(SeuratObject)
options(Seurat.object.assay.version = "v4")
library(SeuratWrappers)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)
library(xlsx)
library(monocle3)
library(parallel)


# Loading data ####

load(file = 'compartments/CAF.RData')
DefaultAssay(s_objs) = 'RNA'

# Reading in pathway signatures ####

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

cxcl14 = 'CXCL14'

iCAF = c('ABCA1','ABI1','ACVR1B','ACVR2A','ADM','ADORA2B','ADRM1','AHR','APLNR','AQP9','ATP2A2','ATP2B1','ATP2C1','AXL','BDKRB1',
           'BEST1','BST2','BTG2','C3AR1','C5AR1','CALCRL','CCL17','CCL2','CCL20','CCL22','CCL24','CCL5','CCL7','CCR7','CCRL2','CD14',
           'CD40','CD48','CD55','CD69','CD70','CD82','CDKN1A','CHST2','CLEC5A','CMKLR1','CSF1','CSF3','CSF3R','CX3CL1','CXCL10',
           'CXCL11','CXCL6','CXCL9','CXCR6','CYBB','DCBLD2','EBI3','EDN1','EIF2AK2','EMP3','ADGRE1','EREG','F3','FFAR2','FPR1','FZD5',
           'GABBR1','GCH1','GNA15','GNAI3','GP1BA','GPC3','GPR132','GPR183','HAS2','HBEGF','HIF1A','HPN','HRH1','ICAM1','ICAM4',
           'ICOSLG','IFITM1','IFNAR1','IFNGR2','IL10','IL10RA','IL12B','IL15','IL15RA','IL18','IL18R1','IL18RAP','IL1A','IL1B','IL1R1',
           'IL2RB','IL4R','IL6','IL7R','CXCL8','INHBA','IRAK2','IRF1','IRF7','ITGA5','ITGB3','ITGB8','KCNA3','KCNJ2','KCNMB2','KIF1B',
           'KLF6','LAMP3','LCK','LCP2','LDLR','LIF','LPAR1','LTA','LY6E','LYN','MARCO','MEFV','MEP1A','MET','MMP14','MSR1','MXD1','MYC',
           'NAMPT','NDP','NFKB1','NFKBIA','NLRP3','NMI','NMUR1','NOD2','NPFFR2','OLR1','OPRK1','OSM','OSMR','P2RX4','P2RX7','P2RY2',
           'PCDH7','PDE4B','PDPN','PIK3R5','PLAUR','PROK2','PSEN1','PTAFR','PTGER2','PTGER4','PTGIR','PTPRE','PVR','RAF1','RASGRP1',
           'RELA','RGS1','RGS16','RHOG','RIPK2','RNF144B','ROS1','RTP4','SCARF1','SCN1B','SELE','SELL','SELENOS','SEMA4D','SERPINE1',
           'SGMS2','SLAMF1','SLC11A2','SLC1A2','SLC28A2','SLC31A1','SLC31A2','SLC4A4','SLC7A1','SLC7A2','SPHK1','SRI','STAB1','TACR1',
           'TACR3','TAPBP','TIMP1','TLR1','TLR2','TLR3','TNFAIP6','TNFRSF1B','TNFRSF9','TNFSF10','TNFSF15','TNFSF9','TPBG','VIP')

# Plotting ####

pdf(file = 'CAF_pathways.pdf', width = 8, height = 8)

#___________________ plotting pCAF ___________________

genes = pCAF
genes = genes[genes %in% rownames(s_objs)]
s_objs@meta.data$avg_pcaf = colMeans(s_objs[['RNA']]@data[genes,,drop = F])

# feature plot

p_ = FeaturePlot(s_objs, features = 'avg_pcaf', order = T,
                 label = T, label.size = 5, repel = F,
                 min.cutoff = 'q60', max.cutoff = NA,
                 pt.size = 1, cols = c('grey90','red3'))+
  theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
  labs(title = 'Axon guidance/axonogenesis')
p_[[1]]$layers[[2]]$aes_params$fontface = 'bold'
plot(p_)

#___________________ plotting myCAF ___________________

genes = myCAF
genes = genes[genes %in% rownames(s_objs)]
s_objs@meta.data$avg_mycaf = colMeans(s_objs[['RNA']]@data[genes,,drop = F])

# feature plot

p_ = FeaturePlot(s_objs, features = 'avg_mycaf', order = T,
                 label = T, label.size = 5, repel = F,
                 min.cutoff = 'q60', max.cutoff = NA,
                 pt.size = 1, cols = c('grey90','red3'))+
  theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
  labs(title = 'Muscle contraction')
p_[[1]]$layers[[2]]$aes_params$fontface = 'bold'
plot(p_)

#___________________ plotting iCAF ___________________

genes = iCAF
genes = genes[genes %in% rownames(s_objs)]
s_objs@meta.data$avg_icaf = colMeans(s_objs[['RNA']]@data[genes,,drop = F])

# feature plot

p_ = FeaturePlot(s_objs, features = 'avg_icaf', order = T,
                 label = T, label.size = 5, repel = F,
                 min.cutoff = 'q60', max.cutoff = NA,
                 pt.size = 1, cols = c('grey90','red3'))+
  theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
  labs(title = 'Inflammatory response')
p_[[1]]$layers[[2]]$aes_params$fontface = 'bold'
plot(p_)

#___________________ plotting CXCL14 ___________________

genes = cxcl14
genes = genes[genes %in% rownames(s_objs)]
s_objs@meta.data$avg_cxcl14 = colMeans(s_objs[['RNA']]@data[genes,,drop = F])

# feature plot

p_ = FeaturePlot(s_objs, features = 'avg_cxcl14', order = T,
                 label = T, label.size = 5, repel = F,
                 min.cutoff = 'q20', max.cutoff = 'q99',
                 pt.size = 1, cols = c('grey90','red3'))+
  theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold', size = 20, family = 'Helvetica'), axis.text = element_blank(), axis.ticks = element_blank())+
  labs(title = 'CXCL14')
p_[[1]]$layers[[2]]$aes_params$fontface = 'bold'
plot(p_)

graphics.off()
