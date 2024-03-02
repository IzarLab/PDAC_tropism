# This scripts plots QC violins

library(Seurat)
library(SeuratObject)
options(Seurat.object.assay.version = "v4")
library(SeuratWrappers)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)


file_ = 'compartments/healthy_cancer.RData'
load(file = file_)
DefaultAssay(s_objs) = 'RNA'


stress_md = c('BAG3','BLOC1S5-TXNDC5','CALU','DNAJB1','DUSP1','EGR1','FOS','FOSB','HIF1A','HSP90AA1','HSP90AB1','HSP90AB2P','HSP90AB3P',
              'HSP90B1','HSPA1A','HSPA1B','HSPA6','HSPB1','HSPH1','IER2','JUN','JUNB','NFKBIA','NFKBIZ','RGS2','SLC2A3','SOCS3','UBC',
              'ZFAND2A','ZFP36','ZFP36L1')
exAM = c('FOS','HSP90AA1','CCL4','CCL4L2','CCL3','CCL18','CCL3L3','CCL3L1','HSPA1A','HSPA1B','ZFP36','GEM','KLF2','IER5','EGR1','NFKBIZ','JUND','RGS1','JUN','DUSP1','TXNIP','RHOB','ATF3','JUNB','NFKBID')
IFN_ = c('APOBEC3A','CCL2','CCL7','CCL8','CXCL10','CXCL11','CXCL9','DDX58','EIF2AK2','EPSTI1','GBP1','GBP4','HERC5','IDO1','IFI16','IFI35',
        'IFI44','IFI44L','IFI6','IFIH1','IFIT1','IFIT2','IFIT3','IFITM1','IL1RN','IRF7','ISG15','ISG20','LGALS9','LY6E','MT2A','MX1','MX2',
        'NT5C3A','OAS1','OAS2','OAS3','OASL','PLSCR1','RSAD2','SAMD9','SAMD9L','STAT1','TNFSF10','UBE2L6','WARS','XAF1')

stress_md = stress_md[stress_md %in% rownames(s_objs)]
exAM = exAM[exAM %in% rownames(s_objs)]
IFN_ = IFN_[IFN_ %in% rownames(s_objs)]

s_objs$sress = colMeans(s_objs[['RNA']]@data[stress_md,])
s_objs$exAM = colMeans(s_objs[['RNA']]@data[exAM,])
s_objs$ifn_avg = colMeans(s_objs[['RNA']]@data[IFN_,])
s_objs$mt_prcnt = PercentageFeatureSet(s_objs, pattern = "^MT-")

# Violin plot ####

p_ = VlnPlot(s_objs, features = c('sress','exAM','ifn_avg','nFeature_RNA', 'mt_prcnt'),
          pt.size = 0,
          group.by = 'subcompartment',
          sort = T,
          ncol = 3,
          combine = T,
          log = F)
p_[[1]]$layers[[1]]$geom$default_aes$colour = NA

ggsave('CUIMC_snRNAseq_QC.pdf', device = 'pdf', width = 24, height = 14, units = 'in', dpi = 600, plot = p_)

