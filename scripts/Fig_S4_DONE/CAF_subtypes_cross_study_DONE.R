# This scripts plots overlap between cross-study CAF subtypes

library(Seurat)
library(SeuratObject)
options(Seurat.object.assay.version = "v4")
library(SeuratWrappers)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)

# Loading data ####

load(file = 'compartments/CAF.RData')
DefaultAssay(s_objs) = 'RNA'
s_objs = subset(s_objs, subset = cell_type %in% c('myofibroblast','CXCL14 fibroblast','proaxonogenic fibroblast','inflammatory fibroblast'))

# Markers from PMID: 35902743 ####

Immod_ = c('SLC22A3','XKR4','ANKRD29','SLCO2B1','LAMA3','ABCC3','LAMC2','GRIN2B','RBM47','NOL4','CP','KEL','ZNF804B','TNC','ACTB','TMEM108',
           'TMEM178B','CCL21','ABCB11','SLCO2A1','IL15','FDCSP','MUSK','PLA2G4C','ATP8A1','ADGRL3','LIFR','NPY1R','ARHGAP15','CTSS','RASGEF1B',
           'BIRC3','NRG2','JUN','LPAR1','EXOC3L4','PTPRF','CR2','CHL1','EGR1','ANO9','SLCO1A2','ZFP36L2','OSMR','EDNRB','PTMA','PLD5','EHBP1L1',
           'TIMM23B','TRAF1','TAGLN','RASGEF1A','EEF1A1','CACNA2D3','ADRA1A','S1PR3','IER3','PLCXD3','SLC26A7','IRF8','NFAM1','PDE4B','SORL1',
           'ACHE','SLC26A3','TPM4','CDH1','CTSH','PAPLN','SDK1','DAPK2','ACTG1','DEPTOR','CYSLTR2','DTNA','COL27A1','CXCL12','TNFAIP2','NR4A1',
           'LINC01197','CR1','CSF2RB','VCAM1','TMSB4X','LMF1','OCA2','RPS9','TNFRSF1B','THBS1','LDLR','PTPRT','MYO16','EBF3','TLR1','C1RL-AS1',
           'FOS','SERPINB9','COL23A1','GNA14','PKP2','FTH1','SAT1','NCEH1','TNFAIP3','JUND','TPT1','KIAA1671','PNISR','PCOLCE2','FGF7','ITIH5',
           'UBC','HDAC9','CYR61','ADAMTSL1','GRIA4','GARNL3','IL4R','SPNS2','NBEAL2','ZFP36','PNRC1','PARP14','TXNIP','STRIP2','SVEP1','FADS1',
           'PLXNA4','SLC2A3','LINGO1','C7','PRRC2C','DAAM1','CLSTN3','CCL19','INO80D','ATP8B4','ALPK1','NOVA1','COL4A4','PITPNC1','PCDH11X',
           'RPS11','APBB1IP','TNFSF10','RPS6KL1','PDE1C','GAPDH','NFKBIA','CD82','TLE2','SRRM2','TMEM176B','C3','PYHIN1','HMCN2','RPS8','RPL13',
           'RPL10','IRF3','SMAP2','ADGRE5','PER1','PLEC','ZDHHC14','PHPT1','EPHB1','FAM189A1','RPLP0','TXNRD1','LINC00092','SPATA6L','EEF2',
           'ABHD17C','ART4','NFKB2','ARRDC3','HNMT','CLU','SPIC','TOM1L1','FTL','MYL6','RPS2','DOCK8','COLQ','GAK','KLHL25','SORCS1','PLXNB2',
           'SLC9A9','AHNAK','DDX24','IL34','EPAS1','FLNB','SLIT3','EPHA2','ETS2','FAM20A')

myo_ = c('ADAMTS12','CASC15','POSTN','NTM','LINC01429','NREP','PDGFC','LEF1','NUAK1','COL1A1','KIF26B','NOX4','FN1','SULF1','COL1A2','WNT5A',
         'COL3A1','COL11A1','CDH11','NKD1','DOCK4','PLPP4','MMP11','ADAMTS14','ADAMTS6','FAP','RUNX2','RUNX1','MGAT5','SNTB1','KIAA1549L',
         'CTHRC1','LINC00578','RNF144A','ENC1','SYTL2','ITGA1','DCBLD1','COL10A1','CALD1','CARMN','CHST11','PDZD2','ANTXR1','GREM1','INHBA',
         'NPR3','GRIP1','SLC6A6','FBXO32','FGD6','SALL4','KCND2','ITGA11','MIR181A1HG','ACTA2','LAMA4','APBB2','EDNRA','FUT8','BICD1','MBOAT2',
         'PALM2-AKAP2','SUGCT','HIP1','VSNL1','ENTPD1','SGIP1','EEPD1','KANK4','FRMD5','PPFIBP1','FOXP1','ADAM19','SIPA1L1','FARP1','PTK7',
         'NHSL1','VCAN','HMGA2','EPSTI1','CDK6','SPATS2L','PALLD','STAMBPL1','RASGRF2','MYH9','ARHGAP31','CDKL5','TPM1','ATXN1','PTPRE','ZEB1',
         'FAM168A','ST6GAL2','COL5A1','FNDC1','PLXNC1','EIF4G3','ANO1','LRIG3','TCF4','HOXB3','APBA2','GULP1','HECW1','TLN2','SPIN1','IRS1',
         'SPON2','NXN','TSC22D1','TENM4','GRIK2','NPAS2','STX7','BCAT1','PRICKLE1','ZNF521','ZNF532','KLHL2','ITGB5','TNFRSF19','ARMC9','TNS3',
         'GFPT1','MRVI1','WNK1','MANBA','TBL1XR1','MIR4435-2HG','DNAJC15','SCN8A','TWIST1','COG6','PCED1B','C9orf3','MAP3K4','SMC6','DIO2',
         'TTC3','ZMYM4','SAMD3','STARD4-AS1','WWC1','UNC5B','LINC01060','DGKI','BBX','SSH1','KCNQ1OT1','GXYLT2','GPR63','F13A1','WLS','PLS3',
         'SLC24A2','EPC1','KIFAP3','COL7A1','VGLL4','PDZRN3','PRDM1','COL8A1','ISM1','ZNF609','CLCN3','ADAM22','ACVR1','TMEM45A','PRDM6','ETV6',
         'AEBP1','TANC2','SRPK2','DDR2','WIPF1','PRKD1','SPATS2','RFX8','MSC-AS1','MMP14','ZNF292','IGFBP5','ANO4','MAML3','RAI14','HECW2',
         'PRR5L','MDFIC','PCCA','ETV1','SIPA1L3','CAMK4','VEZT','UHRF2','GUCY1A2','PHF21A','GPC6','NBAT1')

notrop = c('SCN7A','NFIA','C7','PID1','C1orf21','MAMDC2','CLMN','PREX2','MTUS1','ADAMTS9-AS2','KCNIP1','LAMA2','EBF1','ABCA6','NID1','EPHA3',
           'IL1RAPL1','TMEM132C','SPTBN1','ADAMTSL3','NEGR1','AC016831.7','SLC9A9','MIR99AHG','ZBTB20','SRPX','ABCA8','TGFBR3','ABCA10','PTEN',
           'ZBTB16','RHOBTB3','SLIT2','PDK4','FREM1','SOX6','CACNA1D','ABI3BP','HMGCLL1','AOX1','MAPK10','SSH2','KAZN','ARHGAP10','AFF3',
           'ARHGAP6','ABLIM1','PTPRG','ADGRD1','SPARCL1','FKBP5','ABCA9','ANKS1B','COL21A1','FRMD3','IMMP2L','CELF2','ADD3','CCNH','HAND2-AS1',
           'DSCAML1','TFPI','NR2F2-AS1','BOC','ADGRB3','PDE1A','MKLN1','NFIB','PBX1','FBLN5','CPED1','HIF3A','PIK3R1','TENM2','COL4A4','SESN3',
           'ITPR1','DLG2','FBLN1','SSBP2','GPHN','ADAMTS3','SAMHD1','KCTD3','LINC01088','NEURL1B','RUNX1T1','GPC5','SOX5','ADGRL2','AFF1','NOVA1',
           'PARD3','GFRA1','CCND3','FAM13A','MFSD6','RGL1','SETBP1','GHR','DYNC1I1','CCDC102B','SPATA6','NSF','DCN','LDLRAD3','DCLK1','RNF13',
           'RORA','NLGN1','CECR2','FOXO3','COL4A1','UTRN','PIAS1','COL4A3','CACNB2','COL4A2','FIGN','FMNL2','SOBP','EGFR','TGFBR2','FAM102B','DPYD',
           'FAM135A','PPP1CB','IRAK3','TMEM144','MGST1','LAMB1','ADCY3','PODN','ADH1B','PLSCR4','ABTB2','ALDH1A1','PAK3','NBEA','ITM2B','TNRC6B',
           'TRERF1','STXBP4','PRKCH','SLC8A1','PITPNC1','CALCRL','C1S','MT1X','JADE1','PTPRK','PDE7B','ACACB','OGFRL1','CLIP4','GAB1','PLEKHA5',
           'SMIM14','TXNIP','NCOA1','RBMS1','SDCCAG8','STK24','FOXP2','ABCA9-AS1','COLGALT2','KDM6A','PTPN13','MATN2','CD47','ARHGAP26','MBP','DANT2',
           'KIF5C','CNKSR2','KCNT2','ABCC9','RALGAPA1','IL6ST','ECHDC2','SYNPO2','TBC1D5','ELMO1','GFOD1','MARCH2','FLRT2','KLF12','RBM26','KCNN3',
           'PDE3A','NAALADL2','INSR','AUH','PLCL2','CDON','PTPN12','PRKAG2','ADAMTS15','OAF','KLHL13')

adh_ = c('NFATC2','EMP1','MIR222HG','SAMD4A','LMNA','GPRC5A','MMP19','MEDAG','NFATC1','TSC22D2','LRRFIP1','RFX2','PFKP','PTPRJ','ANKRD28','CAV1',
         'TEX26-AS1','CDH2','ANXA2','CTNNAL1','SLC19A2','CRY1','CNN1','SYN3','ANXA5','TES','LHFPL2','LMCD1','ERRFI1','UGP2','LMCD1-AS1','IQCJ-SCHIP1',
         'ACSL4','ZSWIM6','DDAH1','PTPN1','ABL2','ESYT2','GFPT2','ATP13A3','BAIAP2','GLIS3','ERCC1','CD44','ENAH','SERPINE1','CLIC4','ATP10A','FNIP2',
         'MYOF','NEDD9','FOSL1','RTN4','COBL','MYH10','FOSB','KDM6B','CAPN2','ANXA1','YWHAZ','RGCC','EGFR','HIF1A','SH3RF1','ELL2','KLF6','WEE1','S100A10',
         'P4HA3','HOMER1','TRIB1','ADAM12','ITGA5','SLC7A1','KCNMA1','TUBB6','HRH1','GEM','GPR176','PCGF5','MICAL2','PER2','ST6GALNAC5','DOK5','LOX','COL12A1',
         'TIMP3','ACTN4','FLNB','CRIM1','PMEPA1','EFHD2','S100A6','PXDC1','MYO1B','TAOK3','MLF1','UAP1','CORO1C','TIPARP','ITGB1','DENND5A','MEF2A','RNF149',
         'MKL1','MYO9B','EXT1','GLUD1','FNDC3B','GADD45B','NUP153','HMGA1','FGFR1','FHL2','ARC','CBLB','MBNL2','ACTN1','FAM155A','FSIP1','GATAD2A','CREB5','SIK3',
         'KALRN','CDK17','ITPKC','PDLIM5','PSME4','NCS1','MBNL1','CAMK1D','SGK1','RAB11A','PLAT','RPS6KA3','FLNC','YWHAG','ITGAV','MAP2K3','IL1R1','ADAM17',
         'HIF1A-AS2','PITPNM2','TNFRSF12A','HEG1','IQGAP1','VGLL3','TP53BP2','PLEKHA7','AXL','PHF20','KIF1B','LPP','MCL1','PPP1R12B','ASAP2','GLS','ATP1A1','CLIP1',
         'TRIO','LINC00968','ATP2B4','PTPN14','RYBP','HECTD2','ADAM9','RAI14','UCHL3','DDX21','XYLT1','GREB1L','TGFBR1','KPNA4','DMD','PDZRN4','IPO7','SEPT9','YBX3',
         'RAMP1','IGF2BP2','LDHA','ATP11A','PCDH7','PRKCA-AS1','SERTAD2','DENND4A','WISP1','IL6R','HIPK2','HIVEP2','ACLY','PLAUR','ANO6','CTNNA1','PDLIM3','CAV2',
         'MSRB3','MTHFD1L','COBLL1','DPYSL3')

Immod_ = Immod_[Immod_ %in% rownames(s_objs)]
myo_ = myo_[myo_ %in% rownames(s_objs)]
notrop = notrop[notrop %in% rownames(s_objs)]
adh_ = adh_[adh_ %in% rownames(s_objs)]

s_objs$Immunomodulatory = colMeans(s_objs[['RNA']]@data[Immod_,])
s_objs$Myofibroblastic = colMeans(s_objs[['RNA']]@data[myo_,])
s_objs$Neurotropic = colMeans(s_objs[['RNA']]@data[notrop,])
s_objs$Adhesive = colMeans(s_objs[['RNA']]@data[adh_,])

# Violin plot ####

p_ = VlnPlot(s_objs, features = c('Immunomodulatory','Myofibroblastic','Neurotropic','Adhesive'),
             pt.size = 0,
             group.by = 'cell_type',
             sort = T,
             ncol = 3,
             combine = T,
             log = F)
p_[[1]] = p_[[1]] + labs(y = 'Signature expression', x = '')+
                    stat_summary(fun = mean, geom ='point', size = 10, colour = 'black', shape = 95)
p_[[2]] = p_[[2]] + labs(y = 'Signature expression', x = '')+
                    stat_summary(fun = mean, geom ='point', size = 10, colour = 'black', shape = 95)
p_[[3]] = p_[[3]] + labs(y = 'Signature expression', x = '')+
                    stat_summary(fun = mean, geom ='point', size = 10, colour = 'black', shape = 95)
p_[[4]] = p_[[4]] + labs(y = 'Signature expression', x = '')+
                    stat_summary(fun = mean, geom ='point', size = 10, colour = 'black', shape = 95)

ggsave('CAF_subtypes_cross_study.pdf', device = 'pdf', width = 24, height = 14, units = 'in', dpi = 600, plot = p_)
