# This scripts plots Fig. 5i; it generates a dot plot of canonical pathways of exocrine compartment.

library(ggplot2)
library(gridExtra)
library(patchwork)
library(cowplot)
library(xlsx)
library(Seurat)
options(Seurat.object.assay.version = "v3")

# Reading in data ####

load('compartments/cancer.RData'); cancer_ = s_objs;
DefaultAssay(cancer_) = 'RNA'

load('compartments/exocrine.RData'); exo_ = s_objs; rm(s_objs); gc()
DefaultAssay(exo_) = 'RNA'

# Pathways of interest ####

# Att.: these are taken form Enrichr-curated repos: GO_Biological_Process_2023, KEGG_2021 and MSigDB_Hallmark_2023

pathws_ = list(
TGFB_prod = c('SMAD3','TGFB2','SMAD4','ITGB6','XCL1','SERPINF2','LILRB1','FN1','HIF1A','LTBP1','PTGS2','LRRC32','ITGB8','ITGAV','LGALS9','CD46'),
pan_sec = c('AMY1A','AMY1B','AMY1C','SLC12A2','ATP1B4','ATP1B3','CCK','ATP1B2','ATP1B1','CHRM3','SCT','PRSS2',
                         'PRSS1','ADCY1','ADCY5','ATP1A4','ADCY4','ATP1A3','ADCY3','ATP1A2','ADCY2','ATP1A1','ADCY9',
                         'ADCY8','ADCY7','SLC9A1','ADCY6','PRSS3','CPA3','CELA2A','PRKCG','CPA2','CELA2B','CPA1','CEL',
                         'RAB27B','PRKCB','PRKCA','KCNMA1','GNAQ','GNAS','CFTR','PNLIP','CELA3A','CELA3B','RAB3D','SLC26A3',
                         'CPB2','AMY2A','CPB1','AMY2B','PLA2G3','PLA2G5','CLCA2','CLCA1','CLCA4','RHOA','PLA2G2E','PLA2G2F',
                         'PLA2G2C','PLA2G2D','PLA2G2A','TPCN2','BST1','CTRL','FXYD2','PNLIPRP2','PNLIPRP1','ITPR2','ITPR3',
                         'ITPR1','RAP1A','CA2','RAP1B','CCKAR','TRPC1','CD38','PLA2G12A','PLA2G12B','CTRB2','ATP2B4','ATP2B3',
                         'ATP2B2','SCTR','ATP2B1','RAB11A','PLCB4','KCNQ1','PLA2G10','CTRB1','PLCB2','PLCB3','PLCB1','RYR2',
                         'PLA2G1B','ATP2A3','ATP2A2','SLC4A2','ATP2A1','SLC4A4','RAB8A','RAC1'),
prot_dig_absr = c('COL28A1','COL24A1','ATP1B4','ATP1B3','ATP1B2','ATP1B1','COL1A2','SLC9A3','COL1A1','COL5A2','XPNPEP2','COL5A1',
                  'COL5A3','COL9A2','COL9A1','COL9A3','COL20A1','COL16A1','PRSS2','PRSS1','COL12A1','SLC1A1','ATP1A4','SLC1A5',
                  'ATP1A3','ATP1A2','PRCP','ATP1A1','PRSS3','CPA3','CELA2A','SLC15A1','CPA2','SLC38A2','CELA2B','SLC6A19','CPA1',
                  'COL25A1','SLC16A10','PGA4','PGA3','PGA5','ACE2','COL2A1','COL6A1','MEP1B','COL6A3','COL6A2','COL21A1','COL6A5',
                  'MEP1A','COL6A6','CELA3A','CELA3B','KCNE3','CPB2','COL17A1','CPB1','COL13A1','KCNN4','KCNJ13','COL26A1','SLC7A7',
                  'COL3A1','SLC7A8','SLC7A9','CTRL','FXYD2','COL7A1','COL22A1','KCNK5','COL18A1','COL14A1','ELN','SLC3A1','SLC3A2',
                  'DPP4','COL10A1','SLC36A2','SLC36A1','COL27A1','SLC36A4','SLC36A3','MME','COL23A1','CTRB2','COL4A1','COL4A3','KCNQ1',
                  'COL4A2','COL4A5','COL8A1','COL4A4','CTRB1','COL8A2','COL4A6','COL19A1','COL15A1','COL11A2','SLC8A1','SLC8A2','SLC8A3',
                  'COL11A1'),
IL2_STAT5 = c('IL1R2','TTC39B','RRAGD','ITIH5','CYFIP1','IFITM3','AHNAK','FGL2','PNP','MUC1','FLT3LG','MYC','WLS','GBP4','ST3GAL4',
              'PRKCH','RABGAP1L','BATF3','GABARAPL1','XBP1','GADD45B','CDCP1','UCK2','PLSCR1','PLPP1','CDC42SE2','RNH1','IKZF2','ENO3',
              'IKZF4','GATA1','CD79B','ALCAM','TNFSF11','TNFSF10','NCOA3','IL10','ANXA4','IL13','IL10RA','RHOH','EMP1','LIF','PTCH1','CDC6',
              'FAM126B','RHOB','P4HA1','BHLHE40','ETFBKMT','TLR7','UMPS','IL18R1','CTSZ','PTH1R','SHE','SOCS1','SOCS2','PENK','SH3BGRL2','CTLA4',
              'MAP6','NFIL3','SNX9','CCR4','ICOS','IFNGR1','GSTO1','DHRS3','SERPINB6','ADAM19','CCNE1','PTRH2','F2RL2','CD48','PLIN2','SLC29A2',
              'CD44','EEF1AKMT1','NRP1','ECM1','AHCY','CDKN1C','SERPINC1','PTGER2','SNX14','CAPG','CST7','PDCD2L','NCS1','ENPP1','PLAGL1','DCPS',
              'EOMES','IL4R','PHLDA1','APLP1','S100A1','PUS1','PRAF2','HUWE1','HIPK2','FAH','BATF','CXCL10','MAPKAPK2','CD86','CSF2','CD83','CSF1',
              'ODC1','CD81','RGS16','GPR65','SLC1A5','AHR','LCLAT1','NDRG1','NT5E','CCND3','CCND2','CAPN3','SLC39A8','SPP1','DENND5A','MAP3K8',
              'SWAP70','IGF2R','CISH','PRNP','TRAF1','HOPX','CKAP4','MAFF','SELP','MXD1','COL6A1','AMACR','SELL','BCL2L1','DRC1','GUCY1B1','SLC2A3',
              'AGER','RORA','GPR83','IGF1R','TGM2','GLIPR2','SCN9A','TNFRSF8','TNFRSF9','PIM1','TNFRSF4','POU2F1','GPX4','NFKBIZ','ETV4','TNFRSF1B',
              'TNFRSF21','BMP2','ITGAE','IRF4','IL3RA','IRF8','BCL2','IRF6','ITGA6','PLEC','ARL4A','BMPR2','ABCB1','NOP2','PHTF2','HK2','SYNGR2','CA2',
              'TNFRSF18','CASP3','TIAM1','MYO1E','TWSG1','MYO1C','KLF6','P2RX4','IL2RA','SYT11','IL2RB','LTB','SPRY4','SMPDL3A','GALM','FURIN','COCH',
              'IL1RL1','SPRED2','LRIG1','ITGAV','LRRC8C'),
antigen_proc_pres = c('CTSB','CD74','HSP90AA1','KLRC3','HSPA6','HSPA5','HSPA1L','KLRC4','HSPA8','HSPA2','TAP2','HSPA4','TAP1','TAPBP','CD4','RFX5','IFNG','HSPA1A','HLA-DRA',
                      'HLA-DQB1','KLRD1','CALR','HSPA1B','TNF','KIR2DL1','KIR2DL2','KIR2DL3','KIR2DL4','CTSS','CTSL','B2M','HLA-DOB','KIR2DL5A','CIITA','HLA-DQA2','HLA-DOA',
                      'HLA-DQA1','PDIA3','HLA-DRB5','KIR2DS1','NFYB','KIR2DS2','NFYC','KIR2DS3','KIR2DS4','KIR2DS5','HLA-C','RFXANK','KIR3DL1','CD8B2','KIR3DL2','HLA-A','HLA-B',
                      'KIR3DL3','HLA-G','HLA-E','HLA-F','CREB1','CD8B','RFXAP','CD8A','CANX','LGMN','PSME2','PSME3','HLA-DPB1','PSME1','HLA-DRB4','HLA-DRB3','HLA-DRB1','HSP90AB1',
                      'IFI30','HLA-DMB','HLA-DPA1','KLRC1','KLRC2','NFYA','HLA-DMA'),
IFNg_sign = c('GBP6','SP110','MX1','RSAD2','MX2','PSMB9','PML','NUP93','LAP3','MTHFD2','IFI27','TRIM14','RBCK1','IDO1','IFITM2','IFITM3','FGL2','NOD1',
              'UBE2L6','PSMB10','PNP','APOL6','HERC6','ST3GAL5','SLAMF7','B2M','GBP4','EPSTI1','BATF2','ST8SIA4','PLA2G4A','HLA-A','HLA-B','HLA-G','IFI44',
              'PARP14','LATS2','PARP12','VAMP8','PLSCR1','BANK1','SERPING1','PSME2','SSPN','PSME1','VAMP5','CFB','HLA-DRB1','CFH','IFI35','HIF1A','IFIT2','IFI30',
              'IFIT1','ICAM1','IFIT3','RNF213','MT2A','HELZ2','KLRK1','TNFSF10','CMKLR1','NCOA3','GPR18','PNPT1','IL15','IL10RA','NMI','BST2','RAPGEF6','IL18BP',
              'CD40','BTG1','SECTM1','FPR1','SOCS3','SOCS1','RIPK2','CD38','RIPK1','HLA-DQA1','P2RY14','STAT1','STAT2','GCH1','STAT3','STAT4','PSMA2','PSMA3','CMPK2',
              'PELI1','MARCHF1','CDKN1A','OGFR','PTGS2','IFIH1','PSMB8','NAMPT','PSMB2','METTL7B','HLA-DMA','JAK2','PTPN1','ZBP1','PTPN2','IL4R','CD74','TOR1B','TAP1',
              'VCAM1','CXCL10','CXCL11','IFI44L','ZNFX1','TRAFD1','FAS','LY6E','CD69','PTPN6','CD86','SAMD9L','RNF31','CMTR1','DDX58','AUTS2','ARID5B','TDRD7','ISG15',
              'NFKB1','ISG20','MYD88','SELP','XAF1','IFNAR2','CD274','LGALS3BP','SPPL2A','CSF2RB','DDX60','PIM1','XCL1','EIF2AK2','TAPBP','IRF1','NFKBIA','IL7','IL6',
              'OAS2','IRF4','IRF5','OAS3','IRF2','IRF8','IRF9','IRF7','LCP2','EIF4E3','ARL4A','CXCL9','TNFAIP6','TNFAIP3','NLRC5','TNFAIP2','ADAR','CASP8','CASP7','CASP4',
              'PDE4B','CASP3','ISOC1','CASP1','CIITA','IL15RA','GZMA','LYSMD2','BPGM','WARS1','SOD2','IL2RB','TXNIP','PFKP','C1R','RTP4','C1S','MVP','SRI','USP18','OASL',
              'SAMHD1','SLC25A28','FCGR1A','CCL7','DHX58','CCL5','TRIM25','TRIM26','CCL2','ITGB7','TRIM21','UPP1'),
KRAS = c('ATG10','CIDEA','EREG','MALL','NIN','RBP4','SLPI','ITGBL1','SPARCL1','ZNF277','ADAMDEC1','SATB1','C3AR1','RETN','IL2RG','GLRX','TFPI','ADGRA2','AKT2','RABGAP1L',
         'ST6GAL1','PPBP','GADD45G','NR0B2','AMMECR1','IGF2','BTC','TNNT2','CCSER2','STRN','CFB','SPON1','CFH','PLEK2','IKZF1','GNG11','TSPAN7','EPB41L3','ANKH','SOX9',
         'KCNN4','CMKLR1','TSPAN1','APOD','SCN1B','GABRA3','FUCA1','DUSP6','ACE','JUP','CCL20','IL10RA','EMP1','LIF','PDCD1LG2','SNAP91','MMP11','DCBLD2','WDR33','VWA5A',
         'MMP10','LAT2','IL1B','GPRC5B','CFHR2','ID2','TLR8','SDCCAG8','AVL9','ADGRL4','SERPINA3','TPH1','PTPRR','GYPC','CAB39L','CROT','TOR1AIP2','CTSS','MAP7','CD37',
         'FCER1G','MPZL2','PRELID3B','HSD11B1','IL33','MMD','MAP4K1','TMEM176A','ARG1','G0S2','TMEM176B','WNT7A','LAPTM5','NR1H4','ANO1','TMEM158','ADAM17','F2RL1','ANGPTL4',
         'ENG','HBEGF','NRP1','MTMR10','ERO1A','PIGR','F13A1','HOXD11','PTGS2','PSMB8','AKAP12','CXCR4','GPNMB','CBR4','MAP3K1','PTCD2','GUCY1A1','PRRX1','H2BC3','CXCL10',
         'MYCN','PCP4','CLEC4A','SCG5','ADAM8','SCG3','PLVAP','CSF2','PPP1R15A','RGS16','PRDM1','RELN','CCND2','KIF5C','SPP1','IGFBP3','TRAF1','INHBA','MAFB','DNMBP','YRDC',
         'BIRC3','IL7R','PLAU','PLAT','ETS1','CSF2RA','TMEM100','FGF9','BTBD3','PCSK1N','CDADC1','CBX8','GALNT3','PEG3','GFPT2','ITGA2','ETV1','NAP1L2','MMP9','ETV4','ETV5',
         'TNFRSF1B','BMP2','LCP1','CPE','ZNF639','IRF8','PECAM1','DOCK2','SNAP25','ABCB1','TNFAIP3','LY96','PTBP2','CA2','EVI5','PRKG2','PLAUR','ANXA10','BPGM','NGF','KLF4',
         'ALDH1A3','ALDH1A2','HKDC1','SPRY2','TRIB2','TRIB1','ITGB2','FLT4','USP12','SEMA3B','HDAC9','TSPAN13','CBL','USH1C','RBM4','IL1RL2','FBXO4','EPHB2'),
EMT = c('LUM','CADM1','DKK1','P3H1','COL1A2','SFRP4','COL1A1','MGP','VIM','MYL9','PFN2','MXRA5','PRSS2','TAGLN','LRP1','SERPINE1','COL12A1','SERPINE2','SFRP1','DPYSL3','ABI3BP',
        'SPOCK1','MEST','CCN2','CCN1','GADD45B','FZD8','GADD45A','FERMT2','SLC6A8','SDC1','SNAI2','SPARC','OXTR','ENO2','MYLK','GJA1','LGALS1','ANPEP','SERPINH1','FUCA1','JUN',
        'IL15','TNFRSF12A','MAGEE1','VEGFC','FN1','EMP3','RHOB','VEGFA','CRLF1','COL3A1','MMP14','GAS1','ID2','CALU','FOXC2','SDC4','ELN','TNFRSF11B','AREG','THY1','TIMP3',
        'QSOX1','TIMP1','CAP2','IL32','ADAM12','BGN','COL4A1','COL4A2','FMOD','CD44','LRRC15','ECM1','ECM2','FSTL1','CAPG','FSTL3','COL11A1','CD59','TGFB1','PDGFRB','APLP1',
        'PRRX1','TPM4','TFPI2','WNT5A','TPM2','VCAM1','CXCL12','COL5A2','PLOD1','COL5A1','FAP','FAS','COL5A3','TGFBI','PMP22','SCG2','PPIB','COL16A1','TNC','PLOD3','FBLN2','PVR',
        'PLOD2','FBLN1','LOXL1','LOXL2','FBLN5','COMP','NT5E','CDH6','BASP1','CDH2','GPC1','SLIT2','TPM1','FLNA','MSX1','SPP1','SLIT3','IGFBP2','POSTN','IGFBP4','IGFBP3','INHBA',
        'GEM','LOX','CDH11','COL6A3','COL6A2','MATN2','FBN1','MATN3','FBN2','COPA','NOTCH2','NTM','HTRA1','FGF2','SAT1','TGM2','RGS4','GLIPR1','EFEMP2','CALD1','PMEPA1','EDIL3',
        'ITGA5','MMP1','MMP2','BDNF','GPX7','MMP3','ITGA2','TGFBR3','MFAP5','VCAN','ACTA2','IL6','DAB2','BMP1','COL7A1','MCM7','CXCL8','NNMT','LAMC1','TNFAIP3','LAMC2','CXCL1',
        'NID2','CXCL6','PDLIM4','DST','PCOLCE','PLAUR','DCN','PTHLH','SNTB1','GREM1','COL8A2','PTX3','PCOLCE2','ITGB1','WIPF1','LAMA3','LAMA2','ITGB5','ITGB3','THBS1','LAMA1',
        'THBS2','SGCG','SGCD','SGCB','ITGAV','COLGALT1','CTHRC1'))
names(pathws_) = c('Regulation of TGFB production',
                   'Pancreatic secretion',
                   'Protein digestion and absorption',
                   'IL-2/STAT5 signaling',
                   'Antigen processing and presentation',
                   'IFNG signaling',
                   'KRAS Signaling Up',
                   'EMT')

# Pathway expression score ####

cancer_$cell_type = as.character(cancer_$cell_type)
cancer_$cell_type = as.factor('cancer')                                                   # setting cancer cell types to 'cancer'

cell_types = c('acinar', 'ADM', 'ductal/ductal-like', 'low-grade PanIN', 'PanIN', 'cancer')     # cell types of interest on x-axis (sub-setting?)

exo_@meta.data[,names(pathws_)] = cancer_@meta.data[,names(pathws_)] = 0                  # adding pathways as columns to meta data tables
p_ = list(); j_ = 1
for(c_ in cell_types)
{
  # extracting current cell type
  
  message('processing ', c_)
  
  s_ = tryCatch(subset(exo_, subset = cell_type == c_), error = function(e_){ NULL })                           # is it in exocrine compartment?
  if(is.null(s_)){ s_ = tryCatch(subset(cancer_, subset = cell_type == c_), error = function(e_){ NULL }) }     # or cancer compartment?
  
  # computing pathway (signature) score
  
  for(i_ in 1:length(pathws_))
  {
    pathnm = names(pathws_)[i_]
    message('\t> ', pathnm)
    genes_ = pathws_[[i_]]                                                                                      # genes in the current signature
    
    # signal signature expression
    
    mat_ = matrix(data = 0, ncol = ncol(s_), nrow = length(genes_), dimnames = list(genes_, colnames(s_)))      # a matrix with as many columns as the number of cells in cell type data and as many rows as the signature size (# genes)
    g_ = rownames(s_)[rownames(s_) %in% genes_]                                                                 # signature genes measured in current cell type data
    mat_[g_,] = as.matrix(s_[['RNA']]@data[g_,])                                                                # only measured signature genes are copied, the rest is assigned zero
    pathway_exprs = colMeans(mat_)                                                                              # current signal signature mean expression for each cell
    
    # background signature expression
    
    bkgrnd_exprs = list()
    bkgrnd_genes = rownames(s_[["RNA"]]@data)[! rownames(s_[["RNA"]]@data) %in% genes_]                         # all other genes not present in signal signature comprise background genes
    for(e_ in 1:1e3)
    {
      bkgrnd_genes_e = sample(x = bkgrnd_genes, size = length(genes_))      # sample the same number of genes as in signal signature from background genes
      bkgrnd_exprs[[e_]] = colMeans(s_[["RNA"]]@data[bkgrnd_genes_e,])
    }
    bkgrnd_exprs = do.call(bkgrnd_exprs, what = rbind)                                                          # 1e3 background signature mean expressions for each cell in current cell type data
    
    # identification of significant cells
    
    probs_ = numeric(length = ncol(s_))                                                                         # each cell in current cell type data is assigned a p-value
    names(probs_) = colnames(s_)
    for(cell_ in colnames(s_))
    {
      probs_[cell_] = round(1- ecdf(bkgrnd_exprs[,cell_])(pathway_exprs[cell_]), 5)     # 1- (fraction/frequency of background mean expressions less than or equal to signal mean expression for current cell)
    }
    ###### correction by BH? ######
    # probs_ = p.adjust(probs_, method = 'fdr')
    ###############################
    sig_cells = names(probs_[probs_ <= 0.01])                                                                   # retaining only significant cells (significance level 1%)
    
    if(length(sig_cells) < 10 ){ message('No cell was identified to significantly exprss this program'); next() }
    
    # significant signature expression
    
    pathway_exprs = pathway_exprs[sig_cells] - colMeans(bkgrnd_exprs[, sig_cells, drop = F])                    # pathway (signal signature) score (avoiding batch effect between cells)
    if(s_$histology[1] == 'healthy'){ exo_@meta.data[sig_cells, pathnm] = pathway_exprs }else{ cancer_@meta.data[sig_cells, pathnm] = pathway_exprs }
    
    # visualizing significant cells in current cell type data
    
    p_[[j_]] = DimPlot(s_, cells.highlight = sig_cells)+
               theme(plot.title = element_text(hjust = .5, face = 'bold'), plot.subtitle = element_text(hjust = .5, face = 'bold'),
                     axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())+
               labs(subtitle = pathnm, title = c_) & NoLegend()
    j_ = j_ + 1
  }
}

# Creating dot plot ####

pathway_exprs = rbind(exo_@meta.data[,c("cell_type",names(pathws_))], cancer_@meta.data[,c("cell_type",names(pathws_))])      # concatenating metadata tables by row
pathway_exprs = pathway_exprs[pathway_exprs$cell_type %in% cell_types,]

# formatting data frame for ggplot2

dt_ = list(); i_ = 1
for(c_ in cell_types)
{
  for(pathw_ in names(pathws_))
  {
    exprs_ = pathway_exprs[pathway_exprs$cell_type %in% c_, pathw_]
    dt_[[i_]] = data.frame(cell_type = c_,
                           pathway = pathw_,
                           exprs = mean(exprs_),
                           fraction = round(length(which(0 < exprs_))/length(exprs_)*100, 2)  )
    i_ = i_ + 1
  }
}
dt_ = do.call(dt_, what = rbind)
dt_$cell_type[dt_$cell_type == 'differentiating ductal'] = 'cancer'
dt_$cell_type = factor(dt_$cell_type, levels = unique(dt_$cell_type))
dt_$pathway = factor(dt_$pathway, levels = c('EMT',
                                             'IL-2/STAT5 signaling',
                                             'Regulation of TGFB production',
                                             'Antigen processing and presentation',
                                             'KRAS Signaling Up',
                                             'IFNG signaling',
                                             'Protein digestion and absorption',
                                             'Pancreatic secretion'))

# scaling each pathway score between 0 and 1

for(pathw_ in names(pathws_))
{
  tmp_ = c(0, dt_$exprs[dt_$pathway %in% pathw_])     # 0 is added to make sure min score is not assigned 0 unless it is indeed 0
  tmp_ = (tmp_ - min(tmp_))/diff(range(tmp_))
  dt_$exprs[dt_$pathway %in% pathw_] = tmp_[-1]
}

ggplot(data = dt_, aes(y = pathway, x = cell_type, size = fraction, color = exprs))+
theme(panel.background = element_blank(), panel.grid = element_line(color = 'black',linewidth = .05),
      legend.position = 'right', legend.text = element_text(size = 15), legend.key = element_blank(),
      legend.title = element_text(size = 15), legend.spacing.y = unit(.6, 'cm'),
      axis.text.x = element_text(angle = -20),
      text =  element_text(face = 'bold', size = 18, family = 'Helvetica'))+
labs(y = NULL, x = NULL, color = 'Pathway\nactivity', size = 'Fraction')+
geom_point()+
scale_size_continuous(range = c(1,20), limits=c(0,100), breaks = c(0, 25, 50, 75, 100))+
scale_color_gradient2(low = 'blue4', mid = 'purple3', high = 'red', midpoint = .5, breaks = c(0 ,0.5,1), labels = c('low','mid','high'))+
guides(color = guide_colorbar(ticks.colour = NA, label.position = 'right'),
       size = guide_legend(override.aes = list(size = 1:5)))

ggsave(filename = 'i.pdf', width = 12, height = 10, units = 'in')
