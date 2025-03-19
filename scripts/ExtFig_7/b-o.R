# This script plots Extended Data Fig. 7b-o; it plots FGA and TMB of primary PDAC given SSV, CNV and SV
# mutations of MSK-MET data taken from: https://zenodo.org/records/5801902.
# Downloaded files are data_clinical_sample.txt, data_mutations.txt, data_cna.txt and data_sv.txt

library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)
library(ggpubr)

# Reading in data ####

# reading in sample data (a patient has one sample associated, from either metastatic or primary tumor)

dt_smpl = read.delim('Misc/data_clinical_sample.txt', header = T, comment.char = '#', as.is = T, check.names = F, sep = '\t')
dt_smpl = dt_smpl[dt_smpl$SAMPLE_TYPE %in% 'Primary' & dt_smpl$SUBTYPE %in% 'Pancreatic Adenocarcinoma',]     # only primary samples
rownames(dt_smpl) = dt_smpl$SAMPLE_ID

# subsetting sample data on liver and lung metastatic sites

# dt_smpl = dt_smpl[dt_smpl$MET_SITE_COUNT %in% 1,]     # only one metastasis during the follow-up time? too few lungs will remain (these are primary tumors and should express tropism signatures)
liv_smpls = dt_smpl$SAMPLE_ID[dt_smpl$DMETS_DX_LIVER %in% 'Yes']
lng_smpls = dt_smpl$SAMPLE_ID[dt_smpl$DMETS_DX_LUNG %in% 'Yes']
liv_lng_smpls = intersect(liv_smpls, lng_smpls)
dt_smpl[liv_smpls,"METASTATIC_SITE"] = 'liver'
dt_smpl[lng_smpls,"METASTATIC_SITE"] = 'lung'
dt_smpl[liv_lng_smpls,"METASTATIC_SITE"] = 'liver_lung'

dt_smpl = dt_smpl[union(liv_smpls,lng_smpls),c("SAMPLE_ID","PATIENT_ID","PRIMARY_SITE","SUBTYPE","FGA","TMB_NONSYNONYMOUS","METASTATIC_SITE")]

# reading in mutation data

# point mutation
ssm_ = read.delim('Misc/data_mutations.txt', header = T, comment.char = '#', as.is = T, check.names = F, sep = '\t')
ssm_ = ssm_[ssm_$Tumor_Sample_Barcode %in% dt_smpl$SAMPLE_ID,]

# cnv
cnv_ = read.delim('Misc/data_cna.txt', header = T, comment.char = '#', as.is = T, check.names = F, sep = '\t')      # negative: deletion, positive: amplification
cols_ = c(1,which(colnames(cnv_) %in% dt_smpl$SAMPLE_ID))
cnv_ = cnv_[ ,cols_]

# sv
sv_ = read.delim('Misc/data_sv.txt', header = T, comment.char = '#', as.is = T, check.names = F, sep = '\t')        # in cbioportal, for calculating gene frequency in SV table they considered sites 1 and 2 
sv_ = sv_[sv_$Breakpoint_Type %in% 'PRECISE' & sv_$Sample_ID %in% dt_smpl$SAMPLE_ID,c("Sample_ID","Site1_Hugo_Symbol","Site2_Hugo_Symbol")]

# Gene mutation counts ####

# for each patient sample

genes_ = c('KRAS','TP53','CDKN2A','CDKN2B','SMAD4')     # most frequent oncogenic genes in PDAC according to https://bit.ly/44atD71 
dt_ = list()
for(s_ in unique(dt_smpl$SAMPLE_ID))                    # it's already unique, just in case in the future there are more samples associated with each patient
{
  tmp_ = list()
  
  # for each oncogenic gene
  
  for(g_ in genes_)
  {
    # subsetting on current patient sample and gene
    
    gene_snv = ssm_[ssm_$Tumor_Sample_Barcode %in% s_ & ssm_$Hugo_Symbol %in% g_,, drop = F]
    gene_cnv = cnv_[cnv_$Hugo_Symbol %in% g_, c(1,which(colnames(cnv_) %in% s_)), drop = F]
    gene_sv  = sv_ [sv_$Sample_ID %in% s_ & (sv_$Site1_Hugo_Symbol %in% g_ | sv_$Site2_Hugo_Symbol %in% g_),,drop = F]
    
    tmp_[[g_]] = data.frame(sample = s_,
                            gene = g_,
                            met_site = dt_smpl$METASTATIC_SITE[dt_smpl$SAMPLE_ID %in% s_],
                            FGA = dt_smpl$FGA[dt_smpl$SAMPLE_ID %in% s_],
                            TMB = dt_smpl$TMB_NONSYNONYMOUS[dt_smpl$SAMPLE_ID %in% s_],
                            snv_count = if(0 < nrow(gene_snv)){1}else{0},
                            cnv_count = if(0 < nrow(gene_cnv)){ gene_cnv[,2] }else{0},
                            sv_count  = if(0 < nrow(gene_sv)){1}else{0})
    
    # if the patient has both liver and lung metastasis
    
    if(tmp_[[g_]]$met_site == 'liver_lung')
    {
      tmp_[[g_]] = rbind(tmp_[[g_]], tmp_[[g_]])
      tmp_[[g_]]$met_site = c('liver','lung')     # breaks it into one liver and one lung
    }
  }
  dt_[[s_]] = do.call(what = rbind, args = tmp_)
}
dt_ = do.call(what = rbind, args = dt_)

# Plotting FGA and TMB ####

pdf(file = 'b-o.pdf', width = 15, height = 7.5)

# FGA and TMB given mutated genes

# for each mutation type
for(mut_ in c('SSM','CNV','SV'))
{
  tmp_ = NULL
  for(g_ in genes_)
  {
    if(mut_ == 'SSM'){ tmp_ = dt_[dt_$gene %in% g_ & !(dt_$snv_count %in% 0), c("gene","met_site","snv_count","FGA","TMB")] }
    if(mut_ == 'CNV'){ tmp_ =  dt_[dt_$gene %in% g_ & !(dt_$cnv_count %in% 0), c("gene","met_site","cnv_count","FGA","TMB")] }
    if(mut_ == 'SV') { tmp_ =  dt_[dt_$gene %in% g_ & !(dt_$sv_count %in% 0), c("gene","met_site","sv_count","FGA","TMB")] }
    
    if(nrow(tmp_) == 0 ){ next() }
    if(any(! 5 <= table(tmp_$met_site))){ next() }
    
    p_ = (# FGA
          ggplot(tmp_, aes(x = met_site, y = FGA, fill = met_site))+
          theme(text = element_text(family = 'Helvetica', face = 'bold', size = 20),
                panel.background = element_blank(), panel.grid.major.y = element_line(color = 'grey40'),
                legend.position = 'none')+
          geom_boxplot()+
          stat_compare_means(method = "wilcox.test", fontface = 'bold', size = 9)+  # add p-value using Wilcoxon rank-sum test
          geom_jitter(width = 0, size = 2, color = 'grey40', alpha = 0.4)+
          labs(title = paste0("Fraction genome altered\n(",mut_,':',g_,')'), x = "Metastasis", y = "FGA")+
          scale_fill_manual(values = c('red','blue'))+
          
          # TMB
          ggplot(tmp_, aes(x = met_site, y = TMB, fill = met_site))+
          theme(text = element_text(family = 'Helvetica', face = 'bold', size = 20),
                panel.background = element_blank(), panel.grid.major.y = element_line(color = 'grey40'),
                legend.position = 'none')+
          geom_boxplot()+
          stat_compare_means(method = "wilcox.test", fontface = 'bold', size = 9)+
          geom_jitter(width = 0, size = 2, color = 'grey40', alpha = 0.4)+
          labs(title = paste0("Tumor mutational burden\n(",mut_,':',g_,')'), x = "Metastasis", y = "TMB (%)")+
          scale_fill_manual(values = c('red','blue')))
    plot(p_)
  }
}

# overall FGA and TMB

p_ = (# FGA
      ggplot(dt_, aes(x = met_site, y = FGA, fill = met_site))+
      theme(text = element_text(family = 'Helvetica', face = 'bold', size = 20),
            panel.background = element_blank(), panel.grid.major.y = element_line(color = 'grey40'),
            legend.position = 'none')+
      geom_boxplot()+
      stat_compare_means(method = "wilcox.test", fontface = 'bold', size = 9)+  # add p-value using Wilcoxon rank-sum test
      geom_jitter(width = 0, size = 2, color = 'grey40', alpha = 0.4)+
      labs(title = "Fraction genome altered", x = "Metastasis", y = "FGA")+
      scale_fill_manual(values = c('red','blue'))+
    
      # TMB
      ggplot(dt_, aes(x = met_site, y = TMB, fill = met_site))+
      theme(text = element_text(family = 'Helvetica', face = 'bold', size = 20),
            panel.background = element_blank(), panel.grid.major.y = element_line(color = 'grey40'),
            legend.position = 'none')+
      geom_boxplot()+
      stat_compare_means(method = "wilcox.test", fontface = 'bold', size = 9)+
      geom_jitter(width = 0, size = 2, color = 'grey40', alpha = 0.4)+
      labs(title = "Tumor mutational burden", x = "Metastasis", y = "TMB (%)")+
      scale_fill_manual(values = c('red','blue')))
plot(p_)

graphics.off()

