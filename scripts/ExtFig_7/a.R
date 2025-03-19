# This script plots Extended Data Fig. 7a; it plots types of SSV for primary PDAC tumors in the MSK-MET
# data taken from: https://zenodo.org/records/5801902.
# Downloaded files are data_clinical_sample.txt and data_mutations.txt

library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)
library(ggpubr)

# Reading in data ####

## reading in sample data (a patient has one sample associated, from either metastatic or primary tumor) ####

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

dt_smpl = dt_smpl[union(liv_smpls,lng_smpls),c("SAMPLE_ID","METASTATIC_SITE")]

## Reading in mutation data ####

ssm_ = read.delim('Misc/data_mutations.txt', header = T, comment.char = '#', as.is = T, check.names = F, sep = '\t')      # short somatic mutations
ssm_ = ssm_[ssm_$Tumor_Sample_Barcode %in% dt_smpl$SAMPLE_ID,]

# Gene mutation counts ####

dt_ = list()
for(i_ in 1:nrow(ssm_))     # for each measured gene, collects variant classification per recurrence site (liver vs lung)
{
  gene_ssm = ssm_[i_, c("Hugo_Symbol",'Tumor_Sample_Barcode',"Variant_Classification")]
  
  tmp_ = data.frame(gene = gene_ssm$Hugo_Symbol,
                    mutation_type = gene_ssm$Variant_Classification,
                    sample = gene_ssm$Tumor_Sample_Barcode,
                    met_site = dt_smpl$METASTATIC_SITE[dt_smpl$SAMPLE_ID %in% gene_ssm$Tumor_Sample_Barcode])
  
  # if the patient has both liver and lung metastasis
  if(tmp_$met_site == 'liver_lung')
  {
    tmp_ = rbind(tmp_, tmp_)
    tmp_$met_site = c('liver','lung')     # breaks it into one liver and one lung
  }
  
  dt_[[i_]] = tmp_
}
dt_ = do.call(what = rbind, args = dt_)
dt_ = unique(dt_)     # for some genes like CDKN2A, several loci have been reported for per mutation type and per tumor sample

## calculating frequency of each mutation type per gene per recurrence site

met_freq = table(dt_$met_site)
freq_ = list(); i_ = 1
for(g_ in unique(dt_$gene))
{
  for(t_ in unique(dt_$mutation_type))
  {
    for(met_ in c('liver','lung'))
    {
      tmp_ = dt_[dt_$gene %in% g_ & dt_$mutation_type %in% t_ & dt_$met_site %in% met_,]
      freq_[[i_]] = data.frame(gene = g_, mutation_type = t_, met_site = met_, frequency = round(nrow(tmp_)/met_freq[met_]*100,2))
      i_ = i_ + 1
    }
  }
}
freq_ = do.call(freq_, what = rbind)
freq_ = freq_[0 < freq_$frequency,]

freq_$mutation_type[freq_$mutation_type %in% 'Missense_Mutation'] = 'Missense'
freq_$mutation_type[freq_$mutation_type %in% 'Nonsense_Mutation'] = 'Nonsense'
freq_$mutation_type[freq_$mutation_type %in% 'Nonstop_Mutation'] = 'Nonstop'
freq_$mutation_type[freq_$mutation_type %in% 'Translation_Start_Site'] = 'Trans_Str_Site'

freq_$mutation_type = factor(x = freq_$mutation_type, levels = c("Missense","Frame_Shift_Del","Nonsense","Frame_Shift_Ins","Splice_Site","In_Frame_Del","In_Frame_Ins","Splice_Region","Trans_Str_Site","5'Flank","Intron","Nonstop"))


# Plotting ####

pdf(file = 'a.pdf', width = 20, height = 12)

# labels
freq_$label = NA
freq_$label[freq_$gene %in% c('KRAS','TP53','CDKN2A','CDKN2B','SMAD4')] = freq_$gene[freq_$gene %in% c('KRAS','TP53','CDKN2A','CDKN2B','SMAD4')]      # most frequent oncogenic genes in PDAC according to https://bit.ly/44atD71 

# plotting
ggplot(data = freq_, aes(x = mutation_type, y = frequency, shape = met_site, color = met_site, label = label))+
theme(text = element_text(family = 'Helvetica', face = 'bold', size = 20), axis.text.x = element_text(angle = -75, vjust = 2, hjust = 0),
    panel.background = element_blank(), axis.line = element_line(color = 'black'),
    panel.grid = element_line(color = 'grey90'), panel.grid.major.x = element_blank(),
    legend.position = 'bottom')+
labs(x = 'Varient classification', y = 'Frequency', shape = 'Recurrence site', color = 'Recurrence site')+

geom_point(size = 2, alpha = 0.5, position = position_jitter(seed = 42))+
coord_trans(y = 'log')+                 # dose not change the statistics but visualization
scale_color_manual(values = c('red','blue'))+
scale_y_continuous(breaks = c(min(freq_$frequency), 0.5, 1, 5, 10, 15, max(freq_$frequency)) )+
geom_text_repel(position = position_jitter(seed = 42), size = 5.3, fontface = 'bold', show.legend = F,
                segment.color = "black", 
                box.padding = 0.6)+     # increases space between labels and points
geom_vline(xintercept = (1:length(unique(freq_$mutation_type))+0.5), color = '#66023C', linetype = 'dashed', alpha = 0.7)+
guides(color = guide_legend(override.aes = list(size = 7, alpha = 1)))


graphics.off()
