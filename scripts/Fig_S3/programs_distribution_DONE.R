# This script generates malignant program distribution in PDAC subtypes

library(ggplot2)
library(viridis)
library(gridExtra)
library(Seurat)

# Loading malignant program scores ####

load('Misc/malignant_programs_cancer.RData')

# Plotting bar plots ####

pdf(file = 'malignant_programs_distribution.pdf', width = 42, height = 15)

# forming data frame for ggplot2

dt_ = list()
for(i_ in 1:length(patient_signature_scores))
{
  sgn_ = unlist(patient_signature_scores[[i_]])
  dt_[[i_]] = data.frame(sample = names(sgn_), score = sgn_, program = names(patient_signature_scores[i_]), stringsAsFactors = T)
}
dt_ = do.call(dt_, what = rbind)
dt_$condition = as.factor(conditions_[as.character(dt_$sample)])

# setting up colors

lineage_ = c('Basal-like', 'Classical', 'Progenitor', 'Ciliated')
states_ = as.character(unique(dt_$program[!dt_$program %in% lineage_]))

cols_ = c('#4660D5FF','navy','#FABB39FF','#E7490CFF','chartreuse2','pink3','chocolate4','mistyrose2','#4774ECFF'
          ,'#27BEE9FF','aquamarine3','red','#3E9BFEFF','seashell3','#476DE5FF','violetred2','#CF2E05FF'
          ,'cornsilk4','tan','#EF5911FF','maroon2','purple1','lightcoral','#EECF3AFF','#30123BFF','forestgreen'
          ,'#28BDEAFF','lightsalmon2','mediumseagreen','lightblue4','deeppink4','#A8FB39FF','#4771E9FF')
names(cols_) = states_

# plotting


dt_tmp = lims_ = option_ = NULL
dt_tmp = dt_[!dt_$program %in% lineage_,]; lims_ = c(0, 3.5)

dt_liver = dt_tmp[dt_tmp$condition %in% 'liver',]
p_ =  ggplot(data = dt_liver, aes(fill = program, y = score, x = sample, label = round(score,2)))+
      theme(legend.position = 'bottom', panel.background = element_rect(fill = 'white'), panel.grid = element_line(color = 'grey90'), text = element_text(face = 'bold'))+
      labs(title = 'LIV-MR', x = 'Patient', y = 'Program enrichment', fill = 'Malignant program')+
      geom_bar(position="stack", stat="identity", colour = "black")+
      geom_text(size = 10, position = position_stack(vjust = 0.5), color = 'white', fontface = 2)+
      scale_y_continuous(limits = lims_)
if(i_ == 1){ p_ = p_ + scale_fill_viridis(discrete = T) }else{ p_ = p_ + scale_fill_manual(values = cols_)}

dt_lung = dt_tmp[dt_tmp$condition %in% 'lung',]
q_ =  ggplot(data = dt_lung, aes(fill = program, y = score, x = sample, label = round(score,2)))+
      theme(legend.position = 'bottom', panel.background = element_rect(fill = 'white'), panel.grid = element_line(color = 'grey90'), text = element_text(face = 'bold'))+
      labs(title = 'LUN-MR', x = 'Patient', y = '', fill = 'Malignant program')+
      geom_bar(position="stack", stat="identity", colour="black")+
      geom_text(size = 10, position = position_stack(vjust = 0.5), color = 'white', fontface = 2)+
      scale_y_continuous(limits = lims_)
q_ = q_ + scale_fill_manual(values = cols_)

  plot(p_ + q_)
graphics.off()

# Significance testing ####

p_val = numeric()
for(sgn_ in states_)
{
  liv_scores = dt_$score[dt_$program %in% sgn_ & dt_$condition %in% 'liver']
  lng_scores = dt_$score[dt_$program %in% sgn_ & dt_$condition %in% 'lung']
  if(length(liv_scores) < 4 | length(lng_scores) < 4){ cat('not possible for ',sgn_,'\n'); next() }
  p_val[sgn_] = wilcox.test(liv_scores, lng_scores, alternative = 'two.sided')$p.value
}
print(sort(p_val[p_val <= 0.1]))

