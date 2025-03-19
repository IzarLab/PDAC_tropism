# This script plots Extended Data Fig. 3a; it generates dot plot of malignant program distribution in PDAC subtypes.
# The algorithm for computing signature score per cell population is given in Fig_5/i.R script.

library(ggplot2)
library(viridis)
library(gridExtra)
library(Seurat)
options(Seurat.object.assay.version = "v3")

# Loading malignant program scores ####

load('Misc/malignant_programs_cancer.RData')

# Plotting bar plots ####

pdf(file = 'a.pdf', width = 17, height = 15)

# forming data frame for ggplot2

dt_ = list()
for(i_ in 1:length(patient_signature_scores))     # patient_signature_scores a list with signature (aka program) enrichemnt of malignant cell state programs for each patient
{
  sgn_ = unlist(patient_signature_scores[[i_]])
  dt_[[i_]] = data.frame(sample = names(sgn_), score = sgn_, program = names(patient_signature_scores[i_]), stringsAsFactors = T)
}
dt_ = do.call(dt_, what = rbind)
dt_$condition = as.factor(conditions_[as.character(dt_$sample)])

# making sure each patient has an enrichment for all programs

all_samples = unique(dt_$sample)
for(p_ in unique(dt_$program))
{
  smpls_ = dt_$sample[dt_$program %in% p_]
  absent_smpls = all_samples[!all_samples %in% smpls_]
  if(length(absent_smpls) == 0) next()
  tmp_ = data.frame(sample = absent_smpls, score = 0, program = p_, condition = NA, row.names = absent_smpls)
  for(s_ in absent_smpls){ tmp_[s_,"condition"] = as.character(dt_$condition[dt_$sample %in% s_])[1] }
  dt_ = rbind(dt_, tmp_)
}

# scaling each pathway score between 0 and 1 across patients

for(p_ in unique(dt_$program))
{
  idxs_ = dt_$program %in% p_
  scores_ = c(0, dt_$score[idxs_])     # 0 is added to make sure min score is not assigned 0 unless it is indeed 0
  scores_ = (scores_ - min(scores_))/diff(range(scores_))
  dt_$score[idxs_] = scores_[-1]
}

# plotting

dt_liver = dt_[dt_$condition %in% 'liver',]

p_ =  ggplot(data = dt_liver, aes(y = program, x = sample, size = score, color = score))+
      theme(panel.background = element_blank(), panel.grid = element_line(color = 'black',linewidth = .05),
            legend.position = 'right', legend.text = element_text(size = 15), legend.key = element_blank(),
            legend.title = element_text(size = 15), legend.spacing.y = unit(.6, 'cm'),
            axis.text.x = element_text(angle = -20),
            text =  element_text(face = 'bold', size = 18, family = 'Helvetica'))+
      labs(y = NULL, x = NULL, color = 'Program enrichment', size = 'Program enrichment')+
      geom_point()+
      scale_size_continuous(range = c(1,12), limits = c(0,1), breaks = c(0,.5,1))+
      scale_color_gradient2(low = 'blue4', mid = 'yellow3', high = 'red', midpoint = .5, breaks = c(0 ,0.5,1), labels = c('low','mid','high'))+
      guides(color = guide_colorbar(ticks.colour = NA, label.position = 'right'))
plot(p_)

dt_lung = dt_[dt_$condition %in% 'lung',]
p_ =  ggplot(data = dt_lung, aes(y = program, x = sample, size = score, color = score))+
      theme(panel.background = element_blank(), panel.grid = element_line(color = 'black',linewidth = .05),
            legend.position = 'right', legend.text = element_text(size = 15), legend.key = element_blank(),
            legend.title = element_text(size = 15), legend.spacing.y = unit(.6, 'cm'),
            axis.text.x = element_text(angle = -20),
            text =  element_text(face = 'bold', size = 18, family = 'Helvetica'))+
      labs(y = NULL, x = NULL, color = 'Program enrichment', size = 'Program enrichment')+
      geom_point()+
      scale_size_continuous(range = c(1,12), limits = c(0,1), breaks = c(0,.5,1))+
      scale_color_gradient2(low = 'blue4', mid = 'yellow3', high = 'red', midpoint = .5, breaks = c(0 ,0.5,1), labels = c('low','mid','high'))+
      guides(color = guide_colorbar(ticks.colour = NA, label.position = 'right'))
plot(p_)

graphics.off()

# Significance testing ####

p_val = numeric()
for(sgn_ in unique(dt_$program))
{
  liv_scores = dt_$score[dt_$program %in% sgn_ & dt_$condition %in% 'liver']
  lng_scores = dt_$score[dt_$program %in% sgn_ & dt_$condition %in% 'lung']
  if(length(liv_scores) < 4 | length(lng_scores) < 4){ cat('not possible for ',sgn_,'\n'); next() }
  p_val[sgn_] = wilcox.test(x = liv_scores, y = lng_scores, alternative = 'two.sided', exact = F)$p.value
}
p_val = p_val[p_val <= 0.01]
if(0 < length(p_val)){ print(sort(p_val)) }

