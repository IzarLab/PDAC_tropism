# This script plots Extended Data Fig. 6c; it plots KM curve for OS of MSK-MET data taken from: https://zenodo.org/records/5801902
# Downloaded files are data_clinical_sample.txt and data_clinical_patient.txt

library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)
library(ggplot2)
library(survival)
library(survminer)

# Reading in data ####

# sample data
dt_smpl = read.delim('Misc/data_clinical_sample.txt', header = T, comment.char = '#', as.is = T, check.names = F, sep = '\t')
rownames(dt_smpl) = dt_smpl$PATIENT_ID

# adding metastatic sites as a new column
cols_ = colnames(dt_smpl)[grepl(pattern = 'DMETS_DX_', x = colnames(dt_smpl))][-1]      # DMETS_DX_UNSPECIFIED adds outliers must be removed
dt_smpl$MET_SITE = NA
for(i_ in 1:nrow(dt_smpl))
{
  dt_smpl$MET_SITE[i_] = paste0(cols_[dt_smpl[i_,cols_] %in% 'Yes'], collapse = ',')
}
dt_smpl = dt_smpl[which(dt_smpl$MET_SITE != ''), ]

# patient data

dt_pnt = read.delim('Misc/data_clinical_patient.txt', header = T, comment.char = '#', as.is = T, check.names = F, sep = '\t')
rownames(dt_pnt) = dt_pnt$PATIENT_ID

# initial metastatic site

dt_smpl = dt_smpl[dt_smpl$MET_SITE_COUNT == 1,]     # only one metastatic site during the follow-up time; among MET_COUNT and MET_SITE_COUNT the latter corresponds to number of Yes's
initial_met_lung = dt_smpl[dt_smpl$MET_SITE %in% 'DMETS_DX_LUNG',]
initial_met_lung_os = dt_pnt[initial_met_lung$PATIENT_ID,]

initial_met_others = dt_smpl[!dt_smpl$MET_SITE %in% c('DMETS_DX_LUNG','DMETS_DX_LIVER','DMETS_DX_DIST_LN'),]
initial_met_others_os = dt_pnt[initial_met_others$PATIENT_ID,]

# Plotting KM curve ####

# preparing data

dt_ = data.frame(patients   = c(initial_met_lung$PATIENT_ID,
                                initial_met_others$PATIENT_ID),
                 
                 OS_MONTHS  = c(initial_met_lung_os$OS_MONTHS,
                                initial_met_others_os$OS_MONTHS),
                 
                 OS_STATUS = c(initial_met_lung_os$OS_STATUS,
                               initial_met_others_os$OS_STATUS),
                 
                 age_death = c(initial_met_lung_os$AGE_AT_DEATH,
                               initial_met_others_os$AGE_AT_DEATH),
                 
                 initial_met = rep(x = c('lung',
                                         'other'),
                                   times = c(nrow(initial_met_lung),
                                             nrow(initial_met_others))),
                 stringsAsFactors = T)
dt_ = dt_[! is.na(dt_$OS_MONTHS),]
dt_$OS_STATUS_code = 0      # 0: living; it is required to use a numerical variable to be able to use ggsurvplot, as it does not work with factor
dt_$OS_STATUS_code[dt_$OS_STATUS %in% '1:DECEASED'] = 1

# create a Survival Object

surv_object = Surv(time = dt_$OS_MONTHS, event = dt_$OS_STATUS_code)

# fit a Kaplan-Meier Model

km_fit = survfit(surv_object ~ initial_met, data = dt_)

# plotting

pdf(file = 'c.pdf', width = 7.5, height = 9)
ggsurvplot(km_fit, data = dt_,
           pval = T, pval.method = T, conf.int = T,
           risk.table = T,
           palette = rev(c('red','blue')),
           xlab = "Time in Months", 
           ylab = "OS Probability", 
           title = "Kaplan-Meier Curve of OS")
graphics.off()
