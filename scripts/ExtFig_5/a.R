# This script plots KM curve for OS of MSK-MET data taken from: https://zenodo.org/records/5801902
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

# patient data

dt_pnt = read.delim('Misc/data_clinical_patient.txt', header = T, comment.char = '#', as.is = T, check.names = F, sep = '\t')
rownames(dt_pnt) = dt_pnt$PATIENT_ID

# initial metastatic site

dt_smpl = dt_smpl[dt_smpl$MET_SITE_COUNT == 1,]     # only one metastatic site during the follow-up time
initial_met_liver = dt_smpl[dt_smpl$DMETS_DX_LIVER %in% 'Yes',]
initial_met_liver_os = dt_pnt[initial_met_liver$PATIENT_ID,]
initial_met_others = dt_smpl[!dt_smpl$DMETS_DX_LIVER %in% 'Yes',]
initial_met_others_os = dt_pnt[initial_met_others$PATIENT_ID,]

# Plotting KM curve ####

# preparing data

dt_ = data.frame(patients   = c(initial_met_liver$PATIENT_ID,
                               initial_met_others$PATIENT_ID),
                 
                 OS_MONTHS  = c(initial_met_liver_os$OS_MONTHS,
                               initial_met_others_os$OS_MONTHS),
                 
                 OS_STATUS = c(initial_met_liver_os$OS_STATUS,
                               initial_met_others_os$OS_STATUS),
                 
                 age_death = c(initial_met_liver_os$AGE_AT_DEATH,
                               initial_met_others_os$AGE_AT_DEATH),
                 
                 initial_met = rep(x = c('liver',
                                         'other'),
                                   times = c(nrow(initial_met_liver),
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

pdf(file = 'a.pdf', width = 7.5, height = 9)
  ggsurvplot(km_fit, data = dt_,
             pval = T, pval.method = T, conf.int = T,
             risk.table = T,
             palette = c('red','blue'),
             xlab = "Time in Months", 
             ylab = "Overall Survival Probability", 
             title = "Kaplan-Meier Curve of Overall Survival")
graphics.off()

