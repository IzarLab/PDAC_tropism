# This script plots Extended Data Fig. 6d; it plots KM curve for OS of MSK-MET data taken from: https://zenodo.org/records/5801902
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

# start and end time for DFS
for(i_ in 1:nrow(dt_pnt))
{
  dt_pnt$start[i_] = min(dt_pnt$AGE_AT_SEQUENCING[i_], dt_pnt$AGE_AT_SURGERY[i_], na.rm = T)      # sequencing could be from blood sample at the time of diagnosis for ctDNA or germline WGS 
  dt_pnt$end[i_] = min(dt_pnt$AGE_AT_EVIDENCE_OF_METS[i_], dt_pnt$AGE_AT_DEATH[i_], dt_pnt$AGE_AT_LAST_CONTACT[i_], na.rm = T)
  dt_pnt$DFS_MONTH[i_] = dt_pnt$end[i_] - dt_pnt$start[i_]
  dt_pnt$DFS_STATUS[i_] = if(!is.na(dt_pnt$AGE_AT_EVIDENCE_OF_METS[i_])){ 1 }else{ 0 }
}
dt_pnt = dt_pnt[!is.infinite(dt_pnt$start) & !is.infinite(dt_pnt$end) & (0 <= dt_pnt$DFS_MONTH),]

# initial metastatic site

dt_smpl = dt_smpl[dt_smpl$MET_SITE_COUNT == 1,]     # only one metastatic site during the follow-up time; among MET_COUNT and MET_SITE_COUNT the latter corresponds to number of Yes's
initial_met_lung = dt_smpl[dt_smpl$MET_SITE %in% 'DMETS_DX_LUNG',]
initial_met_lung_os = dt_pnt[dt_pnt$PATIENT_ID %in% initial_met_lung$PATIENT_ID,]

initial_met_others = dt_smpl[!dt_smpl$MET_SITE %in% c('DMETS_DX_LUNG','DMETS_DX_LIVER','DMETS_DX_DIST_LN'),]
initial_met_others_os = dt_pnt[dt_pnt$PATIENT_ID %in% initial_met_others$PATIENT_ID,]

# Plotting KM curve ####

# preparing data

dt_ = data.frame(patients   = c(initial_met_lung_os$PATIENT_ID,
                                initial_met_others_os$PATIENT_ID),
                 
                 DFS_MONTH  = c(initial_met_lung_os$DFS_MONTH,
                                initial_met_others_os$DFS_MONTH),
                 
                 DFS_STATUS = c(initial_met_lung_os$DFS_STATUS,      # it is required to use a numerical variable to be able to use ggsurvplot, as it does not work with factor
                                initial_met_others_os$DFS_STATUS),
                 
                 initial_met = rep(x = c('lung',
                                         'others'),
                                   times = c(nrow(initial_met_lung_os),
                                             nrow(initial_met_others_os))),
                 stringsAsFactors = T)

# create a Survival Object

surv_object = Surv(time = dt_$DFS_MONTH, event = dt_$DFS_STATUS)

# fit a Kaplan-Meier Model

km_fit = survfit(surv_object ~ initial_met, data = dt_)

# plotting

pdf(file = 'd.pdf', width = 7.5, height = 9)
ggsurvplot(km_fit, data = dt_,
           pval = T, pval.method = T, conf.int = T,
           risk.table = T,
           palette = rev(c('red','blue')),
           xlab = "Time in Months", 
           ylab = "DFS Probability", 
           title = "Kaplan-Meier Curve of DFS")
graphics.off()
