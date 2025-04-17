library(survival)
library(ggplot2)
library(grid)
library(gridExtra)
library(survminer)
library(swimplot)
library(car)
library(forestploter)

# Reading cohort data in ####

fcohort_ = read.delim(file = "Misc/fullcohort_cleaned.tsv", header = T, sep = '\t', quote = "", as.is = T, check.names = F, stringsAsFactors = F, na.strings = c('','NA','n/a','na'))
fcohort_$vital_status[!fcohort_$vital_status %in% 1] = 0      # event of interest is dead:1, the rest 0

fcohort_$diagnosis_dt = as.Date(fcohort_$diagnosis_dt)

fcohort_$last_follow_up = as.Date(fcohort_$last_follow_up)
fcohort_$yrs_to_last_followup = (fcohort_$last_follow_up - fcohort_$diagnosis_dt)/365.25

fcohort_$recurrence_dt = as.Date(fcohort_$recurrence_dt)
fcohort_$yrs_to_recurrance = (fcohort_$recurrence_dt - fcohort_$diagnosis_dt)/365.25
fcohort_$yrs_to_recurrance[is.na(fcohort_$yrs_to_recurrance)] = 100000

fcohort_$yrs_to_recurrance_or_death = pmin(fcohort_$yrs_to_recurrance,fcohort_$yrs_to_last_followup)

fcohort_$recur_status = pmax(fcohort_$recur, fcohort_$vital_status)

# Survival analysis by multivaiate Cox Proportional Hazards Regression Model ####

# cohort subset
dt_ = rbind(fcohort_[fcohort_$recur_site %in% "lung",], fcohort_[fcohort_$recur_site %in% "liver",])

# create a survival object
s_obj = Surv(dt_$yrs_to_last_followup, dt_$vital_status)

# adjusted model

fit_adj = coxph(s_obj ~ factor(recur_site) + factor(neoadjuvant) + factor(grade) + age +
                        factor(adjuvant) + factor(sex) + factor(stage) + factor(nodal_status), data = dt_)
s_adj = summary(fit_adj)

# Forest plot ####

## prepare data for forestploter ####

df_ = data.frame(Subgroup = c('Recurrence',
                                'liver',
                                'lung',
                              
                              'Neoadjuvant',
                                'no',
                                'yes',
                              
                              'Grade',
                                '1',
                                '2',
                                '3',
                              
                              'Age',
                              
                              'Adjuvant',
                                'no',
                                'yes',
                              
                              'Sex',
                                'F',
                                'M',
                              
                              'Stage',
                                '2',
                                '3',
                              
                              'Nodal status',
                                'no',
                                'yes'),
                 
                 Number = c(# Recurrence
                            '',
                              nrow(dt_[dt_$recur_site %in% 'liver',]),# liver ref.
                              nrow(dt_[dt_$recur_site %in% 'lung',]),# lung
                            
                            # Neoadjuvant
                            '',
                              nrow(dt_[dt_$neoadjuvant %in% 0,]),# no (0) ref.
                              nrow(dt_[dt_$neoadjuvant %in% 1,]),# yes (1)
                            
                            # Grade
                            '',
                              nrow(dt_[dt_$grade %in% 1,]),# 1 ref.
                              nrow(dt_[dt_$grade %in% 2,]),# 2
                              nrow(dt_[dt_$grade %in% 3,]),# 3
                            
                            # Age
                            nrow(dt_[!is.na(dt_$age),]),
                            
                            # Adjuvant
                            '',
                              nrow(dt_[dt_$adjuvant %in% 0,]),# no (0) ref.
                              nrow(dt_[dt_$adjuvant %in% 1,]),# yes (1)
                            
                            # Sex
                            '',
                              nrow(dt_[dt_$sex %in% 'F',]),# F ref.
                              nrow(dt_[dt_$sex %in% 'M',]),# M
                            
                            # Stage
                            '',
                              nrow(dt_[dt_$stage %in% 2,]),# 2 ref.
                              nrow(dt_[dt_$stage %in% 3,]),# 3
                            
                            # Nodal status
                            '',
                              nrow(dt_[dt_$nodal_status %in% 0,]),# no (0) ref.
                              nrow(dt_[dt_$nodal_status %in% 1,])),# yes (1)
                 
                 est = c(# Recurrence
                         NA,
                         1,                          # liver ref.
                         s_adj$conf.int[1],          # lung
                          
                         # Neoadjuvant
                         NA,
                         1,                          # 0 ref,
                         s_adj$conf.int[2],          # 1
                          
                         # Grade
                         NA,
                         1,                          # 1 ref
                         s_adj$conf.int[3],          # 2
                         s_adj$conf.int[4],          # 3
                          
                         # Age
                         s_adj$conf.int[5],
                          
                         # Adjuvant
                         NA,
                          1,                          # 0 ref
                          s_adj$conf.int[6],          # 1
                          
                         # Sex
                         NA,
                          1,                          # F ref
                          s_adj$conf.int[7],          # M
                          
                         # Stage
                         NA,
                          1,                          # 2 ref
                          s_adj$conf.int[8],          # 3
                          
                        # Nodal status
                        NA,
                         1,                          # 0 ref
                         s_adj$conf.int[9]),         # 1,
                 
                 low = c(# Recurrence
                         NA,
                           1,                          # liver ref.
                           s_adj$conf.int[19],         # lung
                         
                         # Neoadjuvant
                         NA,
                           1,                          # 0 ref,
                           s_adj$conf.int[20],         # 1
                         
                         # Grade
                         NA,
                           1,                          # 1 ref
                           s_adj$conf.int[21],         # 2
                           s_adj$conf.int[22],         # 3
                         
                         # Age
                         s_adj$conf.int[23],
                         
                         # Adjuvant
                         NA,
                          1,                          # 0 ref
                          s_adj$conf.int[24],         # 1
                         
                         # Sex
                         NA,
                          1,                          # F ref
                          s_adj$conf.int[25],         # M
                         
                         # Stage
                         NA,
                          1,                          # 2 ref
                          s_adj$conf.int[26],         # 3
                         
                         # Nodal status
                         NA,
                          1,                          # 0 ref
                          s_adj$conf.int[27]),        # 1,
                 
                 hi = c(# Recurrence
                        NA,
                          1,                          # liver ref.
                          s_adj$conf.int[28],         # lung
                         
                         # Neoadjuvant
                        NA,
                          1,                          # 0 ref,
                          s_adj$conf.int[29],         # 1
                         
                         # Grade
                        NA,
                          1,                          # 1 ref
                          s_adj$conf.int[30],         # 2
                          s_adj$conf.int[31],         # 3
                         
                        # Age
                        s_adj$conf.int[32],
                        
                        # Adjuvant
                        NA,
                         1,                          # 0 ref
                         s_adj$conf.int[33],         # 1
                        
                        # Sex
                        NA,
                         1,                          # F ref
                         s_adj$conf.int[34],         # M
                        
                        # Stage
                        NA,
                         1,                          # 2 ref
                         s_adj$conf.int[35],         # 3
                        
                        # Nodal status
                        NA,
                         1,                          # 0 ref
                         s_adj$conf.int[36]),        # 1,
                 
                 p_value = c(# Recurrence
                             NA,
                                NA,                          # liver ref.
                                s_adj$coefficients[37],      # lung
                           
                            # Neoadjuvant
                             NA,
                                NA,                          # 0 ref,
                                s_adj$coefficients[38],      # 1
                           
                            # Grade
                            NA,
                              NA,                          # 1 ref
                              s_adj$coefficients[39],      # 2
                              s_adj$coefficients[40],      # 3
                           
                            # Age
                            s_adj$coefficients[41],
                            
                            # Adjuvant
                            NA,
                              NA,                          # 0 ref
                              s_adj$coefficients[42],      # 1
                            
                            # Sex
                            NA,
                              NA,                          # F ref
                              s_adj$coefficients[43],      # M
                            
                            # Stage
                            NA,
                              NA,                          # 2 ref
                              s_adj$coefficients[44],      # 3
                            
                            # Nodal status
                            NA,
                              NA,                          # 0 ref
                              s_adj$coefficients[45]),     # 1,
                 
                 FACTOR = c(T, F, F,      # recurrence
                            T, F, F,
                            T, F, F, F,
                            T,            # age
                            T, F, F,
                            T, F, F,
                            T, F, F,
                            T, F, F))

# Indent the subgroup if there is a number in the placebo column
df_$Subgroup = ifelse(df_$FACTOR, yes = df_$Subgroup, no = paste0("   ", df_$Subgroup))

df_$p_value[!is.na(df_$p_value)] = sprintf("%.4f", df_$p_value[!is.na(df_$p_value)])
df_$p_value[is.na(df_$p_value)] = ''

# Create a confidence interval column to display
df_$`HR (95% CI)` = ifelse(df_$p_value %in% '', yes = "", no = sprintf("%.2f (%.2f to %.2f)", df_$est, df_$low, df_$hi))
df_$`HR (95% CI)`[df_$est %in% 1] = 'Reference'

# adding a blank column for the forest plot to display CI.

df_$` ` <- paste(rep(" ", 20), collapse = " ")      # increase the number of spaces below to provide a larger area for drawing the CI

## forest plot ####

p_ = forest(data = df_[, c("Subgroup","Number","HR (95% CI)", " ", "p_value")],
            est = df_$est,
            lower = df_$low,
            upper = df_$hi,
            sizes = 1,
            ci_column = 4,
            ref_line = 1,
            arrow_lab = c("Decreased risk", "Increased risk"),
            xlim = c(0, 3),
            ticks_at = c(0, 0.50, 1, 2, 3),
            footnote = sprintf('\n\n\nConcordance: %.2f,\nLog-rank p: %.4f,\nNo. events: %i',
                               s_adj$concordance[1],
                               s_adj$sctest[3],
                               s_adj$nevent),
            theme = forest_theme(ci_Theight = 0.2,     # T-end CI's
                                 footnote_gp = gpar(cex = 1, family = 'Helvetica', fontface = "bold", col = "black")))

ggsave(filename = "OS_cohort_forestplot_all_variables.pdf", plot = p_, width = 7.5, height = 7.5)
graphics.off()
