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

# converting non-character variables to factors

fcohort_$neoadjuvant = factor(fcohort_$neoadjuvant)                       # 0 (no) ref.
fcohort_$stage = factor(fcohort_$stage)                                   # 2 ref.
fcohort_$sex = factor(fcohort_$sex)                                       # F ref.
fcohort_$grade = factor(fcohort_$grade)                                   # 1 ref.
fcohort_$nodal_status = factor(fcohort_$nodal_status)                     # 0 (no) ref.
fcohort_$adjuvant = factor(fcohort_$adjuvant)                             # 0 (no) ref.

# Survival analysis by multivaiate Cox Proportional Hazards Regression Model ####

# cohort dataset
dt_ = fcohort_[fcohort_$recur_site %in% c("liver","lung"),]

# create a survival object
s_obj = Surv(dt_$yrs_to_recurrance_or_death, dt_$recur_status)

# adjusted model
fit_adj = coxph(s_obj ~ recur_site + neoadjuvant + age, data = dt_)
s_adj = summary(fit_adj)

# Forest plot ####

## prepare data for forestploter ####

df_ = data.frame(Subgroup = c('Recurrence',
                                'liver',
                                'lung',
                              
                              'Neoadjuvant',
                                'no',
                                'yes',
                              
                              'Age'),
                 
                             Number = c(# Recurrence
                                        '',
                                          nrow(dt_[dt_$recur_site %in% 'liver',]),# liver ref.
                                          nrow(dt_[dt_$recur_site %in% 'lung',]),# lung
                                        
                                        # Neoadjuvant
                                        '',
                                          nrow(dt_[dt_$neoadjuvant %in% 0,]),# no ref.
                                          nrow(dt_[dt_$neoadjuvant %in% 1,]),# yes
                                        
                                        # Age
                                        nrow(dt_[!is.na(dt_$age),])),
                             
                             est = c(# Recurrence
                                      NA,
                                        1,                          # liver ref.
                                        s_adj$conf.int[1],          #lung
                                      
                                      # Neoadjuvant
                                      NA,
                                        1,                          # no ref,
                                        s_adj$conf.int[2],          # yes
                                      
                                      # Age
                                      s_adj$conf.int[3]),
                             
                             low = c(# Recurrence
                                     NA,
                                       1,                          # liver ref.
                                       s_adj$conf.int[7],          #lung
                                     
                                     # Neoadjuvant
                                     NA,
                                       1,                          # no ref,
                                       s_adj$conf.int[8],          # yes
            
                                     # Age
                                     s_adj$conf.int[9]),
                             
                             hi = c(# Recurrence
                                    NA,
                                      1,                          # liver ref.
                                      s_adj$conf.int[10],          #lung
                                     
                                     # Neoadjuvant
                                    NA,
                                      1,                          # no ref,
                                      s_adj$conf.int[11],          # yes
                                     
                                    # Age
                                    s_adj$conf.int[12]),
                             p_value = c(# Recurrence
                                         NA,
                                            NA,                          # liver ref.
                                            s_adj$coefficients[13],          #lung
                                       
                                        # Neoadjuvant
                                         NA,
                                            NA,                          # no ref,
                                         s_adj$coefficients[14],          # yes
                                       
                                        # Age
                                        s_adj$coefficients[15]),
                             FACTOR = c(T, F, F,
                                        T, F, F,
                                        T))

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

ggsave(filename = "forestplot_DFS_cohort.pdf", plot = p_, width = 7.5, height = 7.5)
graphics.off()
