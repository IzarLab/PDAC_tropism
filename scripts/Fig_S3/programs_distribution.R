# Purpose: generates stacked bar plots per patient for each malignant program, comparing liver vs lung.
# Algorithm:
# 1 for each cell, it calculates program mean signature expression level by finding the mean of normalized UMI counts across all the genes in a program
# 2 for each cell, it calculates background mean expression level by randomly selecting the same number of genes as in a program and finding their mean, iterating this 5e4 times to generate empirical p-values
# 3 a cell is considered significant if its empirical p-value is less than 1e-3
# 4 for each significant cell, it calculates cell-program-signature score by subtracting the background mean from the signature mean
# 5 for each patient, it calculates patient-program-signature score by finding the mean of all cell-program-signature scores
# 6 it plots bar plots
# 7 it finds programs significantly different between liver and long using Wilcoxon rank sum test

library(ggplot2)
library(viridis)
library(gridExtra)
library(Seurat)
library(parallel)

set.seed(41)

if(!file.exists('malignant_programs'))
{
  # reading in program signatures ####
  
  sgns_ = list()
  sgns_[['Basal-like']] = c('SDC4','ITGA3','MYOF','KRT7','EMP1','TYMP','SCEL','SLC25A37','PLAU','PLEK2','LAMC2','DCBLD2','PODXL','LAMB3','NAV2','CDK14','SEMA7A','PHLDB2','GRB10','F3','CASC15','PMEPA1','MUC16','SERPINE1','HMGA2','BCAR3','ETS1','KRT17','TNFAIP2','ITGA2')
  sgns_[['Classical']] = c('MUC5AC','FER1L6','CEACAM5','CEACAM6','MUCL3','MUC17','MUC3A','TFF1','REG4','SULT1C2','S100A6','FTH1','B2M','TFF2','DHRS9','SYTL2','GPC5','ACTB','BCAS1','AGR2','CDH17','CD55','BTNL8','CLDN18','GCNT3','PSCA','HMGCS1','PLAC8','MUC13','ACTG1')
  sgns_[['Progenitor']] = c("HNF1B", "SOX6", "REG1A", "FGFR2", 'MUC6','OLFM4','CFTR','SOX5','CHRM3','SLC4A4','CWC27','PKHD1','PGC','LRIG3','SCTR','PIK3C2G','PAK3','CFAP221','LGR5','PLEKHS1','SLC12A2','PTPRG','RORA','NFIA','MEIS2','DYNC2H1','KCNJ15','MSR1','ARHGAP24','HOMER2')
  sgns_[['Ciliated']] = c('DNAH12','DNAH3','DNAH7','DNAH6','DNAH11','CFAP54','CFAP69','CFAP43','CFAP47','CFAP77','CFAP157','CFAP70','CFAP52','CFAP44','CFAP46','CFAP206')
  for(f_ in list.files('../data/datasets/pathways/', full.names = T))
  {
    signal_ = read.table(file = f_, sep = '\t', header = F, quote = '', row.names = 1, as.is = T, check.names = F)
    sgns_[[ signal_['STANDARD_NAME',"V2"] ]] = strsplit(signal_['MAPPED_SYMBOLS','V2'], split = ',')[[1]]
  }
  
  # reading in cancer data ####
  
  load(file = '../data/integration/PN30/other/cancer/treated_untreated/integrated_untreated_treated_cancer_data_QC_CST_processed.RData')
  DefaultAssay(s_objs) = 'RNA'
  conditions_ = s_objs$condition
  names(conditions_) = s_objs$orig.ident
  
  # Main part ####
  
  probs_pnts_sgns = bgrnd_means_pnts_sgns = patient_signature_scores = list()
  p_ = list(); i_ = 1
  for(sgn_ in 1:length(sgns_))
  {
    sgn_ = names(sgns_[sgn_])
    message(sgn_,'\n\t', appendLF = F)
    
    sgn_genes = sgns_[[sgn_]]
    rslt_ = mclapply(X = unique(s_objs$orig.ident), FUN = function(pnt_, genes_ = sgn_genes)
                                                          {
                                                            message(pnt_)
                                                            
                                                            cells_ = colnames(s_objs)[ s_objs$orig.ident %in% pnt_]     # cancer cells of patient pnt_
                                                            s_obj = s_objs[['RNA']]@data[,cells_]                       # expression matrix of pnt_
                                                            
                                                            sgn_size = length(genes_)                                   # size of program signature
                                                            genes_ = genes_[genes_ %in% rownames(s_obj)]                # signature markers; those expressed in pnt_
                                                            signature_means = colSums(s_obj[genes_,])/sgn_size          # signature mean of each cancer cell
                                                            
                                                            # finding significant cells expressing this program
                                                            
                                                            genes_bgrnd = rownames(s_obj)[rowSums(s_obj) > 0]           # background genes; those with values for all cells and
                                                            genes_bgrnd = genes_bgrnd[!genes_bgrnd %in% genes_]         # minus markers
                                                            bgrnd_means = list()                                        # background mean for each cell in the form of a matrix with rows containing background means for iterations
                                                            for(itr_ in 1:1e3)
                                                            {
                                                              genes_tmp = sample(genes_bgrnd, size = sgn_size)
                                                              bgrnd_means[[itr_]] = colMeans(s_obj[genes_tmp,])
                                                            }
                                                            bgrnd_means = do.call(bgrnd_means, what = rbind)
                                                            
                                                            probs_ = numeric(length = length(cells_))
                                                            names(probs_) = cells_
                                                            for(cell_ in cells_)
                                                            {
                                                              probs_[cell_] = round(1- ecdf(bgrnd_means[,cell_])(signature_means[cell_]), 5)
                                                            }
                                                            cells_ = names(probs_[probs_ <= 0.001])                     # retaining only significant cancer cells #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                            if(length(cells_) < 5){ cat('No cell was identified with this program'); return(NA) }#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                            
                                                            q_ = DimPlot(s_objs, cells.highlight = cells_)+labs(title = paste0(pnt_,': ',sgn_))
                                                            
                                                            cell_signature_scores = signature_means[cells_] - colMeans(bgrnd_means[, cells_, drop = F])
                                                            
                                                            return(list(patient_signature_score = mean(cell_signature_scores), bgrnd_means_pnts = bgrnd_means, probs_pnts = probs_, plot = q_))
                                                          }, sgn_genes, mc.cores = min(detectCores()-1,length(unique(s_objs$orig.ident))) )
    names(rslt_) = unique(s_objs$orig.ident)
    rslt_ = rslt_[! is.na(rslt_)]
    if(length(rslt_) == 0){ next() }
    
    patient_signature_score = bgrnd_means_pnts = probs_pnts = list()
    for(e_ in 1:length(rslt_))
    {
      pnt_ = names(rslt_[e_])
      e_ = rslt_[[e_]]
      patient_signature_score[[pnt_]] =  e_$patient_signature_score
      bgrnd_means_pnts[[pnt_]] =  e_$bgrnd_means_pnts
      probs_pnts[[pnt_]] = e_$probs_pnts
      p_[[i_]] = e_$plot
      i_ = i_ + 1
    }
    patient_signature_scores[[sgn_]] =  patient_signature_score
    bgrnd_means_pnts_sgns[[sgn_]] =  bgrnd_means_pnts
    probs_pnts_sgns[[sgn_]] = probs_pnts
  }
  save(patient_signature_scores, probs_pnts_sgns, bgrnd_means_pnts_sgns, p_, conditions_, file = 'malignant_programs.RData')
}else
{
  load(file = 'result.RData')
}

# Plotting scatter plots ####

pdf(file = paste0('programs_scatterplots.pdf'), width = 80, height = 7*(length(p_)/10))
  grid.arrange(arrangeGrob(grobs = p_, ncol = 10))
graphics.off()

# Plotting bar plots ####

dt_ = list()
for(i_ in 1:length(patient_signature_scores))
{
  sgn_ = unlist(patient_signature_scores[[i_]])
  dt_[[i_]] = data.frame(sample = names(sgn_), score = sgn_, program = names(patient_signature_scores[i_]), stringsAsFactors = T)
}
dt_ = do.call(dt_, what = rbind)
dt_$condition = as.factor(conditions_[as.character(dt_$sample)])

pdf(file = 'malignant_programs_distribution.pdf', width = 21)
lineage_ = c('Basal-like', 'Classical', 'Progenitor', 'Ciliated')
states_ = as.character(unique(dt_$program[!dt_$program %in% lineage_]))
cols_ = sample(x = viridis(n = length(states_), option = 'turbo'), size = length(states_))
names(cols_) = states_
for(i_ in 1:2)
{
  dt_tmp = lims_ = option_ = NULL
  if(i_ == 1){ dt_tmp = dt_[dt_$program %in% lineage_,]; lims_ = c(0, 3) }
  if(i_ == 2){ dt_tmp = dt_[!dt_$program %in% lineage_,]; lims_ = c(0, 3.5) }

  dt_liver = dt_tmp[dt_tmp$condition %in% 'liver',]
  p_ =  ggplot(data = dt_liver, aes(fill = program, y = score, x = sample, label = round(score,2)))+
        theme(panel.background = element_rect(fill = 'white'), panel.grid = element_line(color = 'grey90'), text = element_text(face = 'bold'))+
        labs(title = 'Liver', x = 'Patient', y = 'Program expression', fill = 'Malignant program')+
        geom_bar(position="stack", stat="identity", colour = "black")+
        geom_text(size = 3, position = position_stack(vjust = 0.5), color = 'white', fontface = 2)+
        scale_y_continuous(limits = lims_)
  if(i_ == 1){ p_ = p_ + scale_fill_viridis(discrete = T) }else{ p_ = p_ + scale_fill_manual(values = cols_)}

  dt_lung = dt_tmp[dt_tmp$condition %in% 'lung',]
  q_ =  ggplot(data = dt_lung, aes(fill = program, y = score, x = sample, label = round(score,2)))+
        theme(panel.background = element_rect(fill = 'white'), panel.grid = element_line(color = 'grey90'), text = element_text(face = 'bold'))+
        labs(title = 'Lung', x = 'Patient', y = '', fill = 'Malignant program')+
        geom_bar(position="stack", stat="identity", colour="black")+
        geom_text(size = 3, position = position_stack(vjust = 0.5), color = 'white', fontface = 2)+
        scale_y_continuous(limits = lims_)
  if(i_ == 1){ q_ = q_ + scale_fill_viridis(discrete = T) }else{ q_ = q_ + scale_fill_manual(values = cols_)}

  plot(p_ + q_)
}
graphics.off()


p_val = numeric()
for(sgn_ in unique(dt_$program))
{
  liv_scores = dt_$score[dt_$program %in% sgn_ & dt_$condition %in% 'liver']
  lng_scores = dt_$score[dt_$program %in% sgn_ & dt_$condition %in% 'lung']
  if(length(liv_scores) < 4 | length(lng_scores) < 4){ cat('not possible for ',sgn_,'\n'); next() }
  p_val[sgn_] = wilcox.test(liv_scores, lng_scores, alternative = 'two.sided')$p.value
}
print(sort(p_val))

# pdf(file = 'programs.pdf', width = 14, height = 7.5)
# lung_cells = colnames(s_objs)[s_objs$condition %in% 'lung']
# liver_cells = colnames(s_objs)[s_objs$condition %in% 'liver']
# for(sgn_ in names(sgns_))
# {
#   cutoff_lung = mean(range(s_objs@meta.data[lung_cells,sgn_]))
#   cutoff_liver = mean(range(s_objs@meta.data[liver_cells,sgn_]))
#   
#   lung_cells_pos = lung_cells[ s_objs@meta.data[lung_cells, sgn_] >= cutoff_lung ]
#   liver_cells_pos = liver_cells[s_objs@meta.data[liver_cells, sgn_] >= cutoff_liver ]
#   
#   lng_t = round(length(lung_cells_pos)/length(lung_cells)*100,2)
#   liv_t = round(length(liver_cells_pos)/length(liver_cells)*100,2)
#   
#   cutoff_ = min(cutoff_lung, cutoff_liver)
#   
#   p_ = FeaturePlot(s_objs, features = sgn_,
#                    cells = lung_cells,
#                    pt.size = .7, cols = c('grey90','red3'), min.cutoff = cutoff_lung-0.3*cutoff_lung, max.cutoff = 'q99', order = F)+
#     theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold'), axis.text = element_blank(), axis.ticks = element_blank())+
#     labs(title = paste0('Lung (',lng_t,'%)\n',sgn_))
#   q_ = FeaturePlot(s_objs, features = sgn_,
#                    cells = liver_cells,
#                    pt.size = .7, cols = c('grey90','red3'), min.cutoff = cutoff_liver-0.3*cutoff_liver, max.cutoff = 'q99', order = F)+
#     theme(legend.position = 'none', axis.title = element_blank(), text = element_text(face = 'bold'), axis.text = element_blank(), axis.ticks = element_blank())+
#     labs(title = paste0('Liver (',liv_t,'%)\n',sgn_))
#   w_ = plot_grid(q_, p_)
#   plot(w_)
# }
# graphics.off()
