# Most updated; Milo constructs a count matrix as input to edgeR. In this matrix column is patient sample, row is neighborhood and values are the number of cells in each neighborhood.
# Neighborhood construction: first it uses KNN to find the neighbors of each cell. Then it subsamples cells, as neighborhood centroids, to construct each neighborhood.

library(ggplot2)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v4")
library(statmod)


# Setting up ####

set.seed(42)
in_dir_ = '../data/integration/healthy/stroma/CAF/'
out_dir_ = '../out/integration/healthy/stroma/CAF/conditions/'
fl_nm = 'Milo_DA_k30_d30_prop0.4_CAF'
pdf(file = paste0(out_dir_,'/',fl_nm,'_stats.pdf'), width = 12, height = 7)
k_ = 30     # kth nearest neighbor to construct KNN
d_ = 30     # PCA dimensions to use

# Step 0: Reading in data ####

load(file = list.files(path = in_dir_, pattern = '_processed.RData', full.names = T))
DefaultAssay(s_objs) = 'integrated'
s_objs = subset(s_objs, subset = cell_type %in% c("activated stellate","CXCL14 fibroblast","inflammatory fibroblast","myofibroblast","proaxonogenic fibroblast"))

# initial status checks

plot(DimPlot(s_objs,
             pt.size = .5, shuffle = T,label = F, label.box = T,
             group.by = 'condition',
             raster = F))

plot(DimPlot(s_objs,
             pt.size = .5, shuffle = F,label = F, label.box = T,
             group.by = 'treatment',
             raster = F))

if(!file.exists(paste0(in_dir_,'/',fl_nm,'.RData')))
{
  # Step 1: Create Milo object from Seurat ####
  
  s_objs = as.SingleCellExperiment(s_objs)
  milo_obj = Milo(s_objs); rm(s_objs); gc()
  
  #Step 2: Constructing a KNN graph ####
  
  milo_obj = buildGraph(milo_obj,
                        k = k_,                   # number of nearest-neighbors (they used 30); 3
                        d = d_,                   # PCA dimensions (they used 30); 50
                        reduced.dim = "PCA")      # it is built from batch-corrected PCA
  
  # Step 3: Finding cell neighborhood size ####
  
  # cell neighborhood size: number of cells connected immediately to a given cell; For efficiency, they don’t use all neighborhoods
  # to cmpute spatial DA, but they sample a subset of cells as representative cells
  milo_obj = makeNhoods(milo_obj,
                        reduced_dims = "PCA",
                        prop = 0.4,         # proportion of cells to start with for sampling (usually 0.1 - 0.2 is sufficient
                        k = k_,             # the k to use for KNN refinement (they recommend using the same k used for KNN graph building)
                        d = d_,             # the number of reduced dimensions to use for KNN refinement (they recommend using the same d used for KNN graph building)
                        refined = TRUE)     # indicates whether you want to use the sampling refinement algorithm, or just pick cells at random. The default
  # and recommended way to go is to use refinement. The only situation in which you might consider using random instead, is if you
  # have batch corrected your data with a graph-based correction algorithm, such as BBKNN, but the results of DA testing will be suboptimal
  # milo_obj@nhoods: is a subset of KNN mtrix with #cells rows by (~prop* #cells) columns; prop is the starting point and could be changed by the algorithm
  
  plotNhoodSizeHist(milo_obj)               # it’s good to take a look at neighborhood sizes (vertex degrees). This affects the power of DA testing. Empirically, they found
  # it’s best to have a distribution peaking between 50 and 100. Otherwise we might consider rerunning makeNhoods increasing k and/or prop.
  
  # Step 4: Computing neighborhood diversity ####
  
  # Milo leverages the variation in the number of cells between replicates (replicate: a set of single cells like a patient sample) of a given condition (liver vs lung mets) to test for differential abundance.
  # They count the number of cells from each replicate (aka sample) in each neighborhood.They use metadata to specify which column contains the sample (replicate) information.
  milo_obj = countCells(milo_obj, meta.data = data.frame(colData(milo_obj)), sample = "orig.ident")     # This adds a m*n matrix, where rows are neighborhood IDs and columns samples (replicates).
  # Values indicate the number of cells from each sample in a neighborhood. This count matrix will be used for DA testing.
  head(as.matrix(nhoodCounts(milo_obj)))                                                                # This matrix, like expression matrix, could be given to edgR/DESeq2 where neighborhoods play genes differentially
  # occupied (differential abundance, DA) between conditions
  
  # Step 5: Defining experimental design ####
  
  # as in DESeq2/edgeR we need to specify a design;
  dsgn_ = data.frame(sample = colData(milo_obj)$"orig.ident",       # sample (replicate)
                     condition = colData(milo_obj)$"condition",     # condition
                     treatment = colData(milo_obj)$"treatment",
                     batch = colData(milo_obj)$"orig.ident",        # The batch effect is different than the one corrected for in Seurat in which we want to rule out systematic differences in the number of cells between two conditions because of library prep.
                     stringsAsFactors = T)                          # Batch here, however, refers to the groups of samples prepared and/or sequenced together; we don't have this info for PDAC; see Milo mouse example
  dsgn_ = distinct(dsgn_)                                           # removing duplicated records
  rownames(dsgn_) = dsgn_$sample
  print(dsgn_)
  
  # Step 6: Computing neighborhood connectivity ####
  
  # Milo uses an adaptation of the Spatial FDR correction, which accounts for the overlap between neighborhoods.
  # Specifically, each P-value is weighted by the reciprocal of the kth nearest neighbor (farthest neighbor) distance (the father the kth neighbor,
  # the smaller the p-value becomes after division). To use this statistic we first need to store the distances between nearest neighbors in the Milo object (time consuming)
  milo_obj =  calcNhoodDistance(milo_obj, d = d_, reduced.dim = "PCA")
  
  save(milo_obj, file = paste0(in_dir_,'/',fl_nm,'.RData'))
}else
{
  load(file = paste0(in_dir_,'/',fl_nm,'.RData'))
  
  plotNhoodSizeHist(milo_obj)
  
  dsgn_ = data.frame(sample = colData(milo_obj)$"orig.ident",       # sample (replicate)
                     condition = colData(milo_obj)$"condition",     # condition
                     treatment = colData(milo_obj)$"treatment",
                     batch = colData(milo_obj)$"orig.ident",        # this batch effect is different the one corrected for in Seurat; here we want to rule out systematic differences in the number of cells between two conditions because of library prep
                     stringsAsFactors = T)                          # batch here refers to groups of samples prepared and/or sequenced together; we don't have this info for PDAC; see Milo mouse example
  dsgn_ = distinct(dsgn_)                                           # removing duplicated records
  rownames(dsgn_) = dsgn_$sample
  print(dsgn_)
}


p_ = q_ = list()
i_ = j_ = 1
for(ds_ in c('condition','treatment'))
{
  # Step 7: DA (differential abundance) test ####
  
  # contrast
  fml_ = as.formula(paste0('~',ds_))
  da_results = testNhoods(milo_obj,
                          design.df = dsgn_,          # design data frame
                          design = fml_)              # design formula; batch effect produces error here; contrast: level 1 (ref) - level 2 (case)
  da_results$nh_size = colSums(nhoods(milo_obj))      # I added this line to add nhood sizes to the da_results; da_result is not sorted by LFc/p-value
  
  View(da_results %>% arrange(SpatialFDR))      # Nhood, like gene ID, is the neighborhood ID
  
  # Step 8: Inspection ####
  
  p_[[i_]] = ggplot(da_results, aes(PValue))+                          # if this shows a uniform distribution, it means no bias was identified; if it shows exponential distribution, it means there's a bias
    labs(title = paste0('Distribution DA p-values (',ds_,')'))+
    geom_histogram(bins = 100)
  
  p_[[i_+1]] = ggplot(da_results, aes(logFC, -log10(SpatialFDR)))+     # shows neighborhoods' adjusted p-values, they must be above 1 to be less than 0.1; FDR is not really useful as they manipulate them
    labs(title = paste0('Adjusted p-values of DA neighborhoods (',ds_,')'))+
    geom_point()+
    geom_hline(yintercept = 1)
  
  # Step 9: Plotting ####
  
  # step 9.1: plotting cell type bias
  
  da_results = annotateNhoods(milo_obj, da_results, coldata_col = "cell_type")                                  # it labels a neighborhood by the most abundant cell type in it (cell_type_fraction)
  head(da_results)
  
  p_[[i_+2]] = ggplot(da_results, aes(cell_type_fraction))+                                                     # Distribution of the fraction of the most abundant cell type in each neighborhood
               labs(title = paste0('Distribution of cell type fractions across all neighborhoods (',ds_,')'))+
               geom_histogram(bins = 50)
  
  da_results$cell_type = ifelse(da_results$cell_type_fraction < .65, "Mixed", da_results$cell_type)            # if the fraction of the most abundant cell type in a given neighborhood is less than cutoff,
  da_results$cell_type = as.factor(da_results$cell_type)                                                       # that neighborhood is not homogeneous
  
  # plotting by me based on L2FC
  
  # setting colors
  max_l2fc = max(abs(da_results$logFC))
  
  pos_cols = colorRampPalette(colors = c('skyblue','blue4'))(10)
  pos_ = seq(from = 0, to = max_l2fc, length.out = 10)     # each value if pos_ is the lower bound of each color interval
  names(pos_cols) = pos_
  
  neg_cols = colorRampPalette(colors = c('red4','pink'))(10)
  neg_ = seq(from = -max_l2fc, to = 0, length.out = 10)     # each value if neg_ is the upper bound of each color interval
  names(neg_cols) = neg_
  
  da_results$cols = NA
  for(r_ in 1:nrow(da_results))
  {
    l2fc_ = da_results$logFC[r_]
    if(0 <= l2fc_)
    {
      da_results$cols[r_] = pos_cols[as.character(max(pos_[pos_ <= l2fc_ ]))]
    }else
    {
      da_results$cols[r_] = neg_cols[as.character(min(neg_[ l2fc_ <= neg_ ]))]
    }
  }
  da_results$cols[0.1 < da_results$PValue] = 'grey88'
  da_results = da_results[order(da_results$cols, decreasing = T),]
  
  # plotting
  p_[[i_+3]] =
    ggplot(data = da_results, aes(x = logFC, y = cell_type))+
    theme(plot.title = element_text(hjust = .5, face = 'bold', family = 'Helvetica', size = 20),
          plot.background = element_blank(),panel.background = element_blank(),
          panel.grid = element_blank(), panel.grid.major.y = element_line(linewidth = .2, color = 'grey80'),
          panel.border = element_blank(),axis.line.x = element_line(color = 'black', linewidth = .5),
          axis.title = element_text(size = 25, hjust = .5, face = 'bold', family = 'Helvetica'),
          axis.text = element_text(face = 'bold', size = 20, family = 'Helvetica', color = 'black'),
          legend.text = element_text(family = 'Helvetica', size = 20, face = 'bold'),
          legend.title = element_text(family = 'Helvetica', size = 20, face = 'bold'),
          axis.line.y = element_line(color = 'black', linewidth = .5))+
    labs(title =  paste0('Cell type DA (',ds_,')'), x = 'L2FC', y = NULL)+
    geom_jitter(height = .25, aes(size = nh_size), color = da_results$cols)+
    geom_violin(scale = 'width', alpha = 0)+
    geom_vline(xintercept = c(-log2(1.5),log2(1.5)), linewidth = .5, linetype = 'dashed')
  
  # original by Milo based on spatialFDR

  # p_[[i_+3]] =
  #   plotDAbeeswarm(da_results, group.by = "cell_type",                                               # distribution of DA Fold Changes in different cell types
  #                  alpha = min(da_results$SpatialFDR)+1e-10)+                                                                      # alpha: !!!! very important: significance level for Spatial FDR (default: 0.1 but works [0.5-1])
  #   theme(plot.title = element_text(hjust = .5, face = 'bold', family = 'Helvetica'),
  #         axis.title = element_text(size = 25, hjust = .5, face = 'bold', family = 'Helvetica'),
  #         axis.text = element_text(face = 'bold', size = 20, family = 'Helvetica'),
  #         panel.grid = element_blank(), panel.grid.major.y = element_line(linewidth = .2, color = 'grey80'),
  #         panel.border = element_blank(),axis.line.x = element_line(color = 'black', linewidth = .5), axis.line.y = element_line(color = 'black', linewidth = .5))+
  #   labs(title =  paste0('Cell type DA (',ds_,')'), subtitle = round(min(da_results$SpatialFDR),2), x = NULL)+
  #   scale_color_gradient2(low = 'red', high = 'blue', mid = 'grey60', midpoint = 0,
  #                         breaks = c(min(da_results$logFC),-log2(1.5),0,log2(1.5), max(da_results$logFC)) )
  # p_[[i_+3]]$layers[[1]]$aes_params$size = 3
  
  
  # step 9.2: plotting condition bias
  
  milo_obj = buildNhoodGraph(milo_obj)      # it builds a graph of neighborhoods with edges representing number of shared neighbors
  
  # Plotting single-cell UMAP
  umap_pl = plotReducedDim(milo_obj, dimred = "UMAP", colour_by= ds_, text_by = "cell_type", text_size = 4, rasterise = F, point_alpha = 1, point_size = .4)+
            guides(color = guide_legend(title = ds_,
            title.theme = element_text(face = 'bold', size = 10),
            label.theme = element_text(face = 'bold', size = 10),
            override.aes = list(size = 5)))
  
  # Plotting neighborhood graph
  
  t_ = if(ds_ == 'condition'){ t_ = paste(levels(dsgn_$condition), collapse = ' vs ') }else{ t_ = paste(levels(dsgn_$treatment), collapse = ' vs ') }
  nh_graph_pl = plotNhoodGraphDA(x = milo_obj, milo_res = da_results,
                                 alpha = 1,
                                 size_range = c(0.4, 5))+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(title = t_)+
    guides(size = guide_legend("Neighborhood size"))
  nh_graph_pl$layers[[1]]$aes_params$edge_width = 0
  
  q_[[j_]] = umap_pl + nh_graph_pl + plot_layout(guides = "collect")
  
  i_ = i_ + 4; j_ = j_ + 1
}

null_ = lapply(X = p_, FUN = function(x_)
{
  plot(x_)
  return(NULL)
})
graphics.off()

for(i_ in 1:length(q_))
{
  ds_ = if(i_ == 1){ 'condition' }else{ 'treatment' }
  # png(filename = paste0(out_dir_,'/',fl_nm,'_',ds_,'_UMAP.png'), width = 14, height = 7, units = 'in', res = 300)
  pdf(file = paste0(out_dir_,'/',fl_nm,'_',ds_,'_UMAP.pdf'), width = 14, height = 7)
  plot(q_[[i_]])
  graphics.off()
}
