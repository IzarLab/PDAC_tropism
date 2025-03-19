# This scripts plots Fig. 3b; it carries out differential cell abundance analysis by Milo on malignant exocrine compartment.

library(ggplot2)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(statmod)

# Loading data ####

load(file = 'Misc/Milo_DA_k30_d30_prop0.4_cancer.RData')

# Designing the test ####

dsgn_ = data.frame(sample = colData(milo_obj)$"orig.ident",       # sample (replicate)
                   condition = colData(milo_obj)$"condition",     # condition
                   batch = colData(milo_obj)$"orig.ident",
                   stringsAsFactors = T)
dsgn_ = distinct(dsgn_)                                           # removing duplicated records
rownames(dsgn_) = dsgn_$sample
print(dsgn_)

# DA (differential abundance) test ####

# test

ds_ = 'condition'
fml_ = as.formula(paste0('~',ds_))
da_results = testNhoods(milo_obj,
                        design.df = dsgn_,                                                            # design data frame
                        design = fml_)
da_results$nh_size = colSums(nhoods(milo_obj))                                                        # I added this line to add nhood sizes to the da_results; da_result is not sorted by LFc/p-value

# cell type assignment

da_results = annotateNhoods(milo_obj, da_results, coldata_col = "cell_type")                          # it labels a neighborhood by the most abundant cell type in it (cell_type_fraction)
head(da_results)

da_results$cell_type = ifelse(da_results$cell_type_fraction < .65, "Mixed", da_results$cell_type)     # if the fraction of the most abundant cell type in a given neighborhood is less than cutoff,
da_results$cell_type = as.factor(da_results$cell_type)                                                # that neighborhood is not homogeneous

# Plotting ####

pdf(file = 'b.pdf', width = 10, height = 10)

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
da_results$cols[0.10 < da_results$PValue] = 'grey88'
da_results = da_results[order(da_results$cols, decreasing = T),]

# plotting violins

p_ =  ggplot(data = da_results, aes(x = logFC, y = cell_type))+
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
plot(p_)

# plotting neighborhood graph

milo_obj = buildNhoodGraph(milo_obj)      # it builds a graph of neighborhoods with edges representing number of shared neighbors

t_ = paste(levels(dsgn_$condition), collapse = ' vs ')
nh_graph_pl = plotNhoodGraphDA(x = milo_obj, milo_res = da_results,
                               alpha = 1,
                               size_range = c(0.4, 5))+
              theme(plot.title = element_text(hjust = 0.5))+
              labs(title = t_)+
              guides(size = guide_legend("Neighborhood size"))
nh_graph_pl$layers[[1]]$aes_params$edge_width = 0

# Plotting single-cell UMAP

umap_pl = plotReducedDim(milo_obj, dimred = "UMAP", colour_by= ds_, text_by = "cell_type", text_size = 4, rasterise = F, point_alpha = 1, point_size = .4)+
          guides(color = guide_legend(title = ds_,
                                      title.theme = element_text(face = 'bold', size = 10),
                                      label.theme = element_text(face = 'bold', size = 10),
                                      override.aes = list(size = 5)))

q_ = umap_pl + nh_graph_pl + plot_layout(guides = "collect")
plot(q_)

graphics.off()
