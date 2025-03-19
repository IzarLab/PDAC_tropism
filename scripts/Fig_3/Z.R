# This scripts checks if DEGA result is affected by a few samples (treated, untreated, size, etc.) with large sizes using leave-one-out (LOO) approach.
# It does not corresponds to any figure in the paper by rebuttal for reviewers.

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratWrappers)
library(ggplot2)
library(gridExtra)
require(patchwork)
require(dbscan)
library(xlsx)
library(monocle3)
library(parallel)

# Reading in data ####
## reading in cancer compartment containing main PDAC subtypes ####

load('compartments/cancer_subtypes.RData')
DefaultAssay(s_objs) = 'RNA'

## reading in PDAC-recurrence signatures (LIV-MR and LUN-MR DE gene sets) ####

pdac_liver_markers = read.delim('Misc/cell_types_markers_main_subtypes_liver.tsv', header = T, sep = '\t')
pdac_liver_markers = pdac_liver_markers[pdac_liver_markers$cluster %in% 'all' & pdac_liver_markers$avg_log2FC > 0,]
pdac_liver_markers = trimws(pdac_liver_markers$gene)

pdac_lung_markers = read.delim('Misc/cell_types_markers_main_subtypes_lung.tsv', header = T, sep = '\t')
pdac_lung_markers = pdac_lung_markers[pdac_lung_markers$cluster %in% 'all' & pdac_lung_markers$avg_log2FC > 0,]
pdac_lung_markers = trimws(pdac_lung_markers$gene)


# Performing LOO ####

## fraction of each patient sample withing each group (LIV-MR vs LUN-MR) ####

lung_smpls = round(table(s_objs$orig.ident[s_objs$condition %in% 'lung'])/sum(table(s_objs$orig.ident[s_objs$condition %in% 'lung']))*100,2)
liver_smpls = round(table(s_objs$orig.ident[s_objs$condition %in% 'liver'])/sum(table(s_objs$orig.ident[s_objs$condition %in% 'liver']))*100,2)

## performing LOO ####

tbl_liv = tbl_lng = list()                  # final table holding the result for LIV0MR and LUN-MR patient groups
for(smpl_ in unique(s_objs$orig.ident))     # for each patient
{
  message('sample out: ',smpl_)
  
  s_ = subset(s_objs, subset = orig.ident != smpl_)       # cancer data except the current patient left out
  
  # performing DGEA
  markers_ = FindMarkers(object = s_,
                         group.by = 'condition', ident.1 = 'lung', ident.2 = 'liver')
  markers_ = markers_[order(abs(markers_$avg_log2FC), decreasing = T),]
  markers_$gene = rownames(markers_)
  
  liver_vs_lung = markers_[markers_$avg_log2FC < 0, ]     # LIV-MR markers recovered
  lung_vs_liver = markers_[0 < markers_$avg_log2FC,]      # LUN-MR markers recovered
  
  # fraction of DEGs recovered after leaving one patient out
  idxs_ = which(pdac_liver_markers %in% liver_vs_lung$gene)
  frac_liver = ceiling(length(idxs_)/length(pdac_liver_markers)*100)
  
  idxs_ = which(pdac_lung_markers %in% lung_vs_liver$gene)
  frac_lung = ceiling(length(idxs_)/length(pdac_lung_markers)*100)
  
  # updating result tables
  if(smpl_ %in% names(liver_smpls))
  {
    tbl_liv[[smpl_]] = data.frame(sample_out = smpl_, condition = 'liver',
                                  fraction_in_liver = liver_smpls[smpl_],
                                  liver_markers_covered = frac_liver,
                                  stringsAsFactors = T)
  }else
  {
    tbl_lng[[smpl_]] = data.frame(sample_out = smpl_, condition = 'lung',
                                  fraction_in_lung = lung_smpls[smpl_],
                                  lung_markers_covered = frac_lung,
                                  stringsAsFactors = T)
  }
}
tbl_liv = do.call(tbl_liv, what = rbind)
tbl_liv = tbl_liv[order(tbl_liv$fraction_in_liver, decreasing = T),]
write.table(x = tbl_liv, file = 'PDAC-liver_signature_leave-one-out.tsv', quote = F, sep = '\t', row.names = F, col.names = T)      # writing the result table for LIV-MR

tbl_lng = do.call(tbl_lng, what = rbind)
tbl_lng = tbl_lng[order(tbl_lng$fraction_in_lung, decreasing = T),]
write.table(x = tbl_lng, file = 'PDAC-lung_signature_leave-one-out.tsv', quote = F, sep = '\t', row.names = F, col.names = T)      # writing the result table for LUN-MR
