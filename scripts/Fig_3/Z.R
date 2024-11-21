# This scripts checks if DEGA result is not affected by a few samples with large sizes using one-leave-out approach

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


load('compartments/cancer_subtypes.RData')
DefaultAssay(s_objs) = 'RNA'

# these markers should be updated with those which overlap with normal liver and lung in single-cell atlases

pdac_liver_markers = read.delim('Misc/cell_types_markers_main_subtypes_liver.tsv', header = T, sep = '\t')
pdac_liver_markers = pdac_liver_markers[pdac_liver_markers$cluster %in% 'all' & pdac_liver_markers$avg_log2FC > 0,]
pdac_liver_markers = trimws(pdac_liver_markers$gene)

pdac_lung_markers = read.delim('Misc/cell_types_markers_main_subtypes_lung.tsv', header = T, sep = '\t')
pdac_lung_markers = pdac_lung_markers[pdac_lung_markers$cluster %in% 'all' & pdac_lung_markers$avg_log2FC > 0,]
pdac_lung_markers = trimws(pdac_lung_markers$gene)

lung_smpls = round(table(s_objs$orig.ident[s_objs$condition %in% 'lung'])/sum(table(s_objs$orig.ident[s_objs$condition %in% 'lung']))*100,2)
liver_smpls = round(table(s_objs$orig.ident[s_objs$condition %in% 'liver'])/sum(table(s_objs$orig.ident[s_objs$condition %in% 'liver']))*100,2)


tbl_liv = tbl_lng = list()
for(smpl_ in unique(s_objs$orig.ident)) 
{
  message('sample out: ',smpl_)
  
  s_ = subset(s_objs, subset = orig.ident != smpl_)
  markers_ = FindMarkers(object = s_,
                         group.by = 'condition', ident.1 = 'lung', ident.2 = 'liver')
  markers_ = markers_[order(abs(markers_$avg_log2FC),decreasing = T),]
  markers_$gene = rownames(markers_)
  
  liver_vs_lung = markers_[markers_$avg_log2FC < 0, ]
  lung_vs_liver = markers_[0 < markers_$avg_log2FC,]
  
  idxs_ = which(pdac_liver_markers %in% liver_vs_lung$gene)
  frac_liver = ceiling(length(idxs_)/length(pdac_liver_markers)*100)
  
  idxs_ = which(pdac_lung_markers %in% lung_vs_liver$gene)
  frac_lung = ceiling(length(idxs_)/length(pdac_lung_markers)*100)
  
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
write.table(x = tbl_liv, file = 'PDAC-liver_signature_leave-one-out.tsv', quote = F, sep = '\t', row.names = F, col.names = T)

tbl_lng = do.call(tbl_lng, what = rbind)
tbl_lng = tbl_lng[order(tbl_lng$fraction_in_lung, decreasing = T),]
write.table(x = tbl_lng, file = 'PDAC-lung_signature_leave-one-out.tsv', quote = F, sep = '\t', row.names = F, col.names = T)
