library(data.table)
library(tidyverse)
library(data.table)
library(Seurat)
library(patchwork)
library(jcolors)
library(viridis)
library(ggthemes)
library(ggpubr)
library(rlist)
library(ggrepel)
library(pheatmap)
library(ggplotify)
library(ggrastr)
library(wesanderson)
library(ggalluvial)
library(cowplot)
library(pals)

`%notin%` <- Negate(`%in%`)
`%notlike%` <- Negate(`%like%`)

# paths 
source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/TOOLBOX/TCR_helper_fxns.R')
source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/TOOLBOX/sc_analyses_custom.R')

tfh_working_dir <-'Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh'
all_cells_path <- paste0( tfh_working_dir, '/10x/objects/TwoYear_all_cells.rds')
Tcells_path <- paste0( tfh_working_dir, '/10x/objects/TwoYear_Tcells.rds')
Bcells_path <- paste0( tfh_working_dir, '/10x/objects/TwoYear_Bcells.rds')
Tfh_lineages_path <- paste0( tfh_working_dir, '/10x/objects/TwoYear_Tfh_lineages.rds')
Tfh_only_path <- paste0( tfh_working_dir, '/10x/objects/TwoYear_Tfh_only.rds')
clone_df_path <- paste0( tfh_working_dir,'/10x/outs/TwoYear_10x_clones.tsv')
IgTcr <- read.delim( "Z:/ResearchHome/Groups/thomagrp/home/sschattg/TOOLBOX/Human_10X_IgTcr.tsv", stringsAsFactors = F)

# color palettes
DonorPal <-jcolors('default')[c(1:3)]
names(DonorPal) <- c('321-04','321-05','321-11')

TissuePal <- c('#005B94', '#D12600')
names(TissuePal) <- c('FNA', 'PBMC')

TimePal <- viridis::viridis(10, option = "D")
names(TimePal) <- c('0','5', '7', '12','14', '28','60','90','120','180')

YearPal <- c("#FF8811","#9DD9D2")
names(YearPal) <- c('1','2')


time_levels <- c( "y1_d0", "y1_d5", "y1_d12", "y1_d28", "y1_d60", "y1_d90",  "y1_d120", "y1_d180", 
                  "y2_d0", "y2_d7", "y2_d14", "y2_d28", "y2_d60", "y2_d90", "y2_d120" )

TimePal2 <- viridis::viridis(15, option = "D")
names(TimePal2) <- time_levels

TfhPal <- jcolors('pal9')[1:3]
names(TfhPal) <- c('pre/memory', 'GC', 'IL10 TFH' )

good_conga_clusters <- c("10,1", "10,4", "8,0",  "8,14", "8,15", "8,18", "8,7",  "8,9" )
conga_cluster_pal <-c(jcolors("rainbow")[c(1,8,7,2,3,6,4,5)])
names(conga_cluster_pal) <- good_conga_clusters

clonePal <- c( jcolors('pal8')[-9], jcolors('pal7')[5] , 'grey90')
picked_ID <-c()
for (i in seq(12)){
  cip <- paste0('Tfh_321-05_clone_', i)
  
  picked_ID <- append(picked_ID, cip)
}
names(clonePal) <- c( picked_ID, 'other')

# Seurat reprocess wrapper

SeuratReprocess <- function(obj,  nfeatures = 2000, dims = 20, dim_reduction = 'umap', cluster_resolution = 0.8){
  
  # wrapper function useful for reprocessing after subsetting
  
  obj <- NormalizeData(object = obj)
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
  obj_var_genes <- VariableFeatures(obj)
  obj_PC_genes <- obj_var_genes[-which(obj_var_genes %in% IgTcr$genes)] 
  obj <- ScaleData(object = obj)
  obj <- RunPCA(object = obj, features = obj_PC_genes , npcs = dims )
  obj <- FindNeighbors(object = obj)
  obj <- FindClusters(object = obj, resolution = cluster_resolution)
  
  if( dim_reduction == 'umap'){
    obj <- RunUMAP(object = obj, dims = 1:dims)
  } else if (dim_reduction == 'tsne') {
    obj <- RunTSNE(object = obj, dims = 1:dims)
  } else if (dim_reduction == 'both') {
    obj <- RunUMAP(object = obj, dims = 1:dims)
    obj <- RunTSNE(object = obj, dims = 1:dims)
  }
  
  
  
  obj@meta.data$ident <- obj@active.ident
  
  return(obj)
  
}


