# set the path to the top directory
tfh_working_dir <-'Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/flu_vaxx_tfh'

library(data.table, quietly =T)
library(tidyverse, quietly =T)
library(Seurat, quietly =T)
library(patchwork, quietly =T)
library(jcolors, quietly =T)
library(viridis, quietly =T)
library(ggthemes, quietly =T)
library(ggpubr, quietly =T)
library(rlist, quietly =T)
library(ggrepel, quietly =T)
library(ggplotify, quietly =T)
library(wesanderson, quietly =T)
library(ggalluvial, quietly =T)
library(cowplot, quietly =T)
library(pals, quietly =T)
library(Peptides, quietly = T)
library(pheatmap, quietly = T)
library(ggplotify, quietly = T)

`%notin%` <- Negate(`%in%`)
`%notlike%` <- Negate(`%like%`)

# paths 
all_cells_path <- paste0( tfh_working_dir, '/data/TwoYear_all_cells.rds')
Tcells_path <- paste0( tfh_working_dir, '/data/TwoYear_Tcells.rds')
Bcells_path <- paste0( tfh_working_dir, '/data/TwoYear_Bcells.rds')
Tfh_lineages_path <- paste0( tfh_working_dir, '/data/TwoYear_Tfh_lineages.rds')
Tfh_only_path <- paste0( tfh_working_dir, '/data/TwoYear_Tfh_only.rds')
clone_df_path <- paste0( tfh_working_dir,'/data/TwoYear_10x_clones.tsv')
IgTcr <- read.delim( "data/Human_10X_IgTcr.tsv", stringsAsFactors = F)

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

flu_pal <- c('grey80', 'red')
names(flu_pal)<- c('unknown', 'yes')

#functions used in the workflow ====

makeClonedf <- function( obj, stripAlleles = FALSE, other_metadata = NULL ) {
  
  # other_metadata: additional metadata can be indicated with a character string
  
  cols <- c('barcode', 'ident', 'clone_id', 'va_gene', 'ja_gene', 'vb_gene', 'jb_gene', 
            'cdr3a', 'cdr3b', 'cdr3a_nucseq', 'cdr3b_nucseq', "UMAP_1", "UMAP_2" )
  
  if ( is.null(other_metadata) == TRUE ){
    
    df <- FetchData(obj, vars = cols ) 
    
  } else {
    
    new_cols <- append( cols, other_metadata )
    df <- FetchData(obj, vars = new_cols ) 
    
  }
  
  if ( stripAlleles == TRUE ) {
    
    print("removing allele information")
    
    df$va_gene <- str_split( df$va_gene , pattern = "[*]", simplify = T)[,1]
    df$ja_gene <- str_split( df$ja_gene , pattern = "[*]", simplify = T)[,1]
    df$vb_gene <- str_split( df$vb_gene , pattern = "[*]", simplify = T)[,1]
    df$jb_gene <- str_split( df$jb_gene , pattern = "[*]", simplify = T)[,1]
    
  }  else{
    
    print( "allele information preserved")
    
  }
  
  return(df)
  
}

countClones <- function(df){
  
  # adds as columns to clone_df contain clone size by cluster and in total
  
  df2 <- df %>% 
    group_by(ident) %>% 
    add_count(clone_id, name = "count_cluster") %>% 
    mutate(freq_cluster = count_cluster /sum(count_cluster))
  
  df3 <- df2 %>% 
    ungroup(ident) %>% 
    add_count(clone_id, name = "count_total") %>% 
    mutate(freq_total = count_total /sum(count_total))
  
  return(df3)
  
}

countTopClones <- function(df, clone_counts){
  
  #clone_counts is string of integers corresponding to the top clones eg 5, 10, etc
  #returns named list split clone_counts input
  
  olist <- list()
  
  for (i in seq_along(clone_counts)) {
    
    clones <- df %>% 
      select(barcode, clone_id) %>% 
      group_by(clone_id) %>% 
      tally() %>% 
      arrange(desc(n)) %>% 
      top_n(clone_counts[i])
    
    olist[[i]] <- clones
    names(olist)[i] <- names(clone_counts[i])
    
  }
  
  return(olist)
  
} 

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

ViridisFeatures <- function(obj, feature_args, size = 1, option = 'B', 
                            combined = FALSE, order = FALSE , reduction = 'umap', slot = 'data', strip_axis = TRUE, remove_legend = TRUE ){
  
  ###
  # feature_args should be a character string of gene names.
  ###
  stripped <- theme( axis.ticks = element_blank(), 
                     axis.line = element_blank(), 
                     axis.text = element_blank(), 
                     axis.title = element_blank())
  NoLegend <- theme(legend.position = 'none')
  
  pre <- list()
  
  for (i in seq_along(feature_args)) {
    
    p <- FeaturePlot(obj, features = feature_args[i], pt.size = size, order = order, reduction = paste(reduction), slot = paste(slot) )  
    p2 <- p + scale_color_viridis(option = paste(option))
    
    if (strip_axis == TRUE ){
      p2 <- p2 + stripped
    } 
    
    if (remove_legend == TRUE ){
      p2 <- p2 + NoLegend
    }
    
    pre[[i]] <- p2
    names(pre)[i] <- feature_args[i]
    
  }
  
  if (combined == TRUE ){
    
    board <- wrap_plots(pre)
    
    return(board)
    
  } else {
    
    return(pre)
    
  }
}

detectUnconventionals <- function(df, species){
  
  
  #calculates iMHC based on conga paper
  
  df2 <- df %>% 
    mutate(cdr3a_length = nchar(cdr3a),
           cdr3b_length = nchar(cdr3b)) %>%
    mutate(fgloop_cdr3a = str_sub(cdr3a, 5, -5),
           fgloop_cdr3b = str_sub(cdr3b, 5, -5)) %>%
    mutate(fgloop_cdr3a_length = nchar(fgloop_cdr3a),
           fgloop_cdr3a_aa_charge = as.numeric(charge(fgloop_cdr3a, pH = 7, pKscale = "Lehninger")),
           fgloop_cdr3a_cys_count =  str_count(fgloop_cdr3a, 'C'),
           fgloop_cdr3a_try_count =  str_count(fgloop_cdr3a, 'W'),
           fgloop_cdr3a_arg_count = str_count(fgloop_cdr3a, 'R'),
           fgloop_cdr3a_lys_count = str_count(fgloop_cdr3a, 'K'),
           fgloop_cdr3a_his_count = str_count(fgloop_cdr3a, 'H'),
           fgloop_cdr3a_asp_count = str_count(fgloop_cdr3a, 'D'),
           fgloop_cdr3a_glu_count = str_count(fgloop_cdr3a, 'E'),
           fgloop_cdr3b_length = nchar(fgloop_cdr3b),
           fgloop_cdr3b_aa_charge = as.numeric(charge(fgloop_cdr3b, pH = 7, pKscale = "Lehninger")),
           fgloop_cdr3b_cys_count =  str_count(fgloop_cdr3b, 'C'),
           fgloop_cdr3b_try_count =  str_count(fgloop_cdr3b, 'W'),
           fgloop_cdr3b_arg_count = str_count(fgloop_cdr3b, 'R'),
           fgloop_cdr3b_lys_count = str_count(fgloop_cdr3b, 'K'),
           fgloop_cdr3b_his_count = str_count(fgloop_cdr3b, 'H'),
           fgloop_cdr3b_asp_count = str_count(fgloop_cdr3b, 'D'),
           fgloop_cdr3b_glu_count = str_count(fgloop_cdr3b, 'E')
    ) #count cysteine and tryptophan 
  
  df2 <- df2 %>%
    mutate(MHCi_score_cdr3a = fgloop_cdr3a_length + 3*fgloop_cdr3a_cys_count + 2*fgloop_cdr3a_try_count + 
             fgloop_cdr3a_arg_count + fgloop_cdr3a_lys_count + 0.5*fgloop_cdr3a_his_count - 
             fgloop_cdr3a_asp_count - fgloop_cdr3a_glu_count, 
           MHCi_score_cdr3b = fgloop_cdr3b_length + 3*fgloop_cdr3b_cys_count + 2*fgloop_cdr3b_try_count + 
             fgloop_cdr3b_arg_count + fgloop_cdr3b_lys_count + 0.5*fgloop_cdr3b_his_count - 
             fgloop_cdr3b_asp_count - fgloop_cdr3b_glu_count 
           ) %>%
    mutate(MHCi = MHCi_score_cdr3a + 2*MHCi_score_cdr3b)

  #detects MAIT and NKT invariant alpha chains
  
  scall <- str_to_lower(species)
  
  if (any(str_detect(df2$va_gene, pattern = "[*]"))) {
    df3 <- df2
    df3$va_gene <- str_split( df2$va_gene , pattern = "[*]", simplify = T)[,1]
    df3$ja_gene <- str_split( df2$ja_gene , pattern = "[*]", simplify = T)[,1]
    df3$vb_gene <- str_split( df2$vb_gene , pattern = "[*]", simplify = T)[,1]
    df3$jb_gene <- str_split( df2$jb_gene , pattern = "[*]", simplify = T)[,1]
    
  } else{
    df3 <- df2
  }     
  
  
  #matching genes for MAITs and NKTs
  hmait_v <- c("TRAV1-2") 
  hmait_j <- c("TRAJ33","TRAJ20", "TRAJ12") 
  hnkt_v <- c("TRAV10") 
  hnkt_j <- c("TRAJ18") 
  hnkt_b <- c("TRBV25") 
  
  mmait_v <- c("TRAV1") 
  mmait_j <- c("TRAJ33") 
  mnkt_v <- c("TRAV11") 
  mnkt_j <- c("TRAJ18") 

  
  if ( scall == "mouse"){
    
    df4 <- df3 %>% 
      mutate( is_NKT =  ifelse( va_gene %in% mnkt_v & ja_gene %in% mnkt_j & cdr3a_length == 15, "yes", "no") ,
              is_MAIT =  ifelse( va_gene %in% mmait_v & ja_gene %in% mmait_j & cdr3a_length == 12, "yes", "no") )
    
  } else if  (scall == "human") {
    
    df4 <- df3 %>% 
      mutate( is_NKT =  ifelse( va_gene %in% hnkt_v & ja_gene %in% hnkt_j & vb_gene %in% hnkt_b, "yes", "no") ,
              is_MAIT =  ifelse( va_gene %in% hmait_v & ja_gene %in% hmait_j & cdr3a_length == 12, "yes", "no") )
    
  } else {
    
    print("Species not indicated or is unsupported")
    
  }
  df2$is_MAIT <- df4$is_MAIT
  df2$is_NKT <- df4$is_NKT
  return(df2)
  
}
