
source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
library(future)
library(harmony)
library(scCustomize)
#plan(multisession, workers = 4)
#session_mem <- 1000 * 1024^2
#options(future.globals.maxSize= session_mem, future.seed=TRUE)

#read data ===

CD4_Tcells <- readRDS(Set2_integrated_CD4_path)


Clonedf <- read.delim( clone_df_path , stringsAsFactors = F) 
Clonedf <- Clonedf%>%
  select(-Tfh_type,-Tfh_lineage)

old_TFH <- readRDS(Tfh_lineages_path)

old_TFH_cells <-colnames(old_TFH)
IL10_cells <- old_TFH@meta.data %>%
  filter(Tfh_type == 'IL10 TFH') %>%
  rownames() %>%
  paste0('TwoYear_',.)


CD4_Tcells <- AddMetaData(CD4_Tcells, meta)
#pull out tfh =====


freq <- CD4_Tcells@meta.data %>%
  select(ident, Tcell_type) %>%
  group_by(ident) %>%
  add_tally(name = 'n_ident') %>%
  ungroup() %>%
  group_by(ident, Tcell_type) %>%
  add_tally(name = 'n_ident_set') %>%
  distinct_all() %>%
  mutate(freq = n_ident_set/ n_ident) %>%
  ggplot(.,aes(ident, freq, fill = Tcell_type)) + geom_col()
(FeaturePlot(CD4_Tcells, c('PDCD1','ICOS','CXCR5','FOXP3'),raster = T, order = T)) | 
  (DimPlot(CD4_Tcells, label = T,raster = T) / 
     DimPlot(CD4_Tcells, group.by = 'tissue', raster = T))
freq

FeaturePlot(CD4_Tcells, c('FOXP3'), split.by = 'ident',pt.size = 1, ncol = 5, order = T)
# ID the Tfh cells ====


# known_TFH <- CD4_Tcells@meta.data %>%  #TFH global
#   filter(Tcell_type == 'Tfh') %>%
#   pull(barcode)

cluster_TFH <- FetchData(CD4_Tcells, c('ident','IL21','CXCL13'))%>% #TFH local
  filter(ident %in% c(13:14,20) ) %>% #| CXCL13 != 0 | IL21 !=0
  rownames()

#TFH_cells <- union(union(cluster_TFH, known_TFH), TFH_cells) #enforce old TFH

Tfh_picked_bc <- Clonedf %>% 
  filter(!is.na(Tfh_clone_id)) %>% 
  pull(barcode)

TFH_cells <- union(union(cluster_TFH, old_TFH_cells), Tfh_picked_bc)
                   
TFH_clones <- CD4_Tcells@meta.data %>%
  filter(barcode %in% TFH_cells & !is.na(clone_id)) %>%
  pull(clone_id)


TFH_clones_cells <- CD4_Tcells@meta.data %>%
  filter(clone_id %in% TFH_clones) %>%
  pull(barcode)

TFH_cells_in <- union(TFH_clones_cells,TFH_cells)

TCR_pos <- CD4_Tcells@meta.data %>%
  filter(!is.na(clone_id)) %>%
  pull(barcode)

TFH_cells_in <- intersect(TFH_cells_in, TCR_pos)

naive_spike1 <- CD4_Tcells@meta.data %>% 
  filter(ident == 12 & tissue =='PBMC') %>% 
  sample_n(size = 200) %>% 
  pull(barcode)

naive_spike2 <- CD4_Tcells@meta.data %>% 
  filter(ident == 12 & tissue =='FNA') %>% 
  sample_n(size = 200) %>% 
  pull(barcode)

naive_spike <- append(naive_spike1, naive_spike2)

TFH_cells_in <- append(naive_spike, TFH_cells_in)
Tfh_cluster <- subset(CD4_Tcells, cells = unique(TFH_cells_in) )

# rerun harmony ====
Tfh_cluster <- Tfh_cluster %>% 
  FindVariableFeatures() 

pc_genes <- VariableFeatures(Tfh_cluster)[
  -which(VariableFeatures(Tfh_cluster) %in% IgTcr$gene)]

Tfh_cluster <- RunPCA(Tfh_cluster, features = pc_genes, npcs = 50) 


Tfh_cluster_harmonized <- RunHarmony(Tfh_cluster, 
                                     assay.use = "SCT", 
                                     group.by.vars = c('dataset','donor','time','tissue'))
ElbowPlot(Tfh_cluster_harmonized, reduction = 'harmony')
Tfh_cluster_harmonized <- RunUMAP(Tfh_cluster_harmonized, dims=1:10, reduction = "harmony")
Tfh_cluster_harmonized <- FindNeighbors(Tfh_cluster_harmonized, reduction = "harmony")
Tfh_cluster_harmonized <- FindClusters(Tfh_cluster_harmonized, resolution = 0.6)
Tfh_cluster_harmonized@meta.data$ident <- Tfh_cluster_harmonized@active.ident


(DimPlot_scCustom(Tfh_cluster_harmonized, label = T) | 
  DimPlot(Tfh_cluster_harmonized, cells.highlight = naive_spike))/ 
  (DimPlot(Tfh_cluster_harmonized, cells.highlight = old_TFH_cells) |
  DimPlot(Tfh_cluster_harmonized, cells.highlight = Tfh_picked_bc))
#Tfh_marksx <- FindAllMarkers(Tfh_cluster_harmonized)
FeaturePlot_scCustom(Tfh_cluster_harmonized, c('CXCR5','PDCD1','IL21','CXCL13','FOXP3','IL10','CCR7','TCF7','TOX2'))
VlnPlot(Tfh_cluster_harmonized, c('CXCR5','PDCD1','IL21','CXCL13','FOXP3','IL10','CCR7','MKI67','TOX2'))

DimPlot(Tfh_cluster_harmonized, cells.highlight = naive_spike)
DimPlot(Tfh_cluster_harmonized, group.by = 'tissue')
FeaturePlot(Tfh_cluster_harmonized, c('MKI67'))



Tfh_cluster_harmonized@meta.data$Tfh_picked_X <- ifelse(colnames(Tfh_cluster_harmonized) %in% Tfh_picked_bc, 'yes',NA)
Tfh_cluster_harmonized@meta.data %>% select(Tfh_picked_X, ident) %>% group_by_all() %>% tally()
# tmp <- FetchData(Tfh_cluster_harmonized, c('clone_id','ident')) 
# tmp_l <- tmp %>% group_by(ident) %>% group_split()
# 
# map_dbl(tmp_l, ~length(.$clone_id[which(tmp_l[[9]]$clone_id %in% .$clone_id )]))
# c8_marks <- FindMarkers(Tfh_cluster_harmonized, ident.1 = 8)

# rerun harmony ====

cluster_TFH2 <- FetchData(Tfh_cluster_harmonized, c('ident','IL21','CXCL13'))%>% #TFH local
  filter(ident %notin% c(7,2) ) %>% #| CXCL13 != 0 | IL21 !=0
  rownames()

#TFH_cells <- union(union(cluster_TFH, known_TFH), TFH_cells) #enforce old TFH

TFH_cells2 <- union(union(cluster_TFH2, naive_spike), Tfh_picked_bc)

Tfh_cluster2 <- subset(Tfh_cluster_harmonized, cells =TFH_cells2 )

Tfh_cluster2 <- Tfh_cluster2 %>% 
  FindVariableFeatures() 

pc_genes <- VariableFeatures(Tfh_cluster2)[
  -which(VariableFeatures(Tfh_cluster2) %in% IgTcr$gene)]

Tfh_cluster2 <- RunPCA(Tfh_cluster2, features = pc_genes, npcs = 50) 


Tfh_cluster_harmonized2 <- RunHarmony(Tfh_cluster2, 
                                     assay.use = "SCT", 
                                     group.by.vars = c('dataset','donor','time','tissue'))

Tfh_cluster_harmonized2 <- RunUMAP(Tfh_cluster_harmonized2, dims=1:15, reduction = "harmony")
Tfh_cluster_harmonized2 <- FindNeighbors(Tfh_cluster_harmonized2, reduction = "harmony") 
Tfh_cluster_harmonized2 <- FindClusters(Tfh_cluster_harmonized2, resolution = 0.6)
Tfh_cluster_harmonized2@meta.data$ident <- Tfh_cluster_harmonized2@active.ident

(DimPlot_scCustom(Tfh_cluster_harmonized2) | 
  DimPlot(Tfh_cluster_harmonized2, cells.highlight = old_TFH_cells))/ 
  (DimPlot(Tfh_cluster_harmonized2, cells.highlight = naive_spike)| 
  DimPlot(Tfh_cluster_harmonized2, cells.highlight = Tfh_picked_bc))

FeaturePlot_scCustom(Tfh_cluster_harmonized2, c('CXCR5','PDCD1','IL21','CXCL13','FOXP3','IL10','CCR7','MKI67','TOX2'))
VlnPlot(Tfh_cluster_harmonized2, c('CXCR5','PDCD1','IL21','CXCL13','FOXP3','IL10','CCR7','MKI67','TOX2'))

# rerun harmony ====
# cluster_TFH3 <- FetchData(Tfh_cluster_harmonized2, c('ident','IL21','CXCL13'))%>% #TFH local
#   filter(ident %notin% c(0,4,11,13) ) %>% #| CXCL13 != 0 | IL21 !=0
#   rownames()
# #TFH_cells <- union(union(cluster_TFH, known_TFH), TFH_cells) #enforce old TFH
# 
# TFH_cells3 <- union(union(cluster_TFH3, naive_spike), Tfh_picked_bc)
# 
# Tfh_cluster3 <- subset(Tfh_cluster_harmonized2, cells =TFH_cells3 )
# 
# Tfh_cluster3 <- Tfh_cluster3 %>% 
#   FindVariableFeatures() 
# 
# pc_genes <- VariableFeatures(Tfh_cluster3)[
#   -which(VariableFeatures(Tfh_cluster3) %in% IgTcr$gene)]
# 
# Tfh_cluster3 <- RunPCA(Tfh_cluster3, features = pc_genes, npcs = 50) 
# 
# 
# Tfh_cluster_harmonized3 <- RunHarmony(Tfh_cluster3, 
#                                       assay.use = "SCT", 
#                                       group.by.vars = c('dataset','donor'))
# 
# 
# Tfh_cluster_harmonized3 <- RunUMAP(Tfh_cluster_harmonized3, dims=1:50, reduction = "harmony")
# Tfh_cluster_harmonized3 <- FindNeighbors(Tfh_cluster_harmonized3, reduction = "harmony") 
# Tfh_cluster_harmonized3 <- FindClusters(Tfh_cluster_harmonized3, resolution = 0.6)
# Tfh_cluster_harmonized3@meta.data$ident <- Tfh_cluster_harmonized3@active.ident
# (DimPlot_scCustom(Tfh_cluster_harmonized3) | 
#   DimPlot_scCustom(Tfh_cluster_harmonized3, group.by = 'tissue')) / 
#   (DimPlot(Tfh_cluster_harmonized3, cells.highlight = old_TFH_cells)| 
#     DimPlot(Tfh_cluster_harmonized3, cells.highlight = naive_spike)) 
#    DimPlot(Tfh_cluster_harmonized3, cells.highlight = Tfh_picked_bc)
# 
# FeaturePlot_scCustom(Tfh_cluster_harmonized3, c('CXCR5','PDCD1','IL21','CXCL13','FOXP3','IL10','CCR7','MKI67','TOX2'))
# VlnPlot(Tfh_cluster_harmonized3, c('CXCR5','PDCD1','IL21','CXCL13','FOXP3','IL10','CCR7','MKI67','TOX2'))

# stats====

stats <- Tfh_cluster_harmonized2@meta.data %>%
  select(dataset, donor, time_point, tissue) %>%
  group_by_all() %>%
  tally()
write.csv(stats, './10x/outs/integrated_TFH_stats.csv')


Tfh_marks <- FindAllMarkers(Tfh_cluster_harmonized2)
write.csv(Tfh_marks, './10x/outs/integrated_TFH_markers.csv')


top_marks <- Tfh_marks %>%
  filter(gene %notin% c('MTRNR2L12','AC016831.7') & grepl('^RP[SL]', gene ) == F) %>%
  group_by(cluster) %>%
  slice_max(n=5, order_by = avg_log2FC ) %>% 
  pull(gene) %>% 
  append(., 'IL10')

DoHeatmap(Tfh_cluster_harmonized2, features = top_marks) + scale_fill_viridis(option = 'B')
ggsave('./10x/outs/Set2_Tfh_lineages_marker_hm.pdf', width = 12, height = 8)



#module scoring ----
# module_genes <- read.delim('./10x/outs/TwoYear_Tfh_subset_module_markers.tsv') %>%
#   group_by(subset) %>%
#   group_split()
# 
# module_gene_sets <-list()
# for(i in seq_along(module_genes)){
#   module_gene_sets[[i]] <- module_genes[[i]]$gene
#   names(module_gene_sets)[i] <- module_genes[[i]]$subset %>% unique()
# }
# 
# Tfh_cluster_harmonized3 <- AddModuleScore(Tfh_cluster_harmonized3, module_gene_sets, assay = 'RNA' )
# 
# module_cols <- which(grepl('Cluster[123]$', colnames(Tfh_cluster_harmonized3@meta.data)))
# colnames(Tfh_cluster_harmonized3@meta.data)[module_cols] <- paste0(gsub(' TFH','', names(module_gene_sets)), '_module')
# 
# module_score_df <- FetchData(Tfh_cluster_harmonized3, 
#                              c(colnames(Tfh_cluster_harmonized3@meta.data)[module_cols],'ident','donor')
# ) %>%
#   pivot_longer(1:3, values_to = 'score', names_to = 'subset')
# 
# 
# ggplot(module_score_df, aes(ident, score, fill = subset) ) + geom_boxplot() 
# 

#classificaiton ----

DimPlot(Tfh_cluster_harmonized2, label = T)
VlnPlot(Tfh_cluster_harmonized2, c('VIM','IL21','FOXP3', 'FOXB1','PVALB','PDCD1', 'CXCR5','TOX2','CXCL13','MKI67','SELL','TMSB10'))

Tfh_cluster_harmonized2@meta.data$Tfh_type <- case_when(
  Tfh_cluster_harmonized2@meta.data$ident %in% c(0) ~ 'naive',
  Tfh_cluster_harmonized2@meta.data$ident %in% c(1,5,7,8) ~ 'GC',
  Tfh_cluster_harmonized2@meta.data$ident %in% c(6) ~ 'IL10 TFH',
  Tfh_cluster_harmonized2@meta.data$ident %in% c(11) ~ 'cycling',
  Tfh_cluster_harmonized2@meta.data$ident %in% c(2,9) ~ 'Treg',
  .default = 'pre/memory'
)

DimPlot(Tfh_cluster_harmonized2, group.by = 'Tfh_type')

#plots ===


freq <- Tfh_cluster_harmonized2@meta.data %>%
  select(ident, Tcell_type) %>%
  group_by(ident) %>%
  add_tally(name = 'n_ident') %>%
  ungroup() %>%
  group_by(ident, Tcell_type) %>%
  add_tally(name = 'n_ident_set') %>%
  distinct_all() %>%
  mutate(freq = n_ident_set/ n_ident) %>%
  ggplot(.,aes(ident, freq, fill = Tcell_type)) + geom_col()

(FeaturePlot_scCustom(Tfh_cluster_harmonized3, c('PDCD1','IL10','FOXP3','CCR7','CXCR5','TOX2')) | freq) 




# mark the Tfh cells in the clonedf====
#Clonedf$Tfh_lineage <- ifelse(Clonedf$barcode %in% rownames(Tfh_cluster_harmonized@meta.data), 'Tfh', NA)

meta_df<- Tfh_cluster_harmonized2@meta.data %>%
  select(barcode, Tfh_type) %>%
  mutate(Tfh_lineage = 'Tfh')

Clonedf <- left_join(Clonedf, meta_df)

write_tsv(Clonedf, clone_df_path )



#save ====

saveRDS(Tfh_cluster_harmonized2 , Set2_integrated_Tfh_path )

