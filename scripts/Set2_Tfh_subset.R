
source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
library(future)
library(harmony)
session_mem <- 1000 * 1024^2
options(future.globals.maxSize= session_mem)

#read data ===

CD4_Tcells <- readRDS(Set2_integrated_CD4_path)
Clonedf <- read.delim( clone_df_path , stringsAsFactors = F)

old_TFH <- readRDS(Tfh_lineages_path)

old_TFH_cells <- paste0('TwoYear_',colnames(old_TFH))
IL10_cells <- old_TFH@meta.data %>%
  filter(Tfh_type == 'IL10 TFH') %>%
  rownames() %>%
  paste0('TwoYear_',.)



meta <- Tcells@meta.data %>%
  filter(barcode %in% colnames(CD4_Tcells)) %>%
  select(Tcell_type)
nrow(meta)==ncol(CD4_Tcells)
  
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
(FeaturePlot(CD4_Tcells, c('PDCD1','ICOS','CXCR5','FOXP3'),raster = T) | freq | DimPlot(CD4_Tcells, label = T,raster = T) )

# ID the Tfh cells ====


# known_TFH <- CD4_Tcells@meta.data %>%  #TFH global
#   filter(Tcell_type == 'Tfh') %>%
#   pull(barcode)

cluster_TFH <- CD4_Tcells@meta.data %>% #TFH local
  filter(ident %in% c(9,10)) %>%
  pull(barcode)

#TFH_cells <- union(union(cluster_TFH, known_TFH), TFH_cells) #enforce old TFH

TFH_cells <- union(cluster_TFH, old_TFH_cells)
                   
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


Tfh_cluster <- subset(CD4_Tcells, cells = TFH_cells_in )

# rerun harmony ====
Tfh_cluster <- Tfh_cluster %>% 
  FindVariableFeatures() 

pc_genes <- VariableFeatures(Tfh_cluster)[
  -which(VariableFeatures(Tfh_cluster) %in% IgTcr$gene)]

Tfh_cluster <- RunPCA(Tfh_cluster, features = pc_genes, npcs = 50) 


Tfh_cluster_harmonized <- RunHarmony(Tfh_cluster, 
                                     assay.use = "SCT", 
                                     group.by.vars = c("dataset","donor"))


Tfh_cluster_harmonized <- RunUMAP(Tfh_cluster_harmonized, dims=1:50, reduction = "harmony")
Tfh_cluster_harmonized <- FindNeighbors(Tfh_cluster_harmonized, reduction = "harmony") %>% 
  FindClusters(., resolution = 1.1)
Tfh_cluster_harmonized@meta.data$ident <- Tfh_cluster_harmonized@active.ident



# rerun harmony ====
Tfh_cluster <- subset(Tfh_cluster_harmonized, subset = UMAP_1 > -10 )

Tfh_cluster <- Tfh_cluster %>% 
  FindVariableFeatures() 

pc_genes <- VariableFeatures(Tfh_cluster)[
  -which(VariableFeatures(Tfh_cluster) %in% IgTcr$gene)]

Tfh_cluster <- RunPCA(Tfh_cluster, features = pc_genes, npcs = 50) 


Tfh_cluster_harmonized <- RunHarmony(Tfh_cluster, 
                                     assay.use = "SCT", 
                                     group.by.vars = c("dataset","donor"))


Tfh_cluster_harmonized <- RunUMAP(Tfh_cluster_harmonized, dims=1:50, reduction = "harmony")
Tfh_cluster_harmonized <- FindNeighbors(Tfh_cluster_harmonized, reduction = "harmony") %>% 
  FindClusters(., resolution = 0.9)
Tfh_cluster_harmonized@meta.data$ident <- Tfh_cluster_harmonized@active.ident




# stats====

stats <- Tfh_cluster_harmonized@meta.data %>%
  select(dataset, donor, time_point, tissue) %>%
  group_by_all() %>%
  tally()
write.csv(stats, './10x/outs/integrated_TFH_stats.csv')


Tfh_marks <- FindAllMarkers(Tfh_obj)
write.csv(Tfh_marks, './10x/outs/integrated_TFH_markers.csv')


top_marks <- Tfh_marks %>%
  filter(gene %notin% c('MTRNR2L12','AC016831.7') & grepl('^RP[SL]', gene ) == F) %>%
  group_by(cluster) %>%
  slice_max(n=5, order_by = avg_log2FC ) %>% 
  pull(gene) %>% 
  append(., 'IL10')

DoHeatmap(Tfh_obj, features = top_marks) + scale_fill_viridis(option = 'B')
ggsave('./10x/outs/Set2_Tfh_lineages_marker_hm.pdf', width = 12, height = 8)


table(TFH_cells %in% colnames(Tfh_cluster_harmonized))

#module scoring ----
module_genes <- read.delim('./10x/outs/TwoYear_Tfh_subset_module_markers.tsv') %>%
  group_by(subset) %>%
  group_split()

module_gene_sets <-list()
for(i in seq_along(module_genes)){
  module_gene_sets[[i]] <- module_genes[[i]]$gene
  names(module_gene_sets)[i] <- module_genes[[i]]$subset %>% unique()
}

Tfh_cluster_harmonized <- AddModuleScore(Tfh_cluster_harmonized, module_gene_sets, assay = 'RNA' )

module_cols <- which(grepl('Cluster[123]$', colnames(Tfh_cluster_harmonized@meta.data)))
colnames(Tfh_cluster_harmonized@meta.data)[module_cols] <- paste0(gsub(' TFH','', names(module_gene_sets)), '_module')

module_score_df <- FetchData(Tfh_cluster_harmonized, 
                             c(colnames(Tfh_cluster_harmonized@meta.data)[module_cols],'ident','donor')
                             ) %>%
  pivot_longer(1:3, values_to = 'score', names_to = 'subset')


ggplot(module_score_df, aes(ident, score, fill = subset) ) + geom_boxplot() 



#classificaiton ----

VlnPlot(Tfh_cluster_harmonized, c('CCR7','IL10','FOXP3', 'VIM','PVALB','PDCD1', 'CXCR5'))

Tfh_cluster_harmonized@meta.data$Tfh_type <- case_when(
  Tfh_cluster_harmonized@meta.data$ident %in% c(2,4,7,14,16,19) ~ 'GC',
  Tfh_cluster_harmonized@meta.data$ident %in% c(13) ~ 'IL10 TFH',
  Tfh_cluster_harmonized@meta.data$ident %in% c(5,9,10) ~ 'Treg',
  .default = 'pre/memory'
)

#plots ===


freq <- Tfh_cluster_harmonized@meta.data %>%
  select(ident, Tcell_type) %>%
  group_by(ident) %>%
  add_tally(name = 'n_ident') %>%
  ungroup() %>%
  group_by(ident, Tcell_type) %>%
  add_tally(name = 'n_ident_set') %>%
  distinct_all() %>%
  mutate(freq = n_ident_set/ n_ident) %>%
  ggplot(.,aes(ident, freq, fill = Tcell_type)) + geom_col()
(FeaturePlot(Tfh_cluster_harmonized, c('PDCD1','IL10','FOXP3','CCR7','CXCR5','CXCR4')) | freq) 




# mark the Tfh cells in the clonedf====
Clonedf$Tfh_lineage <- ifelse(Clonedf$barcode %in% rownames(Tfh_cluster_harmonized@meta.data), 'Tfh', NA)

meta_df<- Tfh_cluster_harmonized@meta.data %>%
  select(barcode, Tfh_type)

Clonedf <- left_join(select(Clonedf, -Tfh_type), meta_df)

write_tsv(Clonedf, clone_df_path )



#save ====

saveRDS(Tfh_cluster_harmonized , Set2_integrated_Tfh_path )



GC_marks <- FindMarkers(Tfh_cluster_harmonized, ident.1 = 'GC', ident.2 = 'pre/memory', group.by = 'Tfh_type')
IL10_marks <- FindMarkers(Tfh_cluster_harmonized, ident.1 = 'IL10 TFH', ident.2 = 'pre/memory', group.by = 'Tfh_type')
Treg_marks <- FindMarkers(Tfh_cluster_harmonized, ident.1 = 'Treg', ident.2 = 'pre/memory', group.by = 'Tfh_type')


fill_missing_vals <- function(df) {
  df %>%
    select(-N, -n) %>%
    pivot_wider(names_from = 'time_point', values_from = 'freq', values_fill = 0) %>%
    pivot_longer(4:ncol(.), names_to= 'time_point', values_to = 'freq' ) %>%
    mutate(time_point = factor(time_point, levels = names(TimePal2)[which(names(TimePal2) %in% time_point)] ))
}
freq_plots <- function(df){
  ggplot(df, aes(time_point, freq, color = Tfh_type, group=Tfh_type)) + 
    geom_point(size = 3) + geom_line(linewidth =1) + 
    scale_color_manual(values = TfhPal) +
    theme_bw() +
    ggtitle(paste(unique(df$donor), unique(df$tissue)))
  
}


TFH_type_freqs <- Tfh_cluster_harmonized@meta.data %>%
  select(donor, tissue, time_point, Tfh_type) %>%
  filter(!is.na(time_point)) %>%
  group_by(donor, tissue, time_point) %>%
  add_tally(name = 'N') %>%
  group_by_all() %>%
  add_tally(name = 'n') %>%
  distinct_all() %>%
  mutate(freq = n/N) %>%
  group_by(donor, tissue) %>%
  group_split()

TFH_type_freqs <- map(TFH_type_freqs, ~fill_missing_vals(.))

freq_plots_all <- map(TFH_type_freqs, ~ freq_plots(.))


clone_df <- makeClonedf(Tfh_cluster_harmonized, other_metadata = c('donor','tissue','time_point','Tfh_type'))


cross_type_clones <- clone_df %>% 
  filter(tissue =='FNA') %>% 
  group_by(donor, clone_id,Tfh_type) %>% 
  tally() %>% 
  distinct_all() %>%
  group_by(clone_id ) %>%
  add_tally(name = 'n_types') %>%
  filter(n_types > 2)
  
  
  pivot_wider(names_from = 'time_point', values_from = 'freq', values_fill = 0) %>%
  pivot_longer(4:ncol(.), names_to= 'time_point', values_to = 'freq' )


