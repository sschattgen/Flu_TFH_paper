# TwoYear Tfh subsetting

source('scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

# import T cells Seurat object and df ----
Tcells <- readRDS(Tcells_path)
Clonedf <- read.delim(clone_df_path, stringsAsFactors = F)

# find all cells that have a matching TCR with clones in the Tfh cluster and subset the object ====

Tfh_clones <- Clonedf %>% 
  filter(Tcell_type == "Tfh") %>%
  pull(clone_id) %>% 
  unique()

Tfh_cl_cells <- Clonedf %>% 
  filter( clone_id %in% Tfh_clones & 
            ident %notlike% 'CD8') %>% 
  pull(barcode)

Clonedf <- Clonedf %>% 
  mutate(Tfh_lineage = 
           ifelse(barcode %in%  Tfh_cl_cells, 'yes', 'no') 
         )

Tcells@meta.data$Tfh_lineage <- ifelse(Tcells@meta.data$barcode %in%  Tfh_cl_cells,
                                       'yes', 'no')

Tfh_lin_Grp <- Tcells@meta.data %>% 
  group_by(ident, Tfh_lineage ) %>% 
  tally() %>% 
  mutate(freq = n / sum(n)) %>% 
  mutate_if(is.character, as.factor)

Tfh_lin_Grp_plt <- ggplot(Tfh_lin_Grp, aes(ident, freq, fill = Tfh_lineage )) + 
  geom_col(position = position_fill(reverse = TRUE)) + 
  theme_minimal() +
  labs(title = "Freq of Tfh lineage cells by cluster") + 
  theme(axis.text.x = element_text(angle = 90))

ggsave('outs/TwoYear_Tcell_Tfh_lineage_cluster_freq.png',
       plot = Tfh_lin_Grp_plt,
       height = 5, width =8)

# subset and rerun seurat workflow ====
Tfh_obj <- subset(Tcells, cells = Tfh_cl_cells) 
Tfh_obj <- subset(Tfh_obj, CD8B == 0)
Tfh_obj <- NormalizeData(object = Tfh_obj)
Tfh_obj <- FindVariableFeatures(object = Tfh_obj)
Tfh_obj_var_genes <- VariableFeatures(Tfh_obj)
Tfh_obj_PC_genes <- Tfh_obj_var_genes[-which(Tfh_obj_var_genes %in% IgTcr$genes)] 
Tfh_obj <- ScaleData(object = Tfh_obj, 
                     features = rownames(Tfh_obj))
Tfh_obj <- RunPCA(object = Tfh_obj, features = Tfh_obj_PC_genes )
Tfh_obj <- FindNeighbors(object = Tfh_obj)
Tfh_obj <- FindClusters(object = Tfh_obj)
Tfh_obj <- RunUMAP(object = Tfh_obj, dims = 1:20)

#stash ident
Tfh_obj@meta.data$ident <- Tfh_obj@active.ident

#mark the cells in the large T cell object
tfh_type <- rep("Tfh", nrow(Tfh_obj@meta.data))
names(tfh_type) <- rownames(Tfh_obj@meta.data)
Tcells <- AddMetaData(Tcells, tfh_type, 'Tfh_type')
Tcells@meta.data$Tfh_type[which( is.na(Tcells@meta.data$Tfh_type) ) ] <- 'other'


# cell count by tissue
Tissue <- FetchData(Tfh_obj , vars = c('tissue')) %>% 
  as.data.frame() %>% 
  group_by_all() %>% 
  tally()

# add the Tfh_clone_id
picked <- Clonedf %>%
  filter(barcode %in% colnames(Tfh_obj)) %>%
  select(barcode, Tfh_clone_id)
pickedc <- picked[['Tfh_clone_id']]
names(pickedc) <-picked[['barcode']]
pickedc[which(is.na(pickedc))] <- 'other'
Tfh_obj <- AddMetaData( Tfh_obj, pickedc, 'Tfh_clone_id')


# annotate the clusters by phenotype ====

# Find marker genes
Tfh_cl_markers <- FindAllMarkers(Tfh_obj)
write_tsv( Tfh_cl_markers, 'outs/TwoYear_Tfh_lineages_markers.tsv')
Tfh_cl_markers <- read.delim('outs/TwoYear_Tfh_lineages_markers.tsv', stringsAsFactors = F)

# marker heatmap
Topmarks_Tfh <- Tfh_cl_markers %>% 
  filter(gene %notin% IgTcr$genes) %>% 
  group_by(cluster) %>% 
  top_n(5, avg_log2FC)

cpal <- cols25(12)
names(cpal) <- seq(12) - 1

DoHeatmap(Tfh_obj, features = Topmarks_Tfh$gene, group.by = "ident") + 
  scale_fill_viridis(option="inferno") 

ggsave('outs/TwoYear_Tfh_lineages_heatmap_markers.pdf', width = 10, height = 8, 
       useDingbats = F)

#annotate each cluster as early or GC Tfh, other CD4 or Treg
(DimPlot(Tfh_obj, label = T) | DimPlot(Tfh_obj, group.by = 'tissue') + scale_color_manual(values = TissuePal)) / 
  FeaturePlot(Tfh_obj, c('PDCD1','CXCR5', 'VIM', 'EMP3', 'TOX2','FOXB1','IL10')) + 
  plot_layout(heights = c(1,3))


#define the cells for each group
Tfh_md_df <- FetchData(Tfh_obj, c('ident', 'clone_id', 'barcode','FOXB1','IL10','CXCL13','TOX2','CXCR5'))

EM_cells <- Tfh_md_df %>%
  filter(ident %in% c(1,2,9,10)) %>%
  pull(barcode)
EM_cells_v <- rep('pre/memory', length(EM_cells))
names(EM_cells_v) <- EM_cells
Tfr_c7 <- Tfh_md_df %>%
  filter(ident %in% 7) %>%
  pull(barcode)
Tfr_c8 <- Tfh_md_df %>%
  filter(ident %in% 8) %>%
  filter(FOXB1 > 0 | IL10 > 0) %>%
  pull(barcode)
Tfr_cells <- c(Tfr_c7, Tfr_c8)
Tfr_cells_v <- rep('IL10 TFH', length(Tfr_cells))
names(Tfr_cells_v) <- Tfr_cells
GC_cells <- Tfh_md_df %>%
  filter(barcode %notin% Tfr_cells) %>%
  filter(ident %notin% c(1,2,9,10)) %>%
  pull(barcode)
GC_cells_v <- rep('GC', length(GC_cells))
names(GC_cells_v) <- GC_cells

# add to metadata
Tfh_type <- c(EM_cells_v, Tfr_cells_v, GC_cells_v)
Tfh_obj <- AddMetaData( Tfh_obj, Tfh_type, 'Tfh_type')

Tfh_obj@meta.data$Tfh_type <- factor(Tfh_obj@meta.data$Tfh_type, levels = c('pre/memory','GC','IL10 TFH'))
# mark the Tfh cells in the clonedf====
tfh_cells <- rownames(Tfh_obj@meta.data)
Clonedf$Tfh_type <- ifelse(Clonedf$barcode %in% tfh_cells, 'Tfh', NA)
write_tsv(Clonedf, clone_df_path )


#save ====
saveRDS(Tfh_obj, Tfh_lineages_path )
saveRDS(Tcells, Tcells_path )

