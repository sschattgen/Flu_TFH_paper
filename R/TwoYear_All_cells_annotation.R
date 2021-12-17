source('Z:/ResearchHome/ResearchHomeDirs/thomagrp/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

Both_sc <-readRDS( all_cells_path)
Both_sc@meta.data$ident <- Both_sc@active.ident


# make umap showing clonotype information====
cinfo_df <- FetchData(Both_sc, c('UMAP_1','UMAP_2','clone_info'))
cinfo_pal <- c('grey80', 'darkblue', 'darkred','green')
names(cinfo_pal) <- c('none','TCR','BCR', 'both')

cinfo_umap <- ggplot(cinfo_df)+
  rasterise( geom_point( data = subset(cinfo_df, clone_info == 'none'),
                         aes(UMAP_1, UMAP_2, color = clone_info), size = 0.3)
  ) +
  rasterise( geom_point( data = subset(cinfo_df, clone_info != 'none'),
                         aes(UMAP_1, UMAP_2, color = clone_info), size = 0.3) 
  ) +
  scale_color_manual(values = cinfo_pal) + 
  theme_minimal()
ggsave('./10x/outs/TwoYear_TCR_IG_clone_umap.pdf',
       plot = cinfo_umap, height = 6, width = 6.5, useDingbats = F)

# mark cell types ====
#dotplot of lineage markers
markers <- list(
  c('BCL6', 'RGS13', 'MEF2B', 'STMN1', 'ELL3', 'SERPINA9'), #GC B cells
  c('XBP1', 'IRF4', 'SEC11C', 'FKBP11', 'JCHAIN', 'PRDM1'), #plasmablasts
  c('TCL1A', 'IL4R', 'CCR7', 'IGHM', 'IGHD'), #naive B cells
  c('TBX21', 'FCRL5', 'ITGAX', 'NKG7', 'ZEB2', 'CR2' ), # activated B cells, CR2 is a negative marker
  c('TNFRSF13B', 'CD27', 'CD24'), #resting memory B cells
  c('MS4A1', 'CD19', 'CD79A','CD2'), # B cells
  c('CD3D', 'CD3E', 'CD3G', 'IL7R', 'CD4', 'CD8A'), # T cells
  c('GZMB', 'GNLY', 'NCAM1'), #NK cells
  c('CD14', 'LYZ', 'CST3', 'MS4A7'), #mono/granulocyte
  c('IL3RA', 'CLEC4C'), #plasmacyotid
  c('PPBP') #platelets
)

names(markers) <- c('GC B cells','plasmablasts','naive B cells', 'active B cells', 'memory B cells',
                    'pan B cells', 'pan T cells', 'NK cells', 'myeloid','plasmacytoid', 'platelets')
Both_sc_DP <- DotPlot( Both_sc , features = markers) 
Both_sc_DP <- Both_sc_DP +  
  scale_color_viridis() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        pan)
ggsave('./10x/outs/TwoYear_all_cells_lineage_markers.pdf',
       plot = Both_sc_DP, height = 9, width = 18, useDingbats = F)


# mark clusters
ident <- Both_sc@meta.data$ident
Both_sc@meta.data$cell_type <- case_when(
     ident == 0 ~'T cell',
     ident == 1 ~'B cell',
     ident == 2 ~'B cell',
     ident == 3 ~'T cell',
     ident == 4 ~'T cell',
     ident == 5 ~'T cell',
     ident == 6 ~'B cell',
     ident == 7 ~'T cell',
     ident == 8 ~'B cell',
     ident == 9 ~'B cell',
     ident == 10 ~'B cell',
     ident == 11 ~'T cell',
     ident == 12 ~'myeloid',
     ident == 13 ~'B cell',
     ident == 14 ~'T cell',
     ident == 15 ~'NK cell',
     ident == 16 ~'T cell',
     ident == 17 ~'B cell',
     ident == 18 ~'T cell',
     ident == 19 ~'B cell',
     ident == 20 ~'doublets',
     ident == 21 ~'B cell',
     ident == 22 ~'B cell',
     ident == 23 ~'B cell',
     ident == 24 ~'myeloid',
     ident == 25 ~'myeloid',
     ident == 26 ~'pDC',
     ident == 27 ~'doublets',
     ident == 28 ~'T cell',
     ident == 29 ~'B cell',
     ident == 30 ~'megakaryocyte',
     ident == 31 ~'myeloid'
     )

DimPlot(Both_sc, group.by = 'cell_type' ) 
ggsave('./10x/outs/TwoYear_all_cells_lineage_umap.pdf',
       height = 6, width = 6.5, useDingbats = F)

# subset the B cells ====

Bcells <- subset( Both_sc, 
                  subset = cell_type == 'B cell' & CD3E == 0 & LYZ == 0 )

Bcells <- SeuratReprocess( Bcells,  cluster_resolution = 0.7)
Bcells@meta.data$ident <- Bcells@active.ident

# umap and dotplot
DimPlot(Bcells, label = T ) / 
  #DimPlot(Bcells, group.by = 'Donor', split.by = 'Tissue') /
  DotPlot( Bcells , features = markers[1:5]) +  
  scale_color_viridis() + 
  NoLegend()

ggsave('./10x/outs/TwoYear_Bcells_clone_info_umap.pdf',
       height = 6, width = 6.5, useDingbats = F)

# plot of cells with clone information
Bcells@meta.data$bcr_clone_info <- ifelse( !is.na(Bcells@meta.data$bcr_clone_id), 
                                           'yes', 'no')

DimPlot(Bcells, group.by = 'bcr_clone_info')

ggsave('./10x/outs/TwoYear_Bcells_clone_info_umap.pdf',
       height = 6, width = 6.5, useDingbats = F)

# annotate B cell types
ident <- Bcells@meta.data$ident 
Bcells@meta.data$B_cell_type <- case_when(
  ident == 0 ~'naive',
  ident == 1 ~'RMB',
  ident == 2 ~'RMB',
  ident == 3 ~'RMB',
  ident == 4 ~'PB',
  ident == 5 ~'PB',
  ident == 6 ~'RMB',
  ident == 7 ~'RMB',
  ident == 8 ~'RMB',
  ident == 9 ~'activated',
  ident == 10 ~'PB',
  ident == 11 ~'GC',
  ident == 12 ~'RMB',
  ident == 13 ~'PB',
  ident == 14 ~'activated',
  ident == 15 ~'GC',
  ident == 16 ~'PB',
  ident == 17 ~'PB',
  ident == 18 ~'activated',
  ident == 19 ~'activated'
)


DimPlot(Bcells, group.by = 'B_cell_type' )
ggsave('./10x/outs/TwoYear_Bcells_types_umap.pdf',
       height = 6, width = 6.5, useDingbats = F)

#save ====
saveRDS(Both_sc, Both_sc_path)
saveRDS(Bcells, Bcells_path)
