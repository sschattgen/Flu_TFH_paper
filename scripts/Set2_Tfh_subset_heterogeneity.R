# BothDonor clonotype analysis 
# the alluvial plots of clonal lineage abundance for the picked clones is in here, too

source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
library(harmony)
library(scCustomize)
library(rlist)
#
Tfh_obj <- readRDS(Set2_integrated_Tfh_path)


calc_time_freq_by_cluster <- function(df){
  
  df %>%
    select(donor, ident, time, year) %>%
    group_by(donor, time) %>%
    add_tally(name = 'N') %>%
    group_by(donor, time, ident) %>%
    add_tally(name = 'n') %>%
    ungroup() %>% 
    distinct_all() %>%
    mutate(freq = n/N) %>% 
    arrange(ident)
}

obj_list <- list()
marker_list <- list()
plot_list <- list()

for (typex in c('pre/memory','GC')){
  
  Tfh_obj2 <- subset(Tfh_obj, subset = Tfh_type == typex)
  
  
  Tfh_obj2 <- FindVariableFeatures(Tfh_obj2) 
  
  pc_genes <- VariableFeatures(Tfh_obj2)[
    -which(VariableFeatures(Tfh_obj2) %in% IgTcr$gene)]
  
  Tfh_obj2 <- RunPCA(Tfh_obj2, features = pc_genes, npcs = 50) 
  
  Tfh_obj2@meta.data$donor <-  as.character(Tfh_obj2@meta.data$donor)
  
  Tfh_cluster_harmonized <- RunHarmony(Tfh_obj2, 
                                       assay.use = "SCT", 
                                       group.by.vars = c("dataset","donor",'time', 'tissue'))
  
  
  Tfh_cluster_harmonized <- Tfh_cluster_harmonized %>%
    RunUMAP(., dims=1:10, reduction = "harmony") %>%
    FindNeighbors(., reduction = "harmony")
  
  Tfh_cluster_harmonized <-  FindClusters(Tfh_cluster_harmonized, resolution = 0.6)
  
  Tfh_cluster_harmonized@meta.data$ident <- Tfh_cluster_harmonized@active.ident
  
  if (typex == 'pre/memory'){
    Tfh_cluster_harmonized <- subset(Tfh_cluster_harmonized, idents = '7', invert = TRUE)
    
    Tfh_cluster_harmonized <- FindVariableFeatures(Tfh_cluster_harmonized)  
      
    pc_genes <- VariableFeatures(Tfh_cluster_harmonized)[
      -which(VariableFeatures(Tfh_cluster_harmonized) %in% IgTcr$gene)]
    
    Tfh_cluster_harmonized <- RunPCA(Tfh_cluster_harmonized, features = pc_genes, npcs = 50) 
    
    Tfh_cluster_harmonized@meta.data$donor <-  as.character(Tfh_cluster_harmonized@meta.data$donor)
    
    Tfh_cluster_harmonized <- RunHarmony(Tfh_cluster_harmonized, 
                                         assay.use = "SCT", 
                                         group.by.vars = c("dataset","donor",'time', 'tissue'))
    
    
    Tfh_cluster_harmonized <- Tfh_cluster_harmonized %>%
      RunUMAP(., dims=1:10, reduction = "harmony") %>%
      FindNeighbors(., reduction = "harmony")
    
    Tfh_cluster_harmonized <-  FindClusters(Tfh_cluster_harmonized, resolution = 0.6)
    
    Tfh_cluster_harmonized@meta.data$ident <- Tfh_cluster_harmonized@active.ident
    
  }
  
  Tfh_cluster_marks <- FindAllMarkers(Tfh_cluster_harmonized)
  
  Clustered_DotPlot(Tfh_cluster_harmonized, Extract_Top_Markers(Tfh_cluster_marks, num_genes = 8, make_unique = T,named_vector = F))
  p1 <- DimPlot(Tfh_cluster_harmonized, group.by = 'time') + scale_color_manual(values = TimePal2) + guides(color=guide_legend(ncol=2))
  p2 <- DimPlot(Tfh_cluster_harmonized, group.by = 'year') + scale_color_manual(values = YearPal) 
  p3 <- DimPlot(Tfh_cluster_harmonized, group.by = 'tissue') + scale_color_manual(values = TissuePal)
  p4 <- DimPlot(Tfh_cluster_harmonized, group.by = 'donor') + scale_color_manual(values = DonorPal)
  
  layout <- (p3 + p4 + p2 + p1) + plot_annotation(title = typex) & NoAxes()
  
  x_out <- DimPlot_scCustom(Tfh_cluster_harmonized,pt.size = 1) + layout
  
  obj_list[[typex]] <- Tfh_cluster_harmonized
  marker_list[[typex]] <- Tfh_cluster_marks
  plot_list[[typex]] <- layout
  
}

DimPlot_scCustom(obj_list[[1]]) + ggtitle('pre/memory')
DimPlot_scCustom(obj_list[[2]]) + ggtitle('GC')

Clustered_DotPlot(obj_list[[1]], Extract_Top_Markers(marker_list[[1]], num_genes = 8, make_unique = T,named_vector = F))
Clustered_DotPlot(obj_list[[2]], Extract_Top_Markers(marker_list[[2]], num_genes = 8, make_unique = T,named_vector = F))


#cluster distrobution across time ====
cpal <- DiscretePalette(10,palette = 'polychrome')
names(cpal) <- seq(10)-1
freq_df1 <- calc_time_freq_by_cluster(obj_list[[1]]@meta.data) %>%
  mutate(time = factor(time, levels = rev(names(TimePal2))))

ggplot( freq_df1,
        aes( freq,time, fill=ident)) +
  geom_col() +
  scale_fill_manual(values = cpal) +
  facet_wrap( ~ donor, scales = 'free_y', nrow = 1)+
  ggtitle('pre/memory TFH')
ggsave('./10x/outs/Set2_Tfh_lineages_premem_subset_distro.png', dpi = 300, width = 10, height = 4)


freq_df2 <- calc_time_freq_by_cluster(obj_list[[2]]@meta.data) %>%
  mutate(time = factor(time, levels = rev(names(TimePal2))))
ggplot( freq_df2,
        aes( freq,time, fill=ident)) +
  geom_col() +
  scale_fill_manual(values = cpal) +
  facet_wrap( ~ donor, scales = 'free_y', nrow = 1) +
  ggtitle('GC TFH')
ggsave('./10x/outs/Set2_Tfh_lineages_GC_subset_distro.png', dpi = 300, width = 10, height = 4)



saveRDS(obj_list,'./10x/objects/intergrated_TFH_PM_GC_split_out.rds')


# density and umaps
Clustered_DotPlot(obj_list[[1]], Extract_Top_Markers(marker_list[[1]], num_genes = 8, make_unique = T,named_vector = F))
premem_density <- Plot_Density_Custom(seurat_object = obj_list[[1]], features = c('CXCR5','CXCR3','CCR6','CCR7','IL7R','CD27','PDCD1','BCL6','RORC','GATA3','TBX21'), pt.size = 1.5, combine = F)
premem_density <-map(premem_density, ~.+NoLegend())
premem_density[[12]] <- DimPlot_scCustom(obj_list[[1]], pt.size = 2, label.size = 7) +ggtitle('GEX clusters') + NoLegend()
wrap_plots(premem_density) + plot_layout(ncol =3)
ggsave('./10x/outs/Set2_Tfh_lineages_premem_subset_feature_plots.pdf', dpi = 300, width = 7.5, height = 10)
DimPlot_scCustom(obj_list[[1]], pt.size = 2, label.size = 7) +ggtitle('GEX clusters')
ggsave('./10x/outs/Set2_Tfh_lineages_premem_subset_umap.pdf',   width = 5, height = 5)

VlnPlot(obj_list[[1]], c('RORC','CXCR3','TBX21','CCR6'), group.by = 'tissue')
# density and umaps
Clustered_DotPlot(obj_list[[2]], Extract_Top_Markers(marker_list[[2]], num_genes = 8, make_unique = T,named_vector = F))

gc_density <- Plot_Density_Custom(seurat_object = obj_list[[2]], 
                                  features = c('IL4','IL21','CXCL13','TNFRSF18','LAG3','PVALB','MAP3K8','GNG4','DUSP4','TOX2','BCL6'), pt.size = 1.5, combine = F)
gc_density <-map(gc_density, ~.+NoLegend())
gc_density[[12]] <- DimPlot_scCustom(obj_list[[2]],  pt.size = 2, label.size = 7) + ggtitle('GEX clusters') + NoLegend()
wrap_plots(gc_density) + plot_layout(ncol =3)
ggsave('./10x/outs/Set2_Tfh_lineages_gc_subset_feature_plots.pdf', dpi = 300, width = 7.5, height = 10)
DimPlot_scCustom(obj_list[[2]],  pt.size = 2, label.size = 7) + ggtitle('GEX clusters') 
ggsave('./10x/outs/Set2_Tfh_lineages_gc_subset_umap.pdf',    width = 5, height = 5)









#psuedotime

cds <- as.cell_data_set(obj_list[[1]])
fData(cds)$gene_short_name <- rownames(fData(cds))

recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions


list.cluster <- obj_list[[1]]@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- obj_list[[1]]@reductions$umap@cell.embeddings



cluster.before.traj <- plot_cells(cds, 
                                  color_cells_by = "cluster", 
                                  label_groups_by_cluster = F, 
                                  group_label_size = 5) + 
  theme(legend.position = "right")

FeaturePlot_scCustom(obj_list[[1]], c('MKI67','CCR7','CXCR5','PDCD1'))

cds <- learn_graph(cds)
cds <- order_cells(cds)
p1 <- plot_cells(cds, color_cells_by = "pseudotime", 
                 label_groups_by_cluster = F,
                 label_branch_points =F, 
                 label_roots =  T, 
                 label_leaves = F, 
                 trajectory_graph_segment_size = 1)

p2 <-plot_cells(cds, color_cells_by = "Tfh_type", 
                label_groups_by_cluster = F,
                label_branch_points =F, 
                label_roots =  F, 
                label_leaves = F,
                group_label_size = 0, 
                trajectory_graph_segment_size = 1) + 
  scale_color_manual(values = TfhPal)


p3 <-plot_cells(cds, color_cells_by = "year", 
                label_groups_by_cluster = F,
                label_branch_points =F, 
                label_roots =  F, 
                label_leaves = F,
                group_label_size = 0, 
                trajectory_graph_segment_size = 1) + 
  scale_color_manual(values = YearPal)
DimPlot(obj_list[[1]], group.by = 'time')+ scale_color_manual(values = TimePal2)
p1 |p2 |p3

