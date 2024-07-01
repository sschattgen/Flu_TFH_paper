
source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
library(harmony)
library(monocle3)
library(SeuratWrappers)
# TFH ====


Tfh_obj <- readRDS(Set2_integrated_Tfh_path )
cds <- readRDS('./10x/objects/intergrated_TFH_lymNode_monocle.rds')

# rerun harmony ====
Tfh_cluster <- subset(Tfh_obj, subset = Tfh_type %notin% c('Treg','high ISG','cycling')
                      & tissue =='FNA')

Tfh_cluster <- Tfh_cluster %>% 
  FindVariableFeatures() 

pc_genes <- VariableFeatures(Tfh_cluster)[
  -which(VariableFeatures(Tfh_cluster) %in% IgTcr$gene)]

Tfh_cluster <- RunPCA(Tfh_cluster, features = pc_genes, npcs = 50) 

Tfh_cluster@meta.data$donor <-  as.character(Tfh_cluster@meta.data$donor)

Tfh_cluster_harmonized <- RunHarmony(Tfh_cluster, 
                                     assay.use = "SCT", 
                                     group.by.vars = c("dataset","donor",'time'))


Tfh_cluster_harmonized <- Tfh_cluster_harmonized %>%
  RunUMAP(., dims=1:10, reduction = "harmony") %>%
  FindNeighbors(., reduction = "harmony")
  
Tfh_cluster_harmonized <-  FindClusters(Tfh_cluster_harmonized, resolution = 0.6)

Tfh_cluster_harmonized@meta.data$ident <- Tfh_cluster_harmonized@active.ident



# monocle. reload below ====

cds <- as.cell_data_set(Tfh_cluster_harmonized)
fData(cds)$gene_short_name <- rownames(fData(cds))

recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions


list.cluster <- Tfh_cluster_harmonized@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- Tfh_cluster_harmonized@reductions$umap@cell.embeddings



cluster.before.traj <- plot_cells(cds, 
                                 color_cells_by = "cluster", 
                                 label_groups_by_cluster = F, 
                                 group_label_size = 5) + 
  theme(legend.position = "right")


cds <- learn_graph(cds)
cds <- order_cells(cds)



# OR skip here and reload ======


#cds <- readRDS('./10x/objects/intergrated_TFH_lymNode_monocle.rds')

#plots====
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


p4 <-plot_cells(cds, color_cells_by = "time_point", 
                 label_groups_by_cluster = F,
                 label_branch_points =F,
                 label_roots =  F, 
                 label_leaves = F,
                 group_label_size = 0, 
                 trajectory_graph_segment_size = 1) + 
  scale_color_manual(values = TimePal2)

p5 <-plot_cells(cds, color_cells_by = "donor", 
                label_groups_by_cluster = F,
                label_branch_points =F, 
                label_roots =  F, 
                label_leaves = F,
                group_label_size = 0, 
                trajectory_graph_segment_size = 1) + 
  scale_color_manual(values = DonorPal)



(p5|p2)/(p3|p4)/ (p1 | plot_spacer())

umap_list <- list(
  'pseudotime' = p1,
  'pseudotime_noLegend' = p1 +theme(legend.position = 'none'),
  'Tfh_type_noLegend' = p2,
  'Tfh_type' = p2+theme(legend.position = 'right'),
  'year_noLegend' = p3,
  'year' = p3+theme(legend.position = 'right'),
  'time_point_noLegend' = p4,
  'time_point' = p4+theme(legend.position = 'right'),
  'donor_noLegend' = p5,
  'donor' = p5+theme(legend.position = 'right')
)
for (i in seq_along(umap_list)){
  ggsave(paste0('./10x/outs/Set2_monocle_Tfh_noTreg_',names(umap_list)[i],'_umap.pdf'),
         plot = umap_list[[i]],
         width = 6,
         height = 6)
  
}

# custom pseudotime plotting =====

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))


example_genes <- c('CCR7','CXCR5','PDCD1','PVALB','TOX2','IL21','MAP3K8','PRDM1','FOXB1','IL10')
gex_data <- FetchData(Tfh_cluster_harmonized, slot = 'data', c('barcode', example_genes))
data.pseudo2 <- left_join(data.pseudo, gex_data)%>%
  filter(Tfh_type != 'naive')


plot_df <- data.pseudo2 %>% 
  select(barcode, Tfh_type, monocle3_pseudotime, all_of(example_genes)) %>%
  pivot_longer(all_of(example_genes), names_to = 'gene', values_to = 'expression') %>%
  mutate(gene = factor(gene, levels = example_genes)) 

z1 <- ggplot(plot_df, aes(monocle3_pseudotime, expression)) + 
  geom_point(aes(color = Tfh_type), alpha = 0.3) + 
  geom_smooth(aes(color = Tfh_type))+
  facet_wrap(~gene, scales = 'free_y', ncol = 2) + 
  scale_color_manual(values = TfhPal) + 
  theme_minimal() +
  xlab('pseudotime') +
  ylab('log normalized expression') + 
  labs(color = 'Tfh subset') + 
  theme(legend.position = 'top',
        panel.grid = element_blank())

ggsave('./10x/outs/Set2_monocle_marker_gene_vs_pseudotime_type_fit.pdf',
       plot = z1,width = 8, height =8)
ggsave('./10x/outs/Set2_monocle_marker_gene_vs_pseudotime_type_fit.png',
       plot = z1,width = 8, height = 8)

z2 <- ggplot(plot_df, aes(monocle3_pseudotime, expression)) + 
  geom_point(aes(color = Tfh_type), alpha = 0.3) + 
  geom_smooth()+
  facet_wrap(~gene, scales = 'free_y', ncol = 2) + 
  scale_color_manual(values = TfhPal) + 
  theme_minimal() +
  xlab('pseudotime') +
  ylab('log normalized expression') + 
  labs(color = 'Tfh subset') + 
  theme(legend.position = 'top',
        panel.grid = element_blank())

ggsave('./10x/outs/Set2_monocle_marker_gene_vs_pseudotime_all_fit.pdf',
       plot = z2,width = 8, height = 8)
ggsave('./10x/outs/Set2_monocle_marker_gene_vs_pseudotime_all_fit.png',
       plot = z2,width = 8, height = 8)





z3 <- ggplot(data.pseudo2, aes(time_point, monocle3_pseudotime, fill = donor)) + 
  geom_boxplot() + 
  scale_fill_manual(values = DonorPal) + 
  theme_bw() +
  ylab('psuedotime') +
  xlab('time point') + 
  theme_bw() +
  facet_wrap(~donor, scales = 'free_x')
  

ggsave('./10x/outs/Set2_monocle_marker_time_point_boxplot.pdf',
       plot = z3, width = 18, height = 6)





data.pseudo4 <- data.pseudo2 %>%
  filter( donor =='321-05') %>%
  filter(time_point %notin% c('y1_d12','y2_d14', 'y1_d90')) %>% #'y2_d90',
  mutate( time_X = factor(case_when(
    time_point %in% c('y1_d0','y2_d0') ~ 'pre',
    time_point %in% c('y1_d5','y2_d7') ~ 'early',
    time_point %in% c('y1_d28','y2_d28') ~ 'mid',
    time_point %in% c('y1_d60','y2_d60') ~ 'late',
    time_point %in% c('y2_d120','y1_d180') ~'very late'
  ), levels = c('pre','early','mid','late','very late')) )


line_point <- data.pseudo4 %>% 
  group_by(Tfh_type) %>%
  summarise(median(monocle3_pseudotime))

line_point_c <- line_point[[2]]
names(line_point_c) <- line_point[[1]]
z4 <- ggplot(data.pseudo4, aes(monocle3_pseudotime)) + 
  geom_vline(xintercept = line_point_c['pre/memory'], color = TfhPal[1]) +
  geom_vline(xintercept = line_point_c['GC'], color = TfhPal[2]) +
  geom_vline(xintercept = line_point_c['IL10 TFH'], color = TfhPal[3]) +
  geom_density(aes(color = year, fill = year), adjust = 2, linewidth = 1, alpha = 0.3) + 
  #geom_histogram(position = position_stack()) + 
  facet_wrap(~time_X, ncol =  1) + 
  theme_minimal() +
  scale_color_manual(values = YearPal) +
  scale_fill_manual(values = YearPal) +
  theme(legend.position = 'bottom', 
        axis.text.y = element_blank(),
        panel.grid = element_blank()) +
  xlab('pseudotime') 

ggsave('./10x/outs/Set2_monocle_d5_picked_time_point_density.pdf',
       plot = z4, width = 262/72, height = 575/72)


z5 <-ggplot(data.pseudo2, aes(monocle3_pseudotime)) + 
  geom_density(aes(color = Tfh_type, fill = Tfh_type), adjust = 1, linewidth = 1, alpha = 0.3) + 
  #geom_histogram(position = position_stack()) + 
  theme_minimal() +
  scale_color_manual(values = TfhPal) +
  scale_fill_manual(values = TfhPal) +
  theme(legend.position = 'top', 
        axis.text.y = element_blank(),
        panel.grid = element_blank()) +
  xlab('pseudotime') +
  labs(fill = 'Tfh type', color = 'Tfh type')

ggsave('./10x/outs/Set2_monocle_all_donor_types_density.pdf',
       plot = z5, width = 262/72, height = 2)



# save monocle object ===

saveRDS(cds,'./10x/objects/intergrated_TFH_lymNode_monocle.rds')


