# BothDonors T cell annotations

source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
library(scCustomize)
# import the T cell object ====

cpal <- cols25(length(levels(Tcells@meta.data$ident)))
names(cpal) <- levels(Tcells@meta.data$ident)

#Tcells <- readRDS('./10x/objects/Ali_bothDonors_Tcells.rds') #This was processed the same way, the cluster names were just a bit different
Tcells <- readRDS( Set2_integrated_Tcells_path )
DimPlot(Tcells, label = T) 
Tcells@meta.data$ident <- Tcells@active.ident

# start by annotating clusters ----
umaps1 <- (DimPlot(Tcells, label = T) + scale_color_manual(values =cpal) |  
  ( DimPlot(Tcells, group.by = 'tissue') + scale_color_manual(values =  TissuePal ) ) )/ 
  (( DimPlot(Tcells, group.by = 'donor') + scale_color_manual(values =  DonorPal ) ) |
     DimPlot(Tcells, group.by = 'Tcell_type'))

#Feature plots
BasicT <- c('CD3E', 'CD8B', 'CD4', 'TRAC', 'TRDC')
UnconT <-c('KLRB1', "TRAV1-2", "ZBTB16", "ZNF683", "GNLY", "TYROBP")
SubsetsT <-c('CCR7', 'FOXP3', "TBX21", "GATA3", 'KLRG1') 
Feats <-c('CD8B', 'CD4', 'CCR7', 'SELL',"PDCD1", 'CXCR5', 'FOXP3', "KLRG1", 'KLRD1',
          "PRF1","CD27",'KLRB1',"KLF2",'IL7R','CD44','LEF1','TCF7' , 'FOXP1','CD248')

FPlots <- ViridisFeatures(Tcells, feature_args = Feats, 0.7, option = "D")
FP_wrap <- wrap_plots( FPlots )
VP_wrap <- VlnPlot(Tcells, c('PDCD1','CXCR5'))

clus <- Tcells@meta.data %>%
  select(ident) %>%
  group_by_all() %>%
  tally()

#new idents
new.cluster.ids <- c("CD4/CD8_naive",#0
                     "CD4/CD8_memory",#1
                     "CD4_effector",#2
                     "CD4/CD8_naive",#3
                     "CD8_CTL", #4
                     "CD4/CD8_naive", #5
                     "CD4_effector", #6
                     "CD4/CD8_memory", #7
                     "CD4/CD8_naive", #8
                     "CD4_effector", #9
                     "CD4/CD8_memory", #10
                     "CD4/CD8_ISG", #11
                     "CD4/CD8_naive", #12
                     "CD4/CD8_memory", #13
                     "Tfh/Treg", #14
                     "CD8_TEMRA", #15
                     "CD4_effector", #16
                     "MAIT/NKT", #17
                     "CD4/CD8_memory", #18
                     "CD4/CD8_memory", #19
                     "cycling"
                     )

#Tcells <- SetIdent(Tcells, value = 'orig.ident')
Tcells@meta.data$orig.ident <- Tcells@active.ident
names(new.cluster.ids) <- as.character(seq(new.cluster.ids) - 1)
Tcells <- RenameIdents(Tcells, new.cluster.ids)
Tcells@meta.data$ident2 <- Tcells@active.ident

# group by T cell types
naives <- new.cluster.ids[which(new.cluster.ids %like% 'naive')]
eff_mem <- new.cluster.ids[which(new.cluster.ids %like% 'eff' |
                                  new.cluster.ids %like% 'cycling'|
                                  new.cluster.ids %like% 'CTL'|
                                  new.cluster.ids %like% 'ISG'|
                                  new.cluster.ids %like% 'memory'|
                                   new.cluster.ids %like% 'MAIT' |
                                   new.cluster.ids %like% 'TEMRA' )]
Tfh <- new.cluster.ids[which(new.cluster.ids %like% 'Tfh/Treg')]


#assign new column
Tcells@meta.data$Tcell_type <- case_when(
  Tcells@meta.data$ident2 %in% naives ~ "naive",
  Tcells@meta.data$ident2 %in% eff_mem ~ "effector/memory",
  Tcells@meta.data$ident2 %in% Tfh ~ "Tfh/Treg"
)

#reorder the ident levels
clusterorder <- c( naives, eff_mem, Tfh) %>% as.character(.) %>% unique() %>% sort()
Tcells@meta.data$Tcell_type <-  factor( Tcells@meta.data$Tcell_type, levels =  c('naive','effector/memory','Tfh/Treg') )
DimPlot(Tcells, group.by = 'Tcell_type')

# Find T cell cluster markers ====


TMarkers <- read.csv('./10x/outs/integrated_markers.csv')

TMarkers5 <- TMarkers %>% 
  group_by(cluster) %>% 
  slice_max( order_by = avg_log2FC , n =5) %>%
  pull(gene) %>%
  as.character(.) %>%
  unique()

Tcells <- SetIdent(Tcells,value = 'ident')
scp_custom <-Clustered_DotPlot(Tcells, TMarkers5, colors_use_exp = viridis(50))

#DoHeatmap(Tcells, TMarkers5 ) + scale_fill_viridis(option = 'c')

ggsave('./10x/outs/Set2_heatmap_markers.png', width = 10, height = 12)
pdf(file = "./10x/outs/Set2_heatmap_markers.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 10)
scp_custom[[2]]
dev.off()
png(file = "./10x/outs/Set2_heatmap_markers.png",   # The directory you want to save the file in
    width = 8,
    units = 'in', # The width of the plot in inches
    height = 10,
    res = 300)
scp_custom[[2]]
dev.off()
#save updated T cell object ----

saveRDS(Tcells , Set2_integrated_Tcells_path )

