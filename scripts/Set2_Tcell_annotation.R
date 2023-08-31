# BothDonors T cell annotations

source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

# import the T cell object ====

#Tcells <- readRDS('./10x/objects/Ali_bothDonors_Tcells.rds') #This was processed the same way, the cluster names were just a bit different
Tcells <- readRDS( Set2_integrated_Tcells_path )
Tcells@meta.data$ident <- Tcells@active.ident

# start by annotating clusters ----
umaps1 <- (DimPlot(Tcells, label = T) |  
  ( DimPlot(Tcells, group.by = 'tissue') + scale_color_manual(values =  TissuePal ) ) )/ 
  (( DimPlot(Tcells, group.by = 'donor') + scale_color_manual(values =  DonorPal ) ) |
     DimPlot(Tcells, group.by = 'Tcell_type'))

#Feature plots
BasicT <- c('CD3E', 'CD8B', 'CD4', 'TRAC', 'TRDC')
UnconT <-c('KLRB1', "TRAV1-2", "ZBTB16", "ZNF683", "GNLY", "TYROBP")
SubsetsT <-c('CCR7', 'FOXP3', "TBX21", "GATA3", 'KLRG1') 
Feats <-c('CD8B', 'CD4', 'CCR7', 'SELL',"PDCD1", 'CXCR5', 'FOXP3', "KLRG1", 
          "PRF1","CD27",'KLRB1',"KLF2",'IL7R','CD44','LEF1','TCF7')
FPlots <- ViridisFeatures(Tcells, feature_args = Feats, 0.7, option = "D")
FP_wrap <- wrap_plots( FPlots )
VP_wrap <- VlnPlot(Tcells, c('PDCD1','CXCR5'))

clus <- Tcells@meta.data %>%
  select(ident) %>%
  group_by_all() %>%
  tally()

#new idents
new.cluster.ids <- c("CD4/CD8_naive",#0
                     "CD4_eff/mem",#1
                     "CD4_naive",#2
                     "CD4_naive",#3
                     "CD8_eff", #4
                     "CD4_eff/mem", #5
                     "CD4_eff/mem", #6
                     "CD4_eff/mem", #7
                     "CD4/CD8_naive", #8
                     "CD4/CD8_naive", #9
                     "TFH/Treg", #10
                     "MAIT/NKT", #11
                     "CD4/CD8_naive", #12
                     "CD4_eff/mem", #13
                     "CD4_eff/mem", #14
                     "CD4_eff/mem", #15
                     "TFH/Treg", #16
                     "CD4/CD8_eff/mem", #17
                     "CD4/CD8_eff/mem", #18
                     "CD4/CD8_naive", #19
                     "CD4/CD8_naive", #20
                     "cycling" #21
                     )

#Tcells <- SetIdent(Tcells, value = 'orig.ident')
Tcells@meta.data$orig.ident <- Tcells@active.ident
names(new.cluster.ids) <- as.character(seq(new.cluster.ids) - 1)
Tcells <- RenameIdents(Tcells, new.cluster.ids)
Tcells@meta.data$ident2 <- Tcells@active.ident

# group by T cell types
naives <- new.cluster.ids[which(new.cluster.ids %like% 'naive')]
eff_mem <- new.cluster.ids[which(new.cluster.ids %like% 'eff' |
                                   new.cluster.ids %like% 'cycling')]
Tfh <- new.cluster.ids[which(new.cluster.ids %like% 'TFH/Treg')]
MAIT <- new.cluster.ids[which(new.cluster.ids %like% 'MAIT')]

#assign new column
Tcells@meta.data$Tcell_type <- case_when(
  Tcells@meta.data$ident2 %in% naives ~ "naive",
  Tcells@meta.data$ident2 %in% eff_mem ~ "effector/memory",
  Tcells@meta.data$ident2 %in% Tfh ~ "Tfh",
  Tcells@meta.data$ident2 %in% MAIT ~ "MAIT/NKT"
)

#reorder the ident levels
clusterorder <- c( naives, MAIT, eff_mem, Tfh) %>% as.character(.) %>% unique() %>% sort()
Tcells@meta.data$Tcell_type <-  factor( Tcells@meta.data$Tcell_type, levels =  c('naive','MAIT/NKT','effector/memory','Tfh') )


# Find T cell cluster markers ====


TMarkers <- read.csv('./10x/outs/integrated_markers.csv')

TMarkers5 <- TMarkers %>% 
  group_by(cluster) %>% 
  slice_max( order_by = avg_log2FC , n =5) %>%
  pull(gene) %>%
  as.character(.)

DoHeatmap(Tcells, TMarkers5 ) + scale_fill_viridis(option = 'c')

ggsave('./10x/outs/Set2_heatmap_markers.png', width = 25, height = 15)



#save updated T cell object ----

Tcells <- SetIdent(Tcells, value = 'ident')
saveRDS(Tcells , Set2_integrated_Tcells_path )

