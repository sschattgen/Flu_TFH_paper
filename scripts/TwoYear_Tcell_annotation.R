# BothDonors T cell annotations

source('scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

# Rebuild the T cell only object ====
gexFiles <- Read10X_h5('data/TwoYear_Tcells_counts.h5')
Tcells <-  CreateSeuratObject( gexFiles )
Tcells@meta.data <- read.csv('data/TwoYear_Tcells_initial_metadata.csv', row.names = 1)
  
# normalize and var features
Tcells <- NormalizeData( Tcells)
Tcells <- FindVariableFeatures(object = Tcells)
Tcells_var_genes <- VariableFeatures(Tcells)
Tcells_PC_genes <- Tcells_var_genes[-which(Tcells_var_genes %in% IgTcr$genes)] 
#Tcells <- ScaleData(object = Tcells,  features = rownames(Tcells) )

## scoring cell cycle 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Tcells <- CellCycleScoring(Tcells, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Tcells$CC.Difference <- Tcells$S.Score - Tcells$G2M.Score

## but not regressing out here 
Tcells <- ScaleData(object = Tcells, 
                    features = rownames(Tcells) )

## dim reductions
Tcells <- RunPCA(object = Tcells, features = Tcells_PC_genes )
Tcells <- FindNeighbors(object = Tcells)
Tcells <- FindClusters(object = Tcells)
Tcells <- RunUMAP(object = Tcells , dims = 1:20)

# stash the idents and umis in metadata for later
Tcells@meta.data$ident <- Tcells@active.ident


# start by annotating clusters ----
umaps1 <- (DimPlot(Tcells, label = T) |  
  ( DimPlot(Tcells, group.by = 'tissue') + scale_color_manual(values =  TissuePal ) ) )/ 
  (( DimPlot(Tcells, group.by = 'donor') + scale_color_manual(values =  DonorPal ) ) |
     FeaturePlot(Tcells, c('CD4','PDCD1','FOXP3','CXCR5')))

#Feature plots
BasicT <- c('CD3E', 'CD8B', 'CD4', 'TRAC', 'TRDC')
UnconT <-c('KLRB1', "TRAV1-2", "ZBTB16", "ZNF683", "GNLY", "TYROBP")
SubsetsT <-c('CCR7', 'FOXP3', "TBX21", "GATA3", 'KLRG1') 
Feats <-c('CD8B', 'CD4', 'CCR7', 'SELL',"PDCD1", 'CXCR5', 'FOXP3', "KLRG1", 
          "PRF1","CD27",'KLRB1',"KLF2",'IL7R','CD44','LEF1','TCF7')
FPlots <- ViridisFeatures(Tcells, feature_args = Feats, 0.7, option = "D")
FP_wrap <- wrap_plots( FPlots[[1]] )
VP_wrap <- VlnPlot(Tcells, c('PDCD1','CXCR5'))

clus <- Tcells@meta.data %>%
  select(ident) %>%
  group_by_all() %>%
  tally()

#new idents
new.cluster.ids <- c("CD4_naive.C1",#0
                     "CD4_naive.C2",#1
                     "CD4_naive.C3",#2
                     "CD8_naive.C1",#3
                     "CD8_eff/mem.C1", #4
                     "CD4_eff/mem.C1", #5
                     "CD4_eff/mem.C2", #6
                     "CD8_naive.C2", #7
                     "CD4_eff/mem.C3", #8
                     "Treg", #9
                     "CD4_naive.C4", #10
                     "CD4_eff/mem.C4", #11
                     "CD4_naive.C5", #12
                     "CD4_naive.C6", #13
                     "CD8_eff/mem.C2", #14
                     "CD8_naive.C3", #15
                     "Tfh", #16
                     "CD8_naive.C4", #17
                     "MAIT/NKT", #18
                     "CD4_naive.C7", #19
                     "CD8_eff/mem.C3", #20
                     "cycling" #21
                     )

#Tcells <- SetIdent(Tcells, value = 'orig.ident')
Tcells@meta.data$orig.ident <- Tcells@active.ident       
names(new.cluster.ids) <- as.character(seq(new.cluster.ids) - 1)
Tcells <- RenameIdents(Tcells, new.cluster.ids)
Tcells@meta.data$ident <- Tcells@active.ident  

# group by T cell types
naives <- new.cluster.ids[which(new.cluster.ids %like% 'naive')] 
eff_mem <- new.cluster.ids[which(new.cluster.ids %like% 'eff/mem' | 
                                   new.cluster.ids %like% 'Treg' | 
                                   new.cluster.ids %like% 'cycling')] 
Tfh <- new.cluster.ids[which(new.cluster.ids %like% 'Tfh')]       
MAIT <- new.cluster.ids[which(new.cluster.ids %like% 'MAIT')]  

#assign new column 
Tcells@meta.data$Tcell_type <- case_when(
  Tcells@meta.data$ident %in% naives ~ "naive",
  Tcells@meta.data$ident %in% eff_mem ~ "effector/memory",
  Tcells@meta.data$ident %in% Tfh ~ "Tfh",
  Tcells@meta.data$ident %in% MAIT ~ "MAIT/NKT"
)

#reorder the ident levels
clusterorder <- c( naives, MAIT, eff_mem, Tfh) %>% as.character(.)
Tcells@active.ident <- factor(Tcells@active.ident, levels = clusterorder )
Tcells@meta.data$ident <- Tcells@active.ident  
Tcells@meta.data$Tcell_type <-  factor( Tcells@meta.data$Tcell_type, levels =  c('naive','MAIT/NKT','effector/memory','Tfh') )

# Find T cell cluster markers ====

TMarkers <- FindAllMarkers(Tcells )

TMarkers10 <- TMarkers %>% 
  group_by(cluster) %>% 
  slice_max( order_by = avg_log2FC , n =5) %>%
  pull(gene) %>%
  as.character(.)

Tcell_DEGS <- DoHeatmap(Tcells, TMarkers10 )
write_tsv(TMarkers, 'outs/TwoYear_Tcell_markers.tsv')


#save updated T cell object ----
saveRDS(Tcells , Tcells_path )
