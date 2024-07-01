
source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
library(future)
library(harmony)
plan(multisession, workers = 6)
session_mem <- 1000 * 1024^2
options(future.globals.maxSize= session_mem)

#read data ===

Tcells <- readRDS(Set2_integrated_Tcells_path)


old_TFH <- readRDS(Tfh_lineages_path)


#subset CD4 clusters ====


VlnPlot(Tcells, 'CD4')
CD4_cells <- WhichCells(Tcells,idents = c(0:3,6:11,14,16,18:20))
TwoYear_keeps <- paste0('TwoYear_', colnames(old_TFH))
in_cells <- union(CD4_cells,TwoYear_keeps )


CD4_obj <- subset(Tcells, cells =  in_cells)


# rerun harmony ====
CD4_obj <- CD4_obj %>% 
  FindVariableFeatures() 

pc_genes <- VariableFeatures(CD4_obj)[
  -which(VariableFeatures(CD4_obj) %in% IgTcr$gene)]

CD4_obj <- RunPCA(CD4_obj, features = pc_genes, npcs = 50) 


CD4_obj_harmonized <- RunHarmony(CD4_obj, 
                                     assay.use = "SCT", 
                                     group.by.vars = c('dataset','donor','time','tissue'))


CD4_obj_harmonized <- RunUMAP(CD4_obj_harmonized, dims=1:50, reduction = "harmony")
CD4_obj_harmonized <- FindNeighbors(CD4_obj_harmonized, reduction = "harmony") %>% 
  FindClusters(., resolution = 1.2)
CD4_obj_harmonized@meta.data$ident <- CD4_obj_harmonized@active.ident



# remove CD8s and cycling====
VlnPlot(CD4_obj_harmonized, 'CD8A')
CD4_cells2 <- CD4_obj_harmonized@meta.data %>% 
  filter(ident %notin% c(11,32)) %>%
  pull(barcode)

noCD8 <- FetchData(CD4_obj_harmonized,  c('CD8A','CD8B')) %>% 
  filter(CD8A == 0 & CD8B ==0) %>% 
  rownames()


CD4_cells2 <- CD4_cells2[which(CD4_cells2 %in% noCD8)]


CD4_cells2 <- union(CD4_cells2, TwoYear_keeps)
CD4_obj2 <- subset(CD4_obj_harmonized, cells = CD4_cells2 )
CD4_obj2 <- subset(CD4_obj2, subset = UMAP_2 < 10)

# rerun harmony ====

CD4_obj2 <- CD4_obj2 %>% 
  FindVariableFeatures() 

pc_genes <- VariableFeatures(CD4_obj2)[
  -which(VariableFeatures(CD4_obj2) %in% IgTcr$gene)]

CD4_obj2 <- RunPCA(CD4_obj2, features = pc_genes, npcs = 50) 


CD4_obj2_harmonized2 <- RunHarmony(CD4_obj2, 
                                      assay.use = "SCT", 
                                      group.by.vars = c('dataset','donor','time','tissue'))


# embedding and clustering ====

CD4_obj2_harmonized2 <- RunUMAP(CD4_obj2_harmonized2, dims=1:50, reduction = "harmony")
CD4_obj2_harmonized2 <- FindNeighbors(CD4_obj2_harmonized2, reduction = "harmony") 
CD4_obj2_harmonized2 <- FindClusters(CD4_obj2_harmonized2, resolution = 0.8)
CD4_obj2_harmonized2@meta.data$ident <-CD4_obj2_harmonized2@active.ident


Feats <- c('CD4', 'CCR7', 'SELL', "KLRG1", 'KLRB1', 'CXCR3',
           "FOXP3", 'ICOS','TOX2','CXCR5', "PDCD1", 'IL21')
FPs_panel <- FeaturePlot_scCustom(CD4_obj2_harmonized2, Feats, num_columns = 6, order = F) & NoLegend() & NoAxes()
DimPlot(CD4_obj2_harmonized2, label = T)


#Find DEGS ====

CD4_marks <- FindAllMarkers(CD4_obj2_harmonized2)
write.csv(CD4_marks, './10x/outs/integrated_CD4_markers.csv')


#save ====
saveRDS(CD4_obj2_harmonized2 , Set2_integrated_CD4_path )


CD4_obj2_harmonized2@meta.data$ident <- CD4_obj2_harmonized2@active.ident
nclust <- length(unique(CD4_obj2_harmonized2@meta.data$ident))

cpal <- append(glasbey() ,alphabet())[1:nclust]
names(cpal) <- seq(length(unique(CD4_obj2_harmonized2@meta.data$ident)))-1


( DimPlot(CD4_obj2_harmonized2, label = T,label.size = 6) + scale_color_manual(values = cpal)  | 
    DimPlot(CD4_obj2_harmonized2, group.by = 'donor') + scale_color_manual(values = DonorPal)) /
  (DimPlot(CD4_obj2_harmonized2, group.by = 'tissue') + scale_color_manual(values = TissuePal) |
     DimPlot(CD4_obj2_harmonized2, group.by = 'Tcell_type'))



freq <- CD4_obj2_harmonized2@meta.data %>%
  select(ident, Tcell_type) %>%
  group_by(ident) %>%
  add_tally(name = 'n_ident') %>%
  ungroup() %>%
  group_by(ident, Tcell_type) %>%
  add_tally(name = 'n_ident_set') %>%
  distinct_all() %>%
  mutate(freq = n_ident_set/ n_ident) %>%
  ggplot(.,aes(ident, freq, fill = Tcell_type)) + geom_col()
(FeaturePlot(CD4_obj2_harmonized2, c('PDCD1','ICOS','CXCR5','FOXP3')) | freq) 