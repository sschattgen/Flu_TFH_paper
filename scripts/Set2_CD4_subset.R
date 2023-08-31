
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
CD4_cells <- WhichCells(Tcells,idents = c(0:3,5:7,9,10,12,14,16:21))
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
                                     group.by.vars = c("dataset","donor"))


CD4_obj_harmonized <- RunUMAP(CD4_obj_harmonized, dims=1:50, reduction = "harmony")
CD4_obj_harmonized <- FindNeighbors(CD4_obj_harmonized, reduction = "harmony") %>% 
  FindClusters(., resolution = 1.2)
CD4_obj_harmonized@meta.data$ident <- CD4_obj_harmonized@active.ident



# remove CD8s and cycling====

CD4_cells2 <- CD4_obj_harmonized@meta.data %>% 
  filter(ident %notin% c(9,17,20,21,26,27)) %>%
  pull(barcode)

noCD8 <- FetchData(CD4_obj_harmonized,  c('CD8A','CD8B')) %>% 
  filter(CD8A == 0 & CD8B ==0) %>% 
  rownames()


CD4_cells2 <- CD4_cells2[which(CD4_cells2 %in% noCD8)]


CD4_cells2 <- union(CD4_cells2, TwoYear_keeps)

CD4_obj2 <- subset(CD4_obj_harmonized, cells = CD4_cells2 )
CD4_obj2 <- subset(CD4_obj2, subset = UMAP_2 > -7)

# rerun harmony ====

CD4_obj2 <- CD4_obj2 %>% 
  FindVariableFeatures() 

pc_genes <- VariableFeatures(CD4_obj2)[
  -which(VariableFeatures(CD4_obj2) %in% IgTcr$gene)]

CD4_obj2 <- RunPCA(CD4_obj2, features = pc_genes, npcs = 50) 


CD4_obj2_harmonized2 <- RunHarmony(CD4_obj2, 
                                      assay.use = "SCT", 
                                      group.by.vars = c("dataset","donor"))


# embedding and clustering ====

CD4_obj2_harmonized2 <- RunUMAP(CD4_obj2_harmonized2, dims=1:50, reduction = "harmony")
CD4_obj2_harmonized2 <- FindNeighbors(CD4_obj2_harmonized2, reduction = "harmony") %>% 
  FindClusters(.,resolution = 1.2)
CD4_obj2_harmonized2@meta.data$ident <-CD4_obj2_harmonized2@active.ident

# remove CD8s ====


CD4_obj3 <- subset(CD4_obj2_harmonized2, subset = UMAP_2 < 8)


# rerun harmony ====

CD4_obj3 <- CD4_obj3 %>% 
  FindVariableFeatures() 

pc_genes <- VariableFeatures(CD4_obj3)[
  -which(VariableFeatures(CD4_obj3) %in% IgTcr$gene)]

CD4_obj3 <- RunPCA(CD4_obj3, features = pc_genes, npcs = 50) 


CD4_obj3_harmonized <- RunHarmony(CD4_obj3, 
                                   assay.use = "SCT", 
                                   group.by.vars = c("dataset","donor"))


# embedding and clustering ====

CD4_obj3_harmonized <- RunUMAP(CD4_obj3_harmonized, dims=1:50, reduction = "harmony")
CD4_obj3_harmonized <- FindNeighbors(CD4_obj3_harmonized, reduction = "harmony") %>% 
  FindClusters(.)
CD4_obj3_harmonized@meta.data$ident <-CD4_obj3_harmonized@active.ident



#Find DEGS ====

CD4_marks <- FindAllMarkers(CD4_obj3_harmonized)
write.csv(CD4_marks, './10x/outs/integrated_CD4_markers.csv')


#save ====
saveRDS(CD4_obj3_harmonized , Set2_integrated_CD4_path )


CD4_obj3_harmonized@meta.data$ident <- CD4_obj3_harmonized@active.ident
nclust <- length(unique(CD4_obj3_harmonized@meta.data$ident))

cpal <- append(glasbey() ,alphabet())[1:nclust]
names(cpal) <- seq(length(unique(CD4_obj3_harmonized@meta.data$ident)))-1


( DimPlot(CD4_obj3_harmonized, label = T,label.size = 6) + scale_color_manual(values = cpal)  | 
    DimPlot(CD4_obj3_harmonized, group.by = 'donor') + scale_color_manual(values = DonorPal)) /
  (DimPlot(CD4_obj3_harmonized, group.by = 'tissue') + scale_color_manual(values = TissuePal) |
     DimPlot(CD4_obj3_harmonized, group.by = 'Tcell_type'))



freq <- CD4_obj3_harmonized@meta.data %>%
  select(ident, Tcell_type) %>%
  group_by(ident) %>%
  add_tally(name = 'n_ident') %>%
  ungroup() %>%
  group_by(ident, Tcell_type) %>%
  add_tally(name = 'n_ident_set') %>%
  distinct_all() %>%
  mutate(freq = n_ident_set/ n_ident) %>%
  ggplot(.,aes(ident, freq, fill = Tcell_type)) + geom_col()
(FeaturePlot(CD4_obj3_harmonized, c('PDCD1','ICOS','CXCR5','FOXP3')) | freq) 