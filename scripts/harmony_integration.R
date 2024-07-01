source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
library(future)
library(harmony)
options(future.globals.maxSize= 4194304000)
plan('multisession', workers = 4)
# import 3 objects of interest ====

TwoYear <- readRDS( Tcells_path )
Set2 <- readRDS(Set2_Tcells_path )
d05 <- readRDS(d05_Tcells_path )

# merge the raw counts together and ignore ====
all_Tcells <- merge(TwoYear, y = c(Set2, d05), 
                    add.cell.ids = c("TwoYear", "Set2", "d05"), 
                    project = "FluTFH")


all_Tcells@meta.data$dataset <- str_split(colnames(all_Tcells),'_', simplify = T)[,1]
all_Tcells@meta.data$barcode <- colnames(all_Tcells)

#normalize, scale, pca ====
all_Tcells <- NormalizeData(all_Tcells) %>% 
  FindVariableFeatures() %>% 
  ScaleData(., features = VariableFeatures(.))

all_Tcells <- SCTransform(all_Tcells)
DefaultAssay(all_Tcells) <- 'SCT'

pc_genes <- VariableFeatures(all_Tcells)[
  -which(VariableFeatures(all_Tcells) %in% IgTcr$gene)]

all_Tcells <- RunPCA(all_Tcells, features = pc_genes, npcs = 50) 

# harmony ====
# all_Tcells_harmonized <- RunHarmony(all_Tcells, 
#                          assay.use = "SCT", 
#                          group.by.vars = "dataset")
# 
# # embedding and clustering ====
# 
# all_Tcells_harmonized <- RunUMAP(all_Tcells_harmonized, dims=1:50, reduction = "harmony")
# all_Tcells_harmonized <- FindNeighbors(all_Tcells_harmonized, reduction = "harmony",) %>% 
#   FindClusters(., resolution = 1.2)
# 
# all_Tcells_harmonized@meta.data$ident <- all_Tcells_harmonized@active.ident
# 
# 
# 
# #save=====
# saveRDS(all_Tcells_harmonized , './10x/objects/intergrated_Tcells_harmony.rds' )



#version without harmonizing ====


all_Tcells <- RunUMAP(all_Tcells, dims=1:50) %>% 
  FindNeighbors(.) %>% 
  FindClusters(.)

all_Tcells@meta.data$ident <- all_Tcells@active.ident


saveRDS(all_Tcells , './10x/objects/intergrated_Tcells_SCTonly.rds' )



# harmonize with donor too ====
all_Tcells <- readRDS('./10x/objects/intergrated_Tcells_SCTonly.rds')

donor_cells <- all_Tcells@meta.data %>%
  filter(!is.na(donor)) %>%
  filter(MULTI_ID %notin% c('Negative','Doublet')) %>%
  pull(barcode)

all_Tcells_X <- subset(all_Tcells, cells = donor_cells)

group_vars <- c('dataset','donor','time','tissue') # 'time','tissue' are new

for(var in group_vars){
  print(any(is.na(all_Tcells_X@meta.data[[var]])))
  print(table(all_Tcells_X@meta.data[[var]]))
}


all_Tcells_harmonized_X <- RunHarmony(all_Tcells_X, reduction = 'pca',
                                    assay.use = "SCT", 
                                    group.by.vars = group_vars) 

all_Tcells_harmonized_X <- RunUMAP(all_Tcells_harmonized_X, dims=1:50, reduction = "harmony")
all_Tcells_harmonized_X <- FindNeighbors(all_Tcells_harmonized_X, reduction = "harmony")
all_Tcells_harmonized_X <-  FindClusters(all_Tcells_harmonized_X, resolution = 0.7)
all_Tcells_harmonized_X@meta.data$ident <- all_Tcells_harmonized_X@active.ident

saveRDS(all_Tcells_harmonized_X , './10x/objects/intergrated_Tcells_harmony_bydonor.rds' )

# DEGs ====

TMarkers <- FindAllMarkers(all_Tcells_harmonized_X)

write.csv(TMarkers, './10x/outs/integrated_markers.csv')

  
 
