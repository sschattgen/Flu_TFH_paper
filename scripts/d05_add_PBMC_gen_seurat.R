# Regenerate seurat object with the new aggreagate data from donors 5 and 4 
# this is version 2. This first might be too stringent. 

source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

# make Seurat object. This is the highest level with all cells ====

# generate Seurat object
gexFiles <- Read10X('./10x/objects/d05_add_PBMC/sample_filtered_feature_bc_matrix/')
Both_sc <-  CreateSeuratObject( gexFiles)


#append metadata for library, day, tissue
Both_sc@meta.data$barcode <- colnames(Both_sc)
Both_sc@meta.data$suffix <- str_split(Both_sc@meta.data$barcode, pattern = '-', simplify = T)[,2] %>% as.numeric()

# Add tcr data

Both_sc <- add_clones_to_metadata(Both_sc, parse_conga_clones_for_seurat('./10x/objects/d05_add_PBMC_TCR_clones.tsv'))
Both_sc@meta.data$subject <- '321-05'
Both_sc@meta.data$tissue <- 'PBMC'


# QC
Both_sc[["percent.mt"]] <- PercentageFeatureSet(Both_sc, pattern = "^MT-")
vln <- VlnPlot( Both_sc, features = c('percent.mt', "nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size = 0)
quantile(Both_sc@meta.data$percent.mt, c(.05, .50, .95)) #adjusted mito up to 10
quantile(Both_sc@meta.data$nFeature_RNA, c(.05, .50, .95)) #looks good
quantile(Both_sc@meta.data$nCount_RNA, c(.05, .50, .95)) #looks good

Both_sc <- subset(Both_sc, subset = 
                  percent.mt < 10 & 
                    nFeature_RNA < 4000 & nFeature_RNA > 300 & 
                    nCount_RNA < 40000 & nCount_RNA > 900 )

#norm and scale
Both_sc <- NormalizeData(Both_sc) %>% 
  FindVariableFeatures() %>% 
  ScaleData(., features = VariableFeatures(.))

#dim reduction
pc_genes <- VariableFeatures(Both_sc)[-which(VariableFeatures(Both_sc) %in% IgTcr$gene)]
Both_sc <- RunPCA(Both_sc, features = pc_genes, verbose = FALSE) %>% 
  FindNeighbors(.) %>% 
  FindClusters(., resolution = 1, verbose = FALSE )%>% 
  RunUMAP(., dims=1:50)

Both_sc@meta.data$barcode <- colnames(Both_sc)
Both_sc@meta.data$ident <- Both_sc@active.ident

# Hashing results results
table(Both_sc@meta.data$MULTI_ID)
table(Both_sc@meta.data$tissue)

# fill in meta.data based on HTOs and sample info ----
Both_sc@meta.data$year <- 1

Both_sc@meta.data$day <- 5
Both_sc@meta.data$donor <- '321-05'
# cleanup and factorize metadata ====
Both_sc@meta.data$time <- ifelse(!is.na(Both_sc@meta.data$day), 
                                 paste0('y',Both_sc@meta.data$year,"_d", Both_sc@meta.data$day), NA )
Both_sc@meta.data$time <- factor(Both_sc@meta.data$time , levels = time_levels)
Both_sc@meta.data$day <- factor(Both_sc@meta.data$day , levels = names(TimePal))
Both_sc@meta.data$year <- factor(Both_sc@meta.data$year , levels = names(YearPal))
Both_sc@meta.data$donor <- factor(Both_sc@meta.data$donor , levels = names(DonorPal))
Both_sc@meta.data$tissue <- factor(Both_sc@meta.data$tissue , levels = names(TissuePal))


# subset out the T cells ====

# identify the T cell clusters
FP <- FeaturePlot(Both_sc, c('CD3E','CD19','CD4','TOX2','CXCR3','CXCR4')) /
  DimPlot(Both_sc, label = T) 
# cutting to the chase and keep only the cells with a paired TCR since there are some many
# 10X data was further parsed with make_10x_clone_file in conga. See ./repdata/10x/new/JCC280_make_10x_clone.py
# subsetting out on cluster and by matched TCR info


## Save the full and T cell objects ==== 
saveRDS(Both_sc , d05_Tcells_path )



