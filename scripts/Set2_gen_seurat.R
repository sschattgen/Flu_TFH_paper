# Regenerate seurat object with the new aggreagate data from donors 5 and 4 
# this is version 2. This first might be too stringent. 

source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
# make Seurat object. This is the highest level with all cells ====

# generate Seurat object
gexFiles <- Read10X('./10x/objects/Set2/filtered_feature_bc_matrix/')
Both_sc <-  CreateSeuratObject( gexFiles$`Gene Expression` )
Both_sc[['HTO']] <- CreateAssayObject(counts = gexFiles$`Antibody Capture`)


# Add tcr data
Both_sc <- add_clones_to_metadata(Both_sc, parse_conga_clones_for_seurat('./10x/objects/Set2_TCR_clones.tsv'))
Both_sc@meta.data$subject <- NULL

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


Both_sc <- NormalizeData(Both_sc, assay = 'HTO', normalization.method = "CLR")


# demulitplex hastags
Both_sc <- MULTIseqDemux(Both_sc, assay = "HTO", autoThresh = T)


#dim reduction
pc_genes <- VariableFeatures(Both_sc)[-which(VariableFeatures(Both_sc) %in% IgTcr$gene)]
Both_sc <- RunPCA(Both_sc, features = pc_genes, verbose = FALSE) %>% 
  FindNeighbors(.) %>% 
  FindClusters(., resolution = 1, verbose = FALSE )%>% 
  RunUMAP(., dims=1:50)


Both_sc@meta.data$barcode <- colnames(Both_sc)
Both_sc@meta.data$suffix <- str_split(Both_sc@meta.data$barcode, pattern = '-', simplify = T)[,2] %>% as.integer()
Both_sc@meta.data$ident <- Both_sc@active.ident

# Hashing results results
table(Both_sc@meta.data$MULTI_ID)
table(Both_sc@meta.data$tissue)

# fill in meta.data based on HTOs and sample info ----
Both_sc@meta.data$tissue <- ifelse(Both_sc@meta.data$suffix %in% c(1,2), 'PBMC', 'FNA')

Both_sc@meta.data$year <- 1

Both_sc@meta.data$day <- case_when(
  Both_sc@meta.data$suffix %in% c(1,2) & Both_sc@meta.data$MULTI_ID =='Hu.Hashtag-9' ~ 0,
  Both_sc@meta.data$suffix %in% c(1,2) & Both_sc@meta.data$MULTI_ID =='Hu.Hashtag-10' ~ 5,
  Both_sc@meta.data$tissue == 'FNA' & Both_sc@meta.data$MULTI_ID =='Hu.Hashtag-9' ~ 26,
  Both_sc@meta.data$tissue == 'FNA' & Both_sc@meta.data$MULTI_ID =='Hu.Hashtag-10' ~ 12,
) 

Both_sc@meta.data$donor <- case_when(
  Both_sc@meta.data$suffix == 1  ~ '321-07',
  Both_sc@meta.data$suffix == 2  ~ '321-08',
  Both_sc@meta.data$suffix == 3 & Both_sc@meta.data$MULTI_ID =='Hu.Hashtag-9' ~ '321-07',
  Both_sc@meta.data$suffix == 3 & Both_sc@meta.data$MULTI_ID =='Hu.Hashtag-10' ~ '321-08',
)
backup <- Both_sc


# cleanup and factorize metadata ====
Both_sc@meta.data$time <- ifelse(!is.na(Both_sc@meta.data$day), 
                                 paste0('y',Both_sc@meta.data$year,"_d", Both_sc@meta.data$day), NA )
Both_sc@meta.data$time <- factor(Both_sc@meta.data$time , levels = time_levels)
Both_sc@meta.data$day <- factor(Both_sc@meta.data$day , levels = names(TimePal))
Both_sc@meta.data$year <- factor(Both_sc@meta.data$year , levels = names(YearPal))
Both_sc@meta.data$donor <- factor(Both_sc@meta.data$donor , levels = names(DonorPal))
Both_sc@meta.data$tissue <- factor(Both_sc@meta.data$tissue , levels = names(TissuePal))


#algo to correct metadata for clone sizes N>1 if hto was off ====
# 
# correct_clone <- function(df){  
#   cdonor <- names(sort(table(df$donor)))[5]
#   
#   if (nrow(df) > 1){
#     df %>%
#       mutate(donor = cdonor) %>%
#       mutate(day = case_when(
#         donor == '321-08' & tissue == 'FNA' ~ 12,
#         donor == '321-07' & tissue == 'FNA' ~ 26,
#         donor == '321-07' & tissue == 'PBMC'& day == 0 ~ 0,
#         donor == '321-07' & tissue == 'PBMC'& day == 5 ~ 5,
#         donor == '321-08' & tissue == 'PBMC'& day == 0 ~ 0,
#         donor == '321-08' & tissue == 'PBMC'& day == 5 ~ 5
#       ) %>% factor())
#   } else if (nrow(df) == 1 | all(is.na(df$clone_id))) {
#     df
#   }
# 
# }
# 
# clones <- Both_sc@meta.data %>%
#   group_by(clone_id) %>%
#   group_split()
# 
# corrected_clones <- furrr::future_map(clones, ~correct_clone(.)) 
# 
# fixed_clones <- do.call(rbind, clones) %>%
#   select(barcode, donor, tissue, day) %>% 
#   as.data.frame()
# 
# rownames(fixed_clones) <- fixed_clones$barcode
# 
# 
# 
# 
# 
# Both_sc <- AddMetaData(Both_sc, select(fixed_clones, -barcode))



# subset out the T cells ====

# identify the T cell clusters
FP <- FeaturePlot(Both_sc, c('CD3E','CD8A','CD4','LYZ','CD19','KLRK1','MKI67')) /
  DimPlot(Both_sc, label = T) /
  DimPlot(Both_sc, split.by = 'MULTI_ID', label = T) 

# cutting to the chase and keep only the cells with a paired TCR since there are some many
# 10X data was further parsed with make_10x_clone_file in conga. See ./repdata/10x/new/JCC280_make_10x_clone.py
# subsetting out on cluster and by matched TCR info
Tcells <- subset(Both_sc, subset = UMAP_1 < 5 )
DimPlot(Tcells)
ncol(Tcells)

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
Tcells <- RunUMAP(object = Tcells , dims = 1:50)

# stash the idents and umis in metadata for later
Tcells@meta.data$ident <- Tcells@active.ident

DimPlot(Tcells)


## Save the full and T cell objects ==== 
saveRDS(Both_sc , Set2_all_cells_path )
saveRDS(Tcells , Set2_Tcells_path )


#souporcell ----

geno_clusters <- read.delim('../../cluster_jobs/ss36_TFH/souporcell/clusters.tsv') %>%
  mutate(barcode = gsub('-1','-3', barcode)) %>%
  filter(barcode %in% colnames(Both_sc))

fna_sample <- Both_sc@meta.data %>% 
  filter(barcode %in%geno_clusters$barcode) %>%
  select(barcode, donor, MULTI_ID)


joint <- left_join(geno_clusters, fna_sample)

doublets <- joint %>% filter(status  == 'unassigned' & MULTI_ID == 'Doublet') %>% pull(barcode)

ggplot(joint, aes(cluster0,cluster1, color = donor, shape = assignment)) + geom_point() + facet_wrap(~donor)



ggplot(joint, aes(assignment,  fill = donor)) + geom_col() + facet_wrap(~donor)

DimPlot(Tcells, cells.highlight = doublets)



VlnPlot( Both_sc, features = c('percent.mt', "nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size = 0,)


doublets <- Both_sc@meta.data %>% filter( MULTI_ID == 'Doublet') %>% pull(barcode) 
table(Both_sc@meta.data$MULTI_ID)