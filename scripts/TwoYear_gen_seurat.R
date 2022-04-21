# Regenerate seurat object with the new aggreagate data from donors 5 and 4 
# this is version 2. This first might be too stringent. 

source('scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

# make Seurat object. This is the highest level with all cells ====

# generate Seurat object
gexFiles <- Read10X('data/JCC280_Year2Aggreg_CR6/filtered_feature_bc_matrix/')
Both_sc <-  CreateSeuratObject( gexFiles )

#append metadata for library, day, tissue
Both_sc@meta.data$barcode <- colnames(Both_sc)
Both_sc@meta.data$suffix <- str_split(Both_sc@meta.data$barcode, pattern = '-', simplify = T)[,2]

md <- read.csv('data/library_tissue_map.csv')
colnames(md)[1:2] <- c('library_id', 'suffix')
df <- FetchData(Both_sc, c('suffix', 'barcode'))
df$suffix <- as.integer(df$suffix)
new_meta_df <- left_join(df, md)

for ( i in c('library_id', 'tissue', 'day', 'donor', 'year')){
  new_md_V <- new_meta_df[[i]]
  names(new_md_V)<- df[['barcode']]
  Both_sc <- AddMetaData(Both_sc, metadata = new_md_V, col.name = i)
}

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

Both_sc <- NormalizeData( Both_sc)
Both_sc <- FindVariableFeatures(object = Both_sc)
Both_sc_var_genes <- VariableFeatures(Both_sc)
# I didnt regress for cell cycle here since isn't necessary to subset
# we can do that with the T cell subset which is much smaller. 
Both_sc <- ScaleData(object = Both_sc, features = Both_sc_var_genes )
Both_sc <- RunPCA(object = Both_sc, features = Both_sc_var_genes )
Both_sc <- FindNeighbors(object = Both_sc)
Both_sc <- FindClusters(object = Both_sc)
Both_sc <- RunUMAP(object = Both_sc , dims = 1:20)
Both_sc@meta.data$ident <- Both_sc@active.ident

# prepare TCR data and map to Both_sc ====

# read metadata
tcrmeta <- read.csv('data/tcr_metadata.csv', stringsAsFactors = F)

#read in mappings
#barcode suffixes in the map_files.list need to be adjusted to match each corresponding library in the aggreagated gex dataset to match metadata
map_files <- 'data/TwoYear_TCR_clones.tsv.barcode_mapping.tsv'
map_files.list <- map( map_files, ~ read.delim( . , stringsAsFactors = F) )
names(map_files.list) <- 'x' 

#read in clone files and add subjects to df
clone_files <- 'data/TwoYear_TCR_clones.tsv'
clone_files.list <- map( clone_files, ~ read.delim( . , stringsAsFactors = F) )
names(clone_files.list) <- 'x'

#let's map in the clone_id into the metadata by the barcode to guide the rest of the metadata input
#split the concatenated barcodes then melt together

melt_barcode <- function(x){
  samp <- x
  df <- data.frame(samp$clone_id,
                   str_split(samp$barcodes, pattern = ",", simplify = T),
                   stringsAsFactors = F)
  colnames(df)[1] <- "clone_id"
  df <- df %>% gather(2:ncol(df), value = barcode, key = col) %>% select(-col) %>% filter(barcode != "")
  return(df)
}
map_files.list <- lapply(map_files.list, melt_barcode)

#left join each pair of map, clone files, and metadata
clone_book <- left_join( map_files.list[[1]] , clone_files.list[[1]] , by = "clone_id")
clone_book$suffix <- str_split( clone_book$barcode , '-', simplify = T )[,2]
tcrmeta$suffix <- as.character( tcrmeta$suffix )
#tcrmeta <- tcrmeta %>% select(suffix, donor)
#tcrmeta$Day <- as.character( tcrmeta$Day )
clone_book <- left_join( clone_book, tcrmeta, 'suffix')

# need to update the subject to match Donor then rename clone_ids to match within each donor
#clone_book$subject <- clone_book$Donor

clone_subjs <- clone_book %>% 
  mutate(strict.match.cdr = paste(va_gene, cdr3a, ja_gene, 
                                  vb_gene, cdr3b, jb_gene,
                                  sep = ";") , 
         orig.clone_id = clone_id ) %>%
  select(-clone_id) %>% 
  group_by(donor) %>%
  group_split()

names(clone_subjs) <- c('321-04', '321-05')

out_clones <- list()
for (i in seq_along(clone_subjs)){
  
  uclones <- data.frame(strict.match.cdr = unique(clone_subjs[[i]]$strict.match.cdr))
  
  sub <- names(clone_subjs)[i]
  
  uclones.df <- uclones %>% 
    mutate(clone_id = paste("clone", rep( 1:length( uclones$strict.match.cdr )), sub, sep = '_') )
  
  out_clones[[i]] <- left_join( clone_subjs[[i]] , uclones.df, by = 'strict.match.cdr') %>% distinct()
  
}
rough_clones_df <- do.call(rbind, out_clones)

# map back to the full object

keeps <- c('clone_id', 'va_gene', 'ja_gene', 'vb_gene', 'jb_gene', 
           'cdr3a', 'cdr3b', 'cdr3a_nucseq', 'cdr3b_nucseq', 'strict.match.cdr', 
           'orig.clone_id') #select the values to be added into metadata

for (i in keeps) {
  value <- rough_clones_df %>% pull(i) 
  names(value) <- rough_clones_df$barcode
  Both_sc  <- AddMetaData(Both_sc , metadata = value, col.name =  i )
} #loops through keep to add in each value into the metadata by barcode

# Add BCR data to big object ====
# read metadata
bcrmeta <- read.csv('data/bcr_metadata.csv', stringsAsFactors = F)

#read in mappings
#barcode suffixes in the map_files.list need to be adjusted to match each corresponding library in the aggreagated gex dataset to match metadata
map_files <- 'data/TwoYear_BCR_clones.tsv.barcode_mapping.tsv'
map_files.list <- map( map_files, ~ read.delim( . , stringsAsFactors = F) )
names(map_files.list) <- 'x' 

#read in clone files and add subjects to df
clone_files <- 'data/TwoYear_BCR_clones.tsv'
clone_files.list <- map( clone_files, ~ read.delim( . , stringsAsFactors = F) )
names(clone_files.list) <- 'x'

# melt barcodes
map_files.list <- lapply(map_files.list, melt_barcode)

#left join each pair of map, clone files, and metadata
clone_book <- left_join( map_files.list[[1]] , clone_files.list[[1]] , by = "clone_id")
clone_book$suffix <- str_split( clone_book$barcode , '-', simplify = T )[,2]
bcrmeta$suffix <- as.character( bcrmeta$suffix )
bcrmeta <- bcrmeta %>% select(suffix, donor)
#tcrmeta$Day <- as.character( tcrmeta$Day )
clone_book <- left_join( clone_book, bcrmeta, 'suffix')

# need to update the subject to match Donor then rename clone_ids to match within each donor
#clone_book$subject <- clone_book$Donor

clone_subjs <- clone_book %>% 
  mutate(strict.match.cdr = paste(va_gene, cdr3a, ja_gene, 
                                  vb_gene, cdr3b, jb_gene,
                                  sep = ";") , 
         orig.clone_id = clone_id ) %>%
  select(-clone_id) %>% 
  group_by(donor) %>%
  group_split()

names(clone_subjs) <- c('321-04', '321-05')

out_clones <- list()
for (i in seq_along(clone_subjs)){
  
  uclones <- data.frame(strict.match.cdr = unique(clone_subjs[[i]]$strict.match.cdr))
  
  sub <- names(clone_subjs)[i]
  
  uclones.df <- uclones %>% 
    mutate(clone_id = paste("clone", rep( 1:length( uclones$strict.match.cdr )), sub, sep = '_') )
  
  out_clones[[i]] <- left_join( clone_subjs[[i]] , uclones.df, by = 'strict.match.cdr') %>% distinct()
  
}
rough_clones_df <- do.call(rbind, out_clones)

# map back to the full object
colnames(rough_clones_df)[c(4:15,23:25)] <- c('vl_gene','jl_gene','vl2_gene','jl2_gene', 'vh_gene', 'jh_gene',
                                     'cdr3l', 'cdr3l_nucseq', 'cdr3l2', 'cdr3l2_nucseq', 'cdr3h', 'cdr3h_nucseq',
                                     "bcr_strict.match.cdr" ,"bcr_orig.clone_id", "bcr_clone_id")

keeps <- c('bcr_clone_id', 'vl_gene', 'jl_gene', 'vh_gene', 'jh_gene', 
           'cdr3l', 'cdr3h', 'cdr3l_nucseq', 'cdr3h_nucseq', 'bcr_strict.match.cdr', 
           'bcr_orig.clone_id') #select the values to be added into metadata

for (i in keeps) {
  value <- rough_clones_df %>% pull(i) 
  names(value) <- rough_clones_df$barcode
  Both_sc  <- AddMetaData(Both_sc , metadata = value, col.name = i)
} #loops through keep to add in each value into the metadata by barcode

# cleanup and factorize metadata ====
Both_sc@meta.data$bcr_orig.clone_id <- gsub('bcr_bcr_', 'bcr_', Both_sc@meta.data$bcr_orig.clone_id)
Both_sc@meta.data$orig.clone_id <- gsub('tcr_tcr_', 'tcr_', Both_sc@meta.data$orig.clone_id)
Both_sc@meta.data$strict.match.cdr <- gsub('[*]..', '', Both_sc@meta.data$strict.match.cdr)
Both_sc@meta.data$bcr_strict.match.cdr <- gsub('[*]..', '', Both_sc@meta.data$bcr_strict.match.cdr)

Both_sc@meta.data$time <- paste0('y',Both_sc@meta.data$year,"_d", Both_sc@meta.data$day )
Both_sc@meta.data$time <- factor(Both_sc@meta.data$time , levels = time_levels)
Both_sc@meta.data$day <- factor(Both_sc@meta.data$day , levels = names(TimePal))
Both_sc@meta.data$year <- factor(Both_sc@meta.data$year , levels = names(YearPal))
Both_sc@meta.data$donor <- factor(Both_sc@meta.data$donor , levels = names(DonorPal))
Both_sc@meta.data$tissue <- factor(Both_sc@meta.data$tissue , levels = names(TissuePal))

Both_sc@meta.data$clone_info <- case_when( 
  !is.na(Both_sc@meta.data$bcr_clone_id) & !is.na(Both_sc@meta.data$clone_id) ~ 'both',
  is.na(Both_sc@meta.data$bcr_clone_id) & !is.na(Both_sc@meta.data$clone_id) ~ 'TCR',
  !is.na(Both_sc@meta.data$bcr_clone_id) & is.na(Both_sc@meta.data$clone_id) ~ 'BCR',
  is.na(Both_sc@meta.data$bcr_clone_id) & is.na(Both_sc@meta.data$clone_id) ~ 'none')

# subset out the T cells ====

# identify the T cell clusters
FP <- FeaturePlot(Both_sc, c('CD3E','CD19','TRDC','MKI67')) /
  DimPlot(Both_sc, reduction = 'umap', label = T)
VP <- VlnPlot(Both_sc,  c('CD3E','CD19','MKI67'))
lin_df <- FetchData(Both_sc, c('ident','CD3E','CD19','MKI67'))


ggplot(lin_df, aes(CD3E, CD19)) +
  geom_point() + 
  facet_wrap(~ident)
DimPlot(Both_sc, group.by = 'donor', split.by = 'year') + 
  scale_color_manual(values = DonorPal)
DimPlot(Both_sc, group.by = 'day', split.by = 'year')

Tcells_clust <- c(0,3,4,5,7,11,14,15,16,18,28)

# cutting to the chase and keep only the cells with a paired TCR since there are some many
# 10X data was further parsed with make_10x_clone_file in conga. See ./repdata/10x/new/JCC280_make_10x_clone.py
# subsetting out on cluster and by matched TCR info
Tcells <- subset(Both_sc, subset = ident %in% Tcells_clust & 
                   !is.na(Both_sc@meta.data$clone_id) &
                   is.na(Both_sc@meta.data$bcr_clone_id))
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
Tcells <- RunUMAP(object = Tcells , dims = 1:20)

# stash the idents and umis in metadata for later
Tcells@meta.data$ident <- Tcells@active.ident

## Save the full and T cell objects ==== 
saveRDS(Both_sc , all_cells_path )
saveRDS(Tcells , Tcells_path )


