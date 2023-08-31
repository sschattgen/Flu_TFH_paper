# BothDonor clonotype analysis 
# the alluvial plots of clonal lineage abundance for the picked clones is in here, too

source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

#
Tfh_obj <- readRDS(Set2_integrated_Tfh_path)
Tfh_obj@meta.data$day <- factor(Tfh_obj@meta.data$day, levels = names(TimePal))
Tfh_obj@meta.data$time_point <- factor(Tfh_obj@meta.data$time, levels = names(TimePal2))

Clonedf <- read.delim( clone_df_path , stringsAsFactors = F)

# import clone_df and Tfh obj ==== 
ccells <- Tfh_obj@meta.data %>%
  filter(!is.na(Tfh_clone_id)) %>%
  rownames()

pc_obj <- subset(Tfh_obj, cells = ccells) 

pc_obj@meta.data %>% filter(Tfh_clone_id =='Tfh_321-05_clone_3') %>% group_by(tissue) %>% tally()
DefaultAssay(pc_obj) <- 'RNA'
pc_obj <- SeuratReprocess(pc_obj)
pc_obj@meta.data$ident <- pc_obj@active.ident

pcid <- names(clonePal)[c(1,3,11,12)]
pc_obj@meta.data$flu_specific <- ifelse(pc_obj@meta.data$Tfh_clone_id %in% pcid, 'yes','unknown')

pc_markers <- FindAllMarkers(pc_obj)

pc_degs <- pc_markers %>%
  group_by(cluster)%>%
  slice_min(n = 10, order_by = p_val_adj)

saveRDS(pc_obj, './10x/objects/Set2_Donor_321-05_pickedClones.rds')
write_csv(pc_markers, './10x/outs/Set2_Donor_321-05_pickedClones_markers.csv')

# plotting
gex_clust_pal <- wes_palette("Darjeeling1")[1:3]
names(gex_clust_pal) <- seq(0,2)

flu_pal <- c("grey90","red")
names(flu_pal) <- c("unknown","yes")

pc_umap <- DimPlot(pc_obj, pt.size = 3) + 
  scale_color_manual(values = gex_clust_pal) +
  labs(title =  'GEX clusters') +
  CenterTitle() +
  theme(axis.line = element_blank())

pc_day_umap <- DimPlot(pc_obj, group.by = 'time', pt.size = 3) + 
  scale_color_manual(values = TimePal2) +
  theme(axis.line = element_blank())

pc_clone <- DimPlot(pc_obj, group.by = 'Tfh_clone_id', pt.size = 3) + 
  scale_color_manual(values = clonePal)+
  labs(title =  'Clone') +
  theme(axis.line = element_blank())

pc_ag_specific_umap <- DimPlot(pc_obj, group.by = 'flu_specific', pt.size = 3) +
  scale_color_manual(values = flu_pal)+
  labs(title =  'Flu specific') +
  theme(axis.line = element_blank())

pc_tissue <- DimPlot(pc_obj, group.by = 'tissue', pt.size = 3) +
  scale_color_manual(values = TissuePal)+
  labs(title =  'Tissue') +
  theme(axis.line = element_blank())

pc_tfh <- DimPlot(pc_obj, group.by = 'Tfh_type', pt.size = 3) +
  scale_color_jcolors('pal9') +
  labs(title =  'Tfh type') +
  theme(axis.line = element_blank())

pc_umap_panel <- pc_umap + 
  pc_day_umap + 
  pc_clone + 
  pc_ag_specific_umap + 
  pc_tissue + 
  pc_tfh

ggsave('./10x/outs/Set2_Donor_321-05_Tfh_pickedClone_umaps.png', plot = pc_umap_panel,
       height = 10, width = 15)

long_pca <- DimPlot(pc_obj, group.by = 'Tfh_clone_id', pt.size = 3, split.by = 'day', reduction = 'pca')+ 
  scale_color_manual(values = clonePal) +
  labs(title =  'Picked clonotypes') 

long_pca_tissue <- DimPlot(pc_obj, group.by = 'tissue', pt.size = 3, split.by = 'day', reduction = 'pca')+ 
  scale_color_manual(values = TissuePal) +
  labs(title =  'Picked clonotypes') 

ggsave('./10x/outs/Set2_Donor_321-05_Tfh_pickedClone_pca_time.png', plot = long_pca,
       height = 4, width = 18)



# fancier pca

pca_df <- FetchData(pc_obj, c('PC_1', 'PC_2','year','day', 'tissue', 'Tfh_clone_id'))

fancy_pca <- ggplot(pca_df, aes(PC_2, PC_1, color = tissue, fill = Tfh_clone_id )) + 
  geom_point(shape = 21, size = 2, stroke = 1.5, alpha = 0.8) +
  scale_fill_manual(values = clonePal) +
  scale_color_manual(values = TissuePal) +
  theme_minimal_grid() +
  facet_grid(year ~ day) 

ggsave('./10x/outs/Set2_Donor_321-05_Tfh_pickedClone_fancy_pca_time.pdf', 
       plot = fancy_pca,
       height = 5, width = 18, useDingbats = F)

#focus on just the PC1 scores
pca_df2 <- pca_df %>% 
  mutate(time = paste0('y',year,'_d', day)) %>% 
  arrange(year, day) %>%
  mutate(time = factor(time, levels = names(TimePal2))) %>%
  mutate(Tfh_clone_id = factor(Tfh_clone_id, levels = names(clonePal))) %>%
  mutate(tissue = as.character(tissue))

pca_df_mean_pc <- pca_df2 %>%
  group_by(time, Tfh_clone_id) %>%
  summarise(across(PC_1:PC_2, ~ mean(.x, na.rm = TRUE))) 

fancy_pca2 <- ggplot(pca_df2, aes(Tfh_clone_id, PC_1, fill = Tfh_clone_id )) + 
  geom_col(data = pca_df_mean_pc, alpha = 0.8, ) +
  geom_point(aes(color = tissue), shape = 21, size = 2, stroke = 1.5, alpha = 0.8) +
  scale_fill_manual(values = clonePal) +
  scale_color_manual(values = TissuePal) +
  theme_minimal_grid() +
  facet_grid(~time) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
ggsave('./10x/outs/Set2_Donor_321-05_Tfh_pickedClone_fancy_pca_time_PC1_only.pdf', 
       plot = fancy_pca2,
       height = 5, width = 18, useDingbats = F)


#fancy heatmaps ====

#color palettes for heatmap
ann_colors = list(
  'Clone ID' = clonePal,
  tissue = TissuePal,
  'Flu specific' = flu_pal,
  day = TimePal,
  year= YearPal,
  Tfh_type =TfhPal
)
# annotation metadata
sample_col <- FetchData(pc_obj, c('Tfh_clone_id', 'flu_specific' , 'year', 'day', 'tissue', 'time_point','Tfh_type'))%>%
  arrange( time_point , Tfh_clone_id)%>% 
  select(-time_point)
colnames(sample_col)[1:2] <- c('Clone ID', 'Flu specific')
sample_md <- sample_col %>% 
  select(`Clone ID`, `Flu specific`, Tfh_type, tissue, day, year)


gaps_col = sample_col %>% filter(year ==1) %>% nrow()
gaps_col_df <- sample_col %>% 
  group_by(year, day) %>% 
  group_split()  


gaps <- c()
for (i in seq_along(gaps_col_df)){
  df <-  gaps_col_df[[i]] %>% 
    mutate_if(is.factor, as.character)

  new_break = nrow(df)
  
  if (i == 1){
    
    gaps <- append(gaps, new_break)
  } else{
    new_break = new_break + gaps[i-1] 
    gaps <- append(gaps, new_break)
    
  }
  
  
}

# making two versions
# one with handpicked genes
picked_tfh_gene_set <-  c('CXCR5', 'PDCD1','CCR7','SELL',
                   'VIM', 'TMSB10', 'EMP3','TCF7',
                   'TOX2','GNG4', 'CXCL13', 'IL21',
                   'FOXB1','DUSP4', 'DDIT4', 'IL10' ,'MAP3K8')

# pull the scaled expression for those genes
scale_data <- pc_obj@assays$RNA@scale.data
picked_degs_scale_data <- scale_data[which(rownames(scale_data) %in% picked_tfh_gene_set), ] %>% 
  as.data.frame() %>%
  select(rownames(sample_col)) # reorder columns to match the df

#plot heatmap
pc_hm <- pheatmap(picked_degs_scale_data, 
                  color = viridis::viridis(50, option = "B"),
                  cluster_cols = FALSE,
                  annotation_col = sample_md, 
                  fontsize_row = 10, border_color= NA,
                  annotation_colors = ann_colors,
                  show_colnames  = FALSE, gaps_col= gaps_col) %>%
  as.ggplot(.)

ggsave('./10x/outs/Set2_Donor_321-05_Tfh_pickedClone_picked_marker_hm.pdf', plot = pc_hm,
       height = 8, width = 13, useDingbats = F)

# now let's use the head an tail loadings for PC1
slice_size = 20
PC_Tfh <- Loadings(pc_obj[['pca']])[,1:20]
PC1_Tfh <- PC_Tfh[,1] %>% sort()
pl_Tfh <- length(PC1_Tfh)
PC_load_Tfh <- c( names(PC1_Tfh)[1:slice_size] , names(PC1_Tfh)[(pl_Tfh-slice_size):pl_Tfh] )
tfh_gene_set <- PC_load_Tfh

degs_scale_data <- scale_data[tfh_gene_set, ] %>% 
  as.data.frame() %>%
  select(rownames(sample_col)) # reorder columns to match the df

pc_hm2 <- pheatmap(degs_scale_data, 
                  color = viridis::viridis(50, option = "B"),
                  cluster_cols = FALSE,
                  annotation_col = sample_md, 
                  fontsize_row = 10, border_color= NA,
                  annotation_colors = ann_colors,
                  show_colnames  = FALSE, gaps_col= gaps_col) %>%
  as.ggplot(.)

ggsave('./10x/outs/Set2_Donor_321-05_Tfh_pickedClone_PC_marker_hm.pdf', plot = pc_hm2,
       height = 8, width = 13, useDingbats = F)


#more complicated heatmaps faceted by lineage 
lineages <- sample_col %>%
  mutate(barcode = rownames(.)) %>%
  group_by(`Clone ID`) %>%
  group_split()

picked_tfh_gene_set2 <-  c('CXCR5', 'PDCD1', 'VIM', 'TMSB10', 
                           'EMP3','TCF7', 'TOX2','GNG4', 
                           'CXCL13', 'IL21','DUSP4', 'DDIT4', 
                           'IL10' ,'MAP3K8')

breaksList <- seq(-4, 4, by = 0.5)

picked_lineage_hm <- list()
PC_lineage_hm <- list()
for ( i in seq_along(lineages)){

  df <- lineages[[i]] %>% as.data.frame()
  colnames(df)[1:2] <- c('Clone ID', 'Flu specific')
  rownames(df) <- df$barcode
  df <- df %>% 
    select(`Clone ID`, `Flu specific`, Tfh_type, tissue, day, year)
  gaps_col <- df %>% filter(year == 1) %>% nrow()
  
  picked_degs_scale_data2 <- picked_degs_scale_data[picked_tfh_gene_set2, which(colnames(picked_degs_scale_data) %in% rownames(df))] %>% 
    as.data.frame() %>%
    select(rownames(df))
  

  #plot heatmap
  pc_hm <- pheatmap(picked_degs_scale_data2, 
                    color = viridis::viridis(length(breaksList), option = "B"),
                    breaks = breaksList,
                    cluster_cols = FALSE, cluster_rows = FALSE,
                    annotation_col = df, 
                    fontsize_row = 10, border_color= NA,
                    annotation_colors = ann_colors,
                    show_colnames  = FALSE, 
                    annotation_legend = FALSE,
                    main = unique(df$`Clone ID`), 
                    gaps_col= gaps_col) %>%
    as.ggplot(.)
  
  picked_lineage_hm <- list.append(picked_lineage_hm , pc_hm)
  
  
  degs_scale_data2 <- degs_scale_data[, which(colnames(degs_scale_data) %in% rownames(df))] %>% 
    as.data.frame() %>%
    select(rownames(df))
  
  pc_hm2 <- pheatmap(degs_scale_data2, 
                     color = viridis::viridis(length(breaksList), option = "B"),
                     breaks = breaksList,
                    cluster_cols = FALSE, cluster_rows = FALSE,
                    annotation_col = df, 
                    fontsize_row = 10, border_color= NA,
                    annotation_colors = ann_colors,
                    show_colnames  = FALSE, 
                    annotation_legend = FALSE,
                    main = unique(df$`Clone ID`), gaps_col= gaps_col) %>%
    as.ggplot(.)
  
  PC_lineage_hm <- list.append(PC_lineage_hm , pc_hm2)
  
  
}

lineage_hm_lay <- wrap_plots(picked_lineage_hm) 
ggsave('./10x/outs/Set2_Donor_321-05_Tfh_pickedClones_lineages_picked_markers_hm.pdf', plot = lineage_hm_lay,
       height = 12.5, width = 17, useDingbats = F)
PC_lineage_hm_lay <- wrap_plots(PC_lineage_hm) 
ggsave('./10x/outs/Set2_Donor_321-05_Tfh_pickedClones_lineages_PC_markers_hm.pdf', plot = PC_lineage_hm_lay,
       height = 12.5, width = 17, useDingbats = F)

# alluvial plot for all the pickedClones with 10x data ====

day_split <- Clonedf %>% 
  filter(donor == '321-05') %>%
  select(donor, tissue, time_point, Tfh_clone_id) %>%
  group_by_all() %>%
  tally() %>%
  ungroup(Tfh_clone_id) %>%
  group_split()

day_split_freq <- list()
for ( i in seq_along(day_split)){
  
   df <- day_split[[i]] %>%
    mutate(freq = n / sum(n))
  
   day_split_freq[[i]] <-df
}

merged_picked_clones <- do.call( rbind, day_split_freq) %>% 
  filter( !is.na(Tfh_clone_id)) %>%
  group_by(tissue) %>%
  group_split()
FNA_df <- merged_picked_clones[[1]]
PBMC_df <- merged_picked_clones[[2]]
PBMC_df$time_point <- factor(PBMC_df$time_point, levels = names(TimePal2))
FNA_df$time_point <- factor(FNA_df$time_point, levels = names(TimePal2))
FNA_df$blood <- ifelse(FNA_df$Tfh_clone_id %in% merged_picked_clones[[2]]$Tfh_clone_id, 'yes', 'no')
new_pc_df <- bind_rows(FNA_df, merged_picked_clones[[2]]) 
new_pc_df$time_point <- factor(new_pc_df$time_point, levels = names(TimePal2))
new_pc_df$tissue <- factor(new_pc_df$tissue, levels = c('PBMC','FNA'))

bloodpal <- c(TissuePal[2], 'black')
names(bloodpal) <- c('yes', 'no')

# alluvial plots
picked_pl_all <- ggplot(data = new_pc_df,
       aes(x = time_point, 
           y = n, 
           alluvium = Tfh_clone_id)) +
  geom_alluvium(aes(fill = Tfh_clone_id, color = blood),
                alpha = .75, 
                width = 1/2,
                decreasing = FALSE,
                size = 0.7) +
  theme_minimal()+
  labs(title = 'Cell numbers of clones picked for screening') +
  ylab('Number of cells') +
  scale_fill_manual(values = clonePal) +
  scale_color_manual(values = bloodpal) + 
  scale_y_continuous(limits = c(0,60), breaks = c(0,10,20,30,40,50,60))+
  facet_wrap(~tissue, ncol = 1)

ggsave('./10x/outs/Set2_BothDonors_321-05_pickedClone_clone_numbers_LN_PBMC.pdf', 
       plot = picked_pl_all,
       height = 12, width = 12, useDingbats = F)


picked_pl_FNA <- ggplot(data = filter(new_pc_df, tissue == 'FNA'),
                        aes(x = time_point, 
                            y = freq, 
                            alluvium = Tfh_clone_id)) +
  geom_alluvium(aes(fill = Tfh_clone_id, color = blood),
                alpha = .75, 
                width = 1/2,
                decreasing = FALSE,
                size = 0.7) +
  theme_minimal()+
  labs(title = 'Frequency of clones picked for cloning in FNA',
       subtitle = paste0(unique(new_pc_df$donor), 
                         " FNA")) +
  ylab('Relative frequency for tissue and day') +
  scale_fill_manual(values = clonePal) +
  scale_color_manual(values = bloodpal)

ggsave('./10x/outs/Set2_BothDonors_321-05_pickedClone_clone_freqs.pdf', 
       plot = picked_pl_FNA,
       height = 6, width = 12, useDingbats = F)

picked_pl_blood <- ggplot(data = PBMC_df,
                        aes(x = time_point, 
                            y = freq, 
                            alluvium = Tfh_clone_id)) +
  geom_alluvium(aes(fill = Tfh_clone_id),
                alpha = .75, 
                width = 1/2,
                decreasing = FALSE,
                size = 0.7) +
  theme_minimal()+
  labs(title = 'Frequency of clones picked for cloning in PBMC',
       subtitle = paste0(unique(PBMC_df$donor), 
                         " ", 
                         unique(PBMC_df$tissue))) +
  ylab('Relative frequency for tissue and day') +
  scale_fill_manual(values = clonePal) +
  scale_color_manual(values = bloodpal)

picked_pl_FNA_n <- ggplot(data = FNA_df,
                          aes(x = time_point, 
                              y = n, 
                              alluvium = Tfh_clone_id)) +
  geom_alluvium(aes(fill = Tfh_clone_id, color = blood),
                alpha = .75, 
                width = 1/2,
                decreasing = FALSE,
                size = 0.7) +
  theme_minimal()+
  labs(title = 'Cell numbers of clones picked for cloning in FNA',
       subtitle = paste0(unique(FNA_df$donor), " ", unique(FNA_df$tissue))) +
  ylab('Number for tissue and day') +
  scale_fill_manual(values = clonePal)+
  scale_color_manual(values = bloodpal)

ggsave('./10x/outs/Set2_BothDonors_321-05_pickedClone_clone_number.pdf', 
       plot = picked_pl_FNA_n,
       height = 6, width = 12, useDingbats = F)





hits <- c("Tfh_321-05_clone_1", "Tfh_321-05_clone_3",
          "Tfh_321-05_clone_11", "Tfh_321-05_clone_12")

FNA_df <-FNA_df %>%
  mutate(flu_specific = ifelse(Tfh_clone_id %in% hits, Tfh_clone_id, 'other'))

picked_pl_hits <- ggplot(data = FNA_df,
                        aes(x = time_point, 
                            y = freq, 
                            alluvium = Tfh_clone_id)) +
  geom_alluvium(aes(fill = flu_specific, color = blood),
                alpha = .75, 
                width = 1/2,
                decreasing = FALSE,
                size = 0.7) +
  theme_minimal()+
  labs(title = 'Frequency of confirmed Flu-specific clones',
       subtitle = paste0(unique(FNA_df$donor), 
                         " ", 
                         unique(FNA_df$tissue))) +
  ylab('Relative frequency for tissue and day') +
  scale_fill_manual(values = clonePal) +
  scale_color_manual(values = bloodpal)

ggsave('./10x/outs/Set2_BothDonors_321-05_pickedClone_flu_specific_clone_freqs.pdf', 
       plot = picked_pl_hits,
       height = 6, width = 12, useDingbats = F)

# alluvial plots with bulk data ====

Clonedf <- read.delim(clone_df_path, stringsAsFactors = F)
pickedClone_seqs <- Clonedf %>%
  filter(!is.na(Tfh_clone_id)) %>%
  select(Tfh_clone_id, va_gene, ja_gene, cdr3a,  vb_gene, jb_gene, cdr3b, match.cdr) %>%
  mutate(match_alpha = paste(va_gene, cdr3a, ja_gene, sep = ","),
         match_beta= paste(vb_gene, cdr3b, jb_gene, sep = ",")
  ) %>% #using vj genes + cdr3 to match
  distinct_all() 

# import bulk data

bulk_clones <- read.delim("./repdata/bulk/clone_tables/Ali_bulk_allDonor_clones.tsv", 
                          stringsAsFactors = F)
#split to alpha beta
bulk_d5_chains <- bulk_clones %>%
  filter(donor == '321-05') %>%
  mutate(day = gsub('d','', time)) %>%
  mutate(day = factor(day, levels = c(5,12,28,60))) %>%
  mutate(match = paste(v, cdr3aa, j, sep = ",")) %>%
  group_by(chain) %>%
  group_split()


alluv_plots <- list()
match_dfs <- list()
for (i in seq_along(bulk_d5_chains)){
  
  chain <- ifelse( i == 1, 'TCRa', 'TCRb')
  match_chain <- ifelse( chain == 'TCRa', 'match_alpha','match_beta')
  label_y_pos <- ifelse( chain == 'TCRa', -0.002,-0.004)
  
  
  #make match df
  patterns <- pickedClone_seqs[[match_chain]]
  d5_match <- bulk_d5_chains[[i]] %>%
    filter(match %in% patterns)
  colnames(d5_match)[17] <- match_chain
  
  d5_match <- left_join( d5_match, select(pickedClone_seqs, c(match_chain, Tfh_clone_id) ) )
  
  #plot
  
  alluv_pl <- ggplot(data = d5_match,
                     aes(x = day, 
                         y = freq, 
                         alluvium = Tfh_clone_id)) +
    geom_alluvium(aes(fill = Tfh_clone_id),
                  alpha = .75, 
                  width = 1/2,
                  decreasing = FALSE,
                  size = 0.7,
                  color='black') +
    theme_minimal()+
    labs(subtitle = paste0(unique(d5_match$donor), " bulk ", 
                           chain)) +
    ylab('Relative frequency for tissue and day') +
    scale_fill_manual(values = clonePal) +
    geom_text(aes(label= tissue), y= label_y_pos, size=3.5, color="black")
  
  alluv_plots[[i]] <- alluv_pl
  match_dfs[[i]] <- d5_match
  filename <- paste0('./10x/outs/Set2_321-05_pickedClone_bulk_', chain,'.pdf')
  ggsave(filename, 
         plot = alluv_pl,
         height = 4, width = 6, useDingbats = F)
}

# make a combined layout of the 10X and bulk plots =
pickedClone_alluvials <- (picked_pl_all/ alluv_plots[[2]]) + plot_layout(guides = 'collect')

ggsave('./10x/outs/Set2_321-05_pickedClone_10x_and_bulk_TCRb.pdf', 
       plot = pickedClone_alluvials,
       height = 13.5, width = 8, useDingbats = F)


### Not run ====
# long version
degs_scale_data_t <- scale_data[which(row.names(scale_data) %in% tfh_gene_set), ] %>% 
  as.data.frame() %>%
  select(rownames(sample_col)) %>% # reorder to match the df
  t()

#plot heatmap
pc_hm_long <- pheatmap(degs_scale_data_t, 
                       color = viridis::viridis(50, option = "B"),
                       cluster_rows = FALSE,
                       annotation_row  = sample_col, 
                       fontsize_row = 10, border_color= NA,
                       annotation_colors = ann_colors,
                       show_rownames  = FALSE
) %>%
  as.ggplot(.)
ggsave('./10x/outs/Set2_Donor_321-05_Tfh_pickedClone_gex_long.pdf', plot = pc_hm_long,
       height = 13, width = 6, useDingbats = F)

# or averaged gene expression
avg_exp_list<- AverageExpression(pc_obj, 
                                 group.by = c('year','day', 'Tfh_clone_id'), 
                                 slot = 'scale.data')
avg_exp_matrix_df <- avg_exp_list$RNA

degs_avg_exp_df1 <- avg_exp_matrix_df[ tfh_gene_set, ] %>% 
  as.data.frame(.) %>%
  mutate(gene = rownames(.)) %>%
  pivot_longer(1:ncol(.)-1, values_to = 'avg_exp', names_to = 'sample') %>%
  mutate( year = str_split(sample, '_', simplify = T)[,1],
          day = str_split(sample, '_', simplify = T)[,2],
          Tfh_clone_id = paste0('Tfh_', str_split(sample, 'Tfh_', simplify = T)[,2]),
          flu_specific = ifelse(Tfh_clone_id %in% pcid, 'yes','unknown'))

degs_avg_exp_mat_df <- degs_avg_exp_df1 %>%
  select(sample, gene, avg_exp) %>%
  pivot_wider(names_from = 'sample', values_from = 'avg_exp') 

degs_avg_exp_mat <- degs_avg_exp_mat_df %>%
  select(2:ncol(.)) %>%
  as.matrix() 
row.names(degs_avg_exp_mat) <- degs_avg_exp_mat_df$gene

sample_col3 <- degs_avg_exp_df1 %>%
  select(year,day, Tfh_clone_id, flu_specific, sample) %>%
  distinct_all() %>%
  as.data.frame()
rownames(sample_col3) <- sample_col3$sample
sample_col3 <- sample_col3[,1:4]
colnames(sample_col3)[3:4] <- c('Clone ID', 'Flu specific')
gaps_col <- sample_col3 %>% filter(year == 1) %>% nrow()
#plot heatmap
ann_colors3 <- ann_colors[ names(ann_colors)!='tissue' ]

pc_hm_avg <- pheatmap(degs_avg_exp_mat, 
                      color = viridis::viridis(50, option = "B"),
                      cluster_cols = FALSE,
                      annotation_col = sample_col3, 
                      fontsize_row = 10, border_color= NA,
                      annotation_colors = ann_colors3,
                      show_colnames  = FALSE, gaps_col = gaps_col) %>%
  as.ggplot(.)

ggsave('./10x/outs/Set2_Donor_321-05_Tfh_pickedClone_gex_averaged.pdf', 
       plot = pc_hm_avg,
       height = 8, width = 13, useDingbats = F)