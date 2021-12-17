source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
library(ggrepel)

# preprocess clones files ====
Clonedf <- read.delim(clone_df_path, stringsAsFactors = F)
sc_df <- read.delim('./repdata/sc/clone_tables/Ali_paired_sc_clones_JCC283.tsv')

all_clones <- read.csv('./10x/outs/AllDonors_Tfh_clones_tcrdist.csv')
all_clones$barcode <- paste0('barcode-',seq(nrow(all_clones)))

sc_clone_count <- sc_df %>%
  select(match.cdr) %>%
  group_by_all() %>%
  tally() 

x_clone_count <- Clonedf  %>%
  select(match.cdr) %>%
  group_by_all() %>%
  tally() 

tally_df <- bind_rows(sc_clone_count, x_clone_count)

tally_df2 <- tally_df %>%
  group_by( match.cdr) %>%
  summarise(count = sum(n))

all_clones <- left_join(all_clones, tally_df2) %>%
  mutate(across(c('picked', 'donor'), as.factor)) %>%
  arrange(count)

# make seurat object with TCRdist ====
#pick a barcode for each clone
empty_matrix <- matrix(, nrow = nrow(all_clones), ncol = nrow(all_clones))
rownames(empty_matrix) <- all_clones$barcode
colnames(empty_matrix) <- all_clones$barcode
#make seurat obj with empy matrix
obj <- CreateSeuratObject(counts = empty_matrix, assay = 'tcr')

for(i in colnames(all_clones)){
  new_md <- all_clones %>% pull(i)
  names(new_md) <- all_clones$barcode
  obj <- AddMetaData(obj,metadata = new_md,col.name = i)
}

#add kPC to tcr_pca
kpc <- read.csv('./10x/outs/AllDonors_Tfh_tcrdist_kpc.csv')
kpc <- kpc[,-1] %>% as.matrix()
rownames(kpc) <- all_clones$barcode
colnames(kpc) <- paste0('kPC_',seq(ncol(kpc)))
dimname = 'tcr_pca'
obj[[dimname]] <- CreateDimReducObject(embeddings = kpc, 
                                           key = dimname, 
                                           assay = DefaultAssay(obj))


obj <- FindNeighbors(obj, dims = 1:50, reduction = dimname) %>% 
  FindClusters(., dims = 1:50) %>% 
  RunUMAP(., 
          dims = 1:50, 
          reduction = dimname, 
          reduction.name = 'tcr_umap', 
          reduction.key = 'tcrUMAP_')
obj@meta.data$tcr_cluster <- obj@active.ident
#prep df for plotting ====
tcrumap_coords <- FetchData(obj, c('clone_id','tcrUMAP_1','tcrUMAP_2','count','donor', 'tcr_cluster'))

sc_df <- Clonedf %>% 
  filter(Tfh_lineage == 'yes') %>%
  filter(tissue == 'FNA') %>%
  select(barcode, time, year, day, tissue, clone_id , Tfh_clone_id_alt )

#expand back out so we can tally by day 
tcr_umap_df <- full_join(tcrumap_coords, sc_df )
tcr_umap_df <- tcr_umap_df %>% 
  filter(!is.na(tcrUMAP_1))%>% 
  filter(!is.na(time))

tcr_umap_df <- tcr_umap_df %>% 
  group_by(donor, clone_id, time) %>%
  add_tally(name ='count_day') %>%
  ungroup() %>%
  select(-barcode) %>%
  distinct_all() %>%
  arrange(year,day)

#labels
picked_clones <- tcr_umap_df$Tfh_clone_id_alt
tcr_umap_df$picked_clone <- ifelse(picked_clones =='other', NA, picked_clones)

tcr_umap_df$label <- case_when(
  picked_clones == 'Tfh_321-05_clone_1' ~ 'M2',
  picked_clones == 'Tfh_321-05_clone_3' ~ 'M1',
  picked_clones == 'Tfh_321-05_clone_12' ~ 'NP',
  picked_clones == 'Tfh_321-05_clone_11' ~ 'HA'
) 

tcr_umap_df_split <- tcr_umap_df %>%
  select(-tissue) %>%
  distinct_all() %>%
  group_by(donor) %>%
  group_split()

# for plots without
plot_df <- tcr_umap_df %>%
  select(-time, -year, -day, -count_day,-tissue) %>%
  distinct_all()


#plotting ====
tpal <- tableau20(20)
odds <- seq_len(length(tpal)) %% 2
pal <- c(tpal[which(odds == 1)], tpal[which(odds != 1)])
names(pal) <- seq(pal)-1
pal <- pal[1:length(unique(plot_df$tcr_cluster))]

#color by tcr_cluster
p1 <- ggplot(plot_df, aes(tcrUMAP_1, tcrUMAP_2, fill = tcr_cluster )) +
  geom_point(aes(size = count), shape = 21, alpha =0.8) +
  theme_minimal()+
  theme(axis.line = element_blank(),
        axis.text = element_text(size = 6)) + 
  labs(subtitle = 'TCR cluster',
       fill = "TCR cluster") +
  scale_fill_manual(values = pal)
#color by donor
p2 <- ggplot(plot_df, aes(tcrUMAP_1, tcrUMAP_2, fill = donor )) +
  geom_point(aes(size = count), shape = 21, alpha =0.8) +
  theme_minimal()+
  theme(axis.line = element_blank(),
        axis.text = element_text(size = 6)) + 
  labs(subtitle = 'Donor',
       fill = "Donor")  +
  scale_fill_manual(values = DonorPal)

layout<- p1+ p2 +plot_layout(guides = 'collect', ncol = 1)
ggsave('./10x/outs/AllDonors_Tfh_tcrdist_umap.pdf', layout, width = 5, height = 8)

# color by picked clones with contours indicating the density of tcr_clusters using all clones
plot_df_d5 <- plot_df %>% filter(donor =='321-05')

p3 <- ggplot(plot_df_d5, aes(tcrUMAP_1, tcrUMAP_2, fill = Tfh_clone_id_alt )) +
  geom_density_2d(data = plot_df,
                  aes(color =tcr_cluster, group = tcr_cluster), 
                  binwidth = 0.01, size =0.01, alpha =0.5)+
  geom_point(data =filter(plot_df_d5, Tfh_clone_id_alt =='other'), aes(size = count), shape = 21, alpha =0.8) +
  geom_point(data =filter(plot_df_d5, Tfh_clone_id_alt !='other'), aes(size = count), shape = 21, alpha =0.8) +
  geom_label_repel(data =plot_df_d5, aes(label = label),
                   size = 5, 
                   point.padding = 0,
                   min.segment.length = 1, 
                   nudge_x =  -1,
                   nudge_y =  1,
                   box.padding = 0.3,
                   alpha = 0.8) +
  theme_minimal()+
  theme(axis.line = element_blank(),
        axis.text = element_text(size = 6),
        legend.key.size = unit(0.5, 'cm')) + 
  labs(subtitle = 'Clones picked for screening',
       fill = "Clone ID")  +
  scale_fill_manual(values = clonePal) +
  scale_color_manual(values = pal) +
  guides(color = "none", fill = "none")

ggsave('./10x/outs/Donor_321_05_Tfh_tcrdist_umap_picked_clones.pdf', p3, width = 5, height = 5)
ggsave('./10x/outs/Donor_321_05_Tfh_tcrdist_umap_picked_clones.png', p3, width = 5, height = 5)


# umap plots laid out by time ====

tcr_umap_plot_fxn2 <- function(df, fill_by){
  
  tpal <- tableau20(20)
  odds <- seq_len(length(tpal)) %% 2
  pal <- c(tpal[which(odds == 1)], tpal[which(odds != 1)])
  names(pal) <- seq(pal)-1
  pal <- pal[1:length(unique(df$tcr_cluster))]
  
  if(unique(df$donor) == '321-05'){
    df <- df %>% filter(time != 'y1_d120')
  }
  if(fill_by == 'Tfh_clone_id_alt'){
    fpal = clonePal
    fpal[13] <- 'grey70'
    stitle = 'Clone ID'
    new_title = paste0('Donor ', unique(df$donor), ' , screened clones')
  } else {
    fpal = pal
    stitle = 'TCR cluster'
    new_title = paste0('Donor ', unique(df$donor))
  }
  
  time_levels <- unique(df$time)
  df$time <- factor(df$time , levels = time_levels)
  top_point_size <- max(df$count_day)

  
  #df <- df %>% arrange(desc(count_day))
  
  time_dfs <- df %>% 
    group_by(time) %>% 
    group_split()
  names(time_dfs) <- time_levels
  
  plots <- list()
  for(i in seq_along(time_dfs)){
    
    p <- ggplot(df, aes(tcrUMAP_1, tcrUMAP_2)) +
      geom_density_2d(data =df, 
                      aes(color =tcr_cluster, 
                          group = tcr_cluster), 
                      binwidth = 0.01, size =0.01, alpha =0.1)+
      geom_point(data = time_dfs[[i]], 
                 aes_string(fill = fill_by, size = 'count_day' ), 
                 shape = 21, alpha = 0.5) + 
      scale_size(limits=c(1,top_point_size), breaks = c(1,3,6,top_point_size))+
      theme_minimal() +
      scale_fill_manual(values = fpal) +
      scale_color_manual(values = pal) +
      theme(axis.title = element_blank(),
            axis.text = element_blank()) +
      labs(title = paste0('Year ',unique(time_dfs[[i]]$year),', Day ',unique(time_dfs[[i]]$day)),
           fill= stitle,
           size = 'Number of cells',
           color = 'TCR cluster') +
      guides(color = "none", fill = "none")
    
    
    if (fill_by == 'Tfh_clone_id_alt'){
      p <- p + geom_label_repel(data =time_dfs[[i]],
                      aes(label = label, fill = Tfh_clone_id_alt),
                      size = 3, 
                      point.padding = 0,
                      min.segment.length = 1, 
                      nudge_x =  -1,
                      nudge_y =  1,
                      box.padding = 0.3,
                      alpha = 0.8)
    }
    plots[[i]] <- p
    
  }
  
  return(plots)
}

time_umaps <- map(tcr_umap_df_split, ~tcr_umap_plot_fxn2(., fill = 'tcr_cluster' ))
pcid_umaps <- map(tcr_umap_df_split, ~tcr_umap_plot_fxn2(., fill = 'Tfh_clone_id_alt'))


pcid_d5 <- wrap_plots(pcid_umaps[[2]]) + plot_layout(guides = 'collect')
ggsave('./10x/outs/Donor_321_05_Tfh_tcrdist_umap_picked_clones_time_series.pdf', pcid_d5, width = 13, height = 9)

time_d5 <- wrap_plots(time_umaps[[2]]) + plot_layout(guides = 'collect')
ggsave('./10x/outs/Donor_321_05_Tfh_tcrdist_umap_time_series.pdf', time_d5, width = 13, height = 9)

time_d4 <- wrap_plots(time_umaps[[1]]) + plot_layout(guides = 'collect')
ggsave('./10x/outs/Donor_321_04_Tfh_tcrdist_umap_time_series.pdf', time_d4, width = 10, height = 9)


# cluster distribution by day ===
cluster_distro <- function(df){
  
  df <- df %>% 
    select(donor, time, year, day, tcr_cluster, count_day) %>%
    group_by(time, tcr_cluster ) %>% 
    add_tally(wt= count_day , name='cluster_count_day') %>%
    select(-count_day) %>%
    distinct_all() %>%
    ungroup(tcr_cluster) %>%
    add_tally(wt= cluster_count_day, name='total_count_day') %>%
    mutate(cluster_freq_day  = cluster_count_day / total_count_day) %>%
    arrange(year, day, tcr_cluster) 
  
  
  df$time <- gsub('y','Year ', df$time)
  df$time <- gsub('_d',', Day ', df$time)
  df$time <- factor(df$time, unique(df$time))

  p1 <- ggplot(df , aes(time, cluster_freq_day, fill =  tcr_cluster)) + 
    geom_col() + 
    scale_fill_manual(values = pal) +
    theme_minimal()+
    ylab('Frequency')+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
          axis.text.y = element_text( size = 8),
          axis.title.y = element_text( size = 12)) +
    labs(fill = 'TCR cluster',
         title = unique(df$donor))
  
  return(p1)
}
cluster_distro_plots  <- map(tcr_umap_df_split, ~cluster_distro(.))

ggsave('./10x/outs/Donor_321_04_Tfh_tcrdist_tcr_cluster_distro.pdf', 
       cluster_distro_plots[[1]], width = 12, height = 5)
ggsave('./10x/outs/Donor_321_05_Tfh_tcrdist_tcr_cluster_distro.pdf', 
       cluster_distro_plots[[2]], width = 12, height = 5)

#save tcrdist_obj ====

saveRDS(obj, './10x/objects/AllDonors_Tfh_tcrdist.rds')