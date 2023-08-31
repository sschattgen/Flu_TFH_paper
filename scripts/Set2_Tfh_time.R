# TwoYear Tfh phenotypes over time

source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
library(future)
session_mem <- 1000 * 1024^2
options(future.globals.maxSize= session_mem)
rename <- dplyr::rename
# import T cells Seurat object and df ----
Tfh_obj <- readRDS(Tfh_lineages_path)
Clonedf <- read.delim(clone_df_path, stringsAsFactors = F)

Bcell_freqs <- read.csv('./10x/outs/TwoYear_LN_Bcell_freq_time.csv') %>%
  mutate(year = ifelse(grepl('y1',time_point),'1','2'),
         day = gsub('d','',str_extract(time_point, 'd[0-9]{1,3}$'))) %>%
  rename(cell_type = subset) %>%
  filter(cell_type == 'GC' & tissue =='FNA') %>%
  mutate(cell_type = 'GC B cell') %>%
  select(-n)
  
tfh_gene_set_c <- c('CXCR5', 'PDCD1','CCR7','SELL',
                    'VIM', 'TMSB10', 'EMP3','TCF7',
                    'TOX2','GNG4', 'CXCL13', 'IL21',
                    'FOXB1','DUSP4', 'DDIT4', 'IL10')

# split object by donor ====
Donors <- c('321-04','321-05')
donor_objs <-list()
for (i in seq(Donors)){
  
  Tfh_d <- subset(Tfh_obj, subset = donor == Donors[i] & tissue == 'FNA')
  DefaultAssay(Tfh_d) <- 'RNA'
  Tfh_d<- SeuratReprocess(Tfh_d)
  donor_objs <-list.append( donor_objs , Tfh_d  )
  names(donor_objs)[i] <- Donors[i]
  
}


# what are the DEGs by time  ====
donor_marks <-list()
for (i in seq_along(donor_objs)){
  
  donor_objs[[i]] <- SetIdent(donor_objs[[i]], value = 'time')
  d5_time_markers <- FindAllMarkers(donor_objs[[i]])
  donor_objs[[i]] <- SetIdent(donor_objs[[i]], value = 'ident')
  fn = paste0('./10x/outs/Set2_Tfh_lineages_',names(donor_objs)[i],'_time_markers.tsv' )
  write_tsv( d5_time_markers, fn)
  
  donor_marks[[i]] <- d5_time_markers
  names(donor_marks)[i] <- names(donor_objs)[i]
}


time_marks <- c('CXCR5', 'PDCD1','SELL','CCR7', 'EMP3','KLF2', 
                'NR4A2', 'IRF1','NFATC1', 'ID3','DUSP4', 'MAP3K8','CXCR4')


# plots by donor ====

#color palettes for heatmap
ann_colors = list(
  #tissue = TissuePal,
  day = TimePal,
  #time = TimePal2,
  donor = DonorPal,
  year = YearPal
)

donor_HM <- list()
donor_Vln <- list()
individual_Vln <- list()
GOI_str <- time_marks
for (i in seq_along(donor_objs)){
  
  Tfh_d <- donor_objs[[i]]
  
  sample_col <- FetchData(Tfh_d, c('time','year','day')) %>%
    arrange(time) %>%
    distinct_all()
  rownames(sample_col) <- sample_col$time
  
  
  time_goi <- donor_marks[[i]] %>%
    filter(p_val_adj <0.001 & 
             gene %notin% IgTcr$genes ) %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 5) %>% 
    pull(gene)
  
  
  
  # pull theaveraged, normalized expression for those genes
  avg_exp <- AverageExpression(Tfh_d, features = time_goi, group.by = 'time')
  
  #split by year
  years <- sample_col %>% 
    group_by(year) %>%
    group_split()
  
  sample_col <- sample_col %>% select( -time)
  
  #plot heatmap
  
  gaps_col <- sample_col %>% filter(year ==1) %>% nrow()
  
  day_hm <- pheatmap(avg_exp$RNA, 
                     color = rev(brewer.rdbu(50)),
                     cluster_cols = FALSE,
                     annotation_col = sample_col, 
                     fontsize_row = 10, border_color= NA,
                     annotation_colors = ann_colors,
                     show_colnames  = FALSE, scale = 'row', gaps_col = gaps_col) %>%
    as.ggplot(.)
  
  
  donor_HM1 <- wrap_plots(day_hm, nrow = 1) + plot_layout(guides = 'collect')
  donor_HM <- list.append( donor_HM, day_hm )
  
  
  plot_list_Tfh3 <- list()
  for (z in seq(GOI_str)){
    
    Tfh_Vln <- VlnPlot(Tfh_d, 
                       GOI_str[z], 
                       group.by = 'time',
                       pt.size = 0.2) + 
      scale_fill_manual(values = TimePal2) 
    
    if (z < length(GOI_str)){
      
      
      Tfh_Vln <- Tfh_Vln + 
        ylab(paste(GOI_str[z])) +
        theme(axis.title.x = element_blank(),
              axis.line = element_blank(),
              plot.title = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 8)) 
      
    }
    
    else {
      
      Tfh_Vln <- Tfh_Vln + 
        theme(axis.line = element_blank(),
              plot.title = element_blank(),
              axis.text.y = element_text(size = 8),
              axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5))+
        ylab(paste(GOI_str[z])) + 
        xlab("Time point post-vaccination") 
      
    }
    
    plot_list_Tfh3 <- list.append(plot_list_Tfh3, Tfh_Vln)
  }
  
  Vln_layout3 <- wrap_plots(plot_list_Tfh3) + plot_layout(guides = 'collect')
  donor_Vln <- list.append(donor_Vln, Vln_layout3)
  individual_Vln[[i]] <- plot_list_Tfh3
  
}

d4_vln <- wrap_plots(individual_Vln[[1]]) + 
  plot_layout(ncol = 1, guides = 'collect') & 
  theme(legend.position = "none")

d5_vln <- wrap_plots(individual_Vln[[2]]) + 
  plot_layout(ncol = 1, guides = 'collect') & 
  theme(legend.position = "none")

ggsave('./10x/outs/Set2_Tfh_lineages_donor4_time_heatmap.pdf', donor_HM[[1]], width = 6, height = 10, useDingbats = F)
ggsave('./10x/outs/Set2_Tfh_lineages_donor5_time_heatmap.pdf', donor_HM[[2]], width = 6, height = 10, useDingbats = F)

ggsave('./10x/outs/Set2_Tfh_lineages_donor4_time_vln.pdf', d4_vln, width = 939/72, height = 879/72, useDingbats = F)
ggsave('./10x/outs/Set2_Tfh_lineages_donor5_time_vln.pdf', d5_vln, width = 939/72, height = 879/72, useDingbats = F)

# chi-squared test and plotting for phenotype enrichment ====

Tfh_ident_df <- FetchData( Tfh_obj, vars = c('donor', 'time_point', 'year','day', 'tissue', 'Tfh_type', 'clone_id')) %>%
  filter(donor %in% c('321-04','321-05')) %>%
  mutate(day = factor(day, levels = names(TimePal)))

Tfh_ident_donors <- Tfh_ident_df %>%
  group_by(donor) %>%
  group_split()

tfh_freq_time <- list()
for ( i in seq_along(Tfh_ident_donors)){
  
  this_d <- unique(Tfh_ident_donors[[i]]$donor)
  df <- Tfh_ident_donors[[i]]  %>%
    select(-clone_id) %>%
    filter(tissue =='FNA') %>%
    select(time_point, Tfh_type) %>%
    group_by_all() %>%
    tally() %>% 
    pivot_wider(values_from = 'n', names_from = 'Tfh_type') %>%
    mutate_if(is.integer, ~nafill(.x, fill = 0)) %>%
    ungroup() %>%
    as.data.frame()
  rownames(df) <- df$time_point
  df <- df %>% select(-time_point)
  
  chisq <- chisq.test(df)
  
  label <-paste0("chisq(df = ", 
                 chisq$parameter,
                 " , N = ",
                 nrow(Tfh_ident_donors[[i]]), 
                 ") = ", 
                 round(chisq$statistic, digits = 2) , 
                 ", p ", 
                 format.pval(chisq$p.value)
  )
  
  dt <- Tfh_ident_donors[[i]] %>%
    select(-clone_id) %>%
    filter(tissue =='FNA') %>%
    select(time_point, Tfh_type) %>%
    table() 
  dt <- dt[which(rowSums(dt) != 0 ), ]
  
  
  Tfh_freq <- Tfh_ident_donors[[i]]%>%
    select(-clone_id) %>%
    filter(tissue =='FNA') %>%
    group_by(time_point, Tfh_type) %>%
    add_count() %>%
    distinct_all() %>%
    ungroup(Tfh_type) %>%
    add_count(wt = n, name = 'N') %>%
    mutate(freq = n/N) %>%
    select(-n,-N) %>%
    rename(cell_type = Tfh_type)
  
  
  time_breaks <- c(0,5,12,28,60,90,120,180,220,227,234,248,280,310,340)
  tdf <- data.frame(time_point = factor(names(TimePal2)[c(1:3,5:length(names(TimePal2)))], levels = names(TimePal2)[c(1:3,5:length(names(TimePal2)))]),
                    time_encode = as.integer(time_breaks))
  Tfh_freq <- bind_rows(Tfh_freq, filter(Bcell_freqs, donor == this_d)) #add B cells
  Tfh_freq <- left_join(Tfh_freq, tdf) 
  

  
  tfh_freq_time[[this_d]]$table <- Tfh_freq
  tfh_freq_time[[this_d]]$label <- label
  # top <- ggplot( year_dfs[[1]], aes(x= day, y= freq, color = Tfh_type )) +
  #   geom_point(size = 4) +
  #   geom_line(aes(group= Tfh_type ), size =2) +
  #   scale_color_manual(values = TfhPal) +
  #   theme_minimal() +
  #   xlab("Days post vaccination") +
  #   ylab("Frequency of Tfh cells by phenotype") +
  #   ggtitle('Year 1')
  # bottom <- ggplot( year_dfs[[2]], aes(x= day, y= freq, color = Tfh_type )) +
  #   geom_point(size = 4) +
  #   geom_line(aes(group= Tfh_type ), size =2) +
  #   scale_color_manual(values = TfhPal) +
  #   theme_minimal() +
  #   xlab("Days post vaccination") +
  #   ylab("Frequency of Tfh cells by phenotype") +
  #   ggtitle('Year 2')
  # 
  # tfh_freq_time[[i]] <- top/ bottom + 
  #   plot_layout(guides = 'collect') + 
  #   plot_annotation(title = this_d, subtitle = label)
  
  
  
}

TfhPal_ <- append(TfhPal, 'grey70')
names(TfhPal_)[5] <- 'GC B cell'

d4_ <- ggplot( tfh_freq_time[[1]]$table , aes(x= time_encode, y= freq, color = cell_type, group= cell_type)) +
  geom_point(size = 4) +
  geom_line(size =2) +
  scale_color_manual(values = TfhPal_) +
  theme_minimal() +
  scale_x_continuous( breaks=time_breaks[which(time_breaks %in% tfh_freq_time[[1]]$table$time_encode)], 
                      labels=c("0", "5", "12", "120", "0","14","90","120"), 
                      limits=c(-2, 342)) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.title.x = element_blank()) +
  xlab("Days post vaccination") +
  ylab("Frequency of Tfh cells by phenotype") +
  labs(title = names(tfh_freq_time)[1], subtitle = tfh_freq_time[[1]]$label)

d5_ <- ggplot( tfh_freq_time[[2]]$table , aes(x= time_encode, y= freq, color = cell_type, group= cell_type)) +
  geom_point(size = 4) +
  geom_line(size =2) +
  scale_color_manual(values = TfhPal_) +
  theme_minimal() +
  scale_x_continuous( breaks=time_breaks[which(time_breaks %in% tfh_freq_time[[2]]$table$time_encode)], 
                      labels=c("0", "5", "12","28","60","90","180", "0","7","28","60","120"), 
                      limits=c(-2, 342)) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.title.x = element_blank()) +
  xlab("Days post vaccination") +
  ylab("Frequency of Tfh cells by phenotype") +
  labs(title = names(tfh_freq_time)[2], subtitle = tfh_freq_time[[2]]$label)

ggsave('./10x/outs/Set2_Tfh_lineages_Tfh_type_freq_donor4_v2.pdf', d4_, width = 18, height = 5, useDingbats = F)
ggsave('./10x/outs/Set2_Tfh_lineages_Tfh_type_freq_donor5_v2.pdf', d5_, width = 18, height = 5, useDingbats = F)




# table of similar time points and phenotypes in FNA ====


Tfh_obj@meta.data %>%
  select(donor, Tfh_type , time_point, tissue) %>%
  group_by(donor, time_point, tissue) %>%
  add_tally(name = 'N') %>%
  group_by(donor, time_point, Tfh_type, tissue)  %>%
  add_tally(name = 'n') %>%
  mutate(freq = n/N) %>%
  distinct_all() %>%
  select(donor,time_point, Tfh_type, freq, tissue) %>%
  pivot_wider(names_from = 'Tfh_type', values_from = 'freq') %>% 
  arrange(donor, time_point) %>%
  write_csv(.,'./10x/outs/Set2_Tfh_lineages_tfh_type_time_freqs.csv')
  
  





