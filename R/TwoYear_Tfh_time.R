# TwoYear Tfh phenotypes over time

source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)


# import T cells Seurat object and df ----
Tfh_obj <- readRDS(Tfh_lineages_path)
Clonedf <- read.delim(clone_df_path, stringsAsFactors = F)

tfh_gene_set_c <- c('CXCR5', 'PDCD1','CCR7','SELL',
                    'VIM', 'TMSB10', 'EMP3','TCF7',
                    'TOX2','GNG4', 'CXCL13', 'IL21',
                    'FOXB1','DUSP4', 'DDIT4', 'IL10')

# split object by donor ====
Donors <- c('321-04','321-05')
donor_objs <-list()
for (i in seq(Donors)){
  
  Tfh_d <- subset(Tfh_obj, subset = donor == Donors[i] & tissue == 'FNA')
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
  fn = paste0('./10x/outs/TwoYear_Tfh_lineages_',names(donor_objs)[i],'_time_markers.tsv' )
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

ggsave('./10x/outs/TwoYear_Tfh_lineages_donor4_time_heatmap.pdf', donor_HM[[1]], width = 6, height = 10, useDingbats = F)
ggsave('./10x/outs/TwoYear_Tfh_lineages_donor5_time_heatmap.pdf', donor_HM[[2]], width = 6, height = 10, useDingbats = F)

ggsave('./10x/outs/TwoYear_Tfh_lineages_donor4_time_vln.pdf', d4_vln, width = 939/72, height = 879/72, useDingbats = F)
ggsave('./10x/outs/TwoYear_Tfh_lineages_donor5_time_vln.pdf', d5_vln, width = 939/72, height = 879/72, useDingbats = F)

# chi-squared test and plotting for phenotype enrichment ====

Tfh_ident_df <- FetchData( Tfh_obj, vars = c('donor', 'time', 'year','day', 'tissue', 'Tfh_type', 'clone_id'))

Tfh_ident_donors <- Tfh_ident_df %>%
  group_by(donor) %>%
  group_split()

tfh_freq_time <- list()
for ( i in seq_along(Tfh_ident_donors)){
  
  this_d <- unique(Tfh_ident_donors[[i]]$donor)
  df <- Tfh_ident_donors[[i]]  %>%
    select(-clone_id) %>%
    filter(tissue =='FNA') %>%
    select(time, Tfh_type) %>%
    group_by_all() %>%
    tally() %>% 
    pivot_wider(values_from = 'n', names_from = 'Tfh_type') %>%
    mutate_if(is.integer, ~nafill(.x, fill = 0)) %>%
    ungroup() %>%
    as.data.frame()
  rownames(df) <- df$time
  df <- df %>% select(-time)
  
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
    select(time, Tfh_type) %>%
    table() 
  dt <- dt[which(rowSums(dt) != 0 ), ]

  
  Tfh_freq <- Tfh_ident_donors[[i]]%>%
    select(-clone_id) %>%
    filter(tissue =='FNA') %>%
    group_by(time, Tfh_type) %>%
    add_count() %>%
    distinct_all() %>%
    ungroup(Tfh_type) %>%
    add_count(wt = n, name = 'N') %>%
    mutate(freq = n/N) %>%
    select(-n,-N) 
    
  year_dfs <- Tfh_freq %>% 
    group_by(year) %>% 
    group_split()
  
  top <- ggplot( year_dfs[[1]], aes(x= day, y= freq, color = Tfh_type )) +
    geom_point(size = 4) +
    geom_line(aes(group= Tfh_type ), size =2) +
    scale_color_manual(values = TfhPal) +
    theme_minimal() +
    xlab("Days post vaccination") +
    ylab("Frequency of Tfh cells by phenotype") +
    ggtitle('Year 1')
  bottom <- ggplot( year_dfs[[2]], aes(x= day, y= freq, color = Tfh_type )) +
    geom_point(size = 4) +
    geom_line(aes(group= Tfh_type ), size =2) +
    scale_color_manual(values = TfhPal) +
    theme_minimal() +
    xlab("Days post vaccination") +
    ylab("Frequency of Tfh cells by phenotype") +
    ggtitle('Year 2')
  
  tfh_freq_time[[i]] <- top/ bottom + 
    plot_layout(guides = 'collect') + 
    plot_annotation(title = this_d, subtitle = label)

  
  
}

ggsave('./10x/outs/TwoYear_Tfh_lineages_Tfh_type_freq_donor4.pdf', tfh_freq_time[[1]], width = 9.5, height =  9.5, useDingbats = F)
ggsave('./10x/outs/TwoYear_Tfh_lineages_Tfh_type_freq_donor5.pdf', tfh_freq_time[[2]], width =  9.5, height =  9.5, useDingbats = F)

