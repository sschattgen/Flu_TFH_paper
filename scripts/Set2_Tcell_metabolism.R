source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
library(cerebroApp)
library(org.Hs.eg.db)
library(GSVA)
select <- dplyr::select
rename <- dplyr::rename

#import Tfh_object ====
Tcells <- readRDS(Set2_integrated_Tcells_path)
Tfh_obj <- readRDS(Set2_integrated_Tfh_path)
LN_TFH <- Tfh_obj@meta.data %>% 
  filter(Tfh_type != 'Treg') %>%
  pull(barcode)

tfh_call <- ifelse(colnames(Tcells) %in% LN_TFH, 'TFH', 'other')
names(tfh_call) <- colnames(Tcells)
AddMetaData(Tcells, tfh_call, 'GS_TFH')
# functions ====
pull_cerebro_enrichr_results <- function( obj, slot.name ){
  
  enriched_pw_list <- obj@misc$enriched_pathways[[slot.name]]
  
  results_all <- list()
  for (i in seq_along(enriched_pw_list)){
    
    enriched_pw_list[[i]]$group <- names(enriched_pw_list)[i]
    results_all[[i]] <-enriched_pw_list[[i]]
    names(results_all)[i] <- names(enriched_pw_list)[i]
    
  }
  
  return(results_all)
  
}

#subsetting on donor 5 ====


Tcells_d5 <- subset(Tfh_obj, subset = tissue == 'FNA')

groups <- c('Tfh_type')

# find positive markers by group
Tcells_d5 <- getMarkerGenes(Tcells_d5, 
                         organism = 'hg',
                         groups = groups, 
                         name = 'cerebro',
                         assay = 'RNA', only_pos = T)

# scan GO biological processes for enriched pathways
Tcells_d5 <- getEnrichedPathways(Tcells_d5, 
                              marker_genes_input = 'cerebro',
                              databases = c('GO_Biological_Process_2018')
)

# pull results ====
enriched_by_group <- pull_cerebro_enrichr_results(Tcells_d5, 'cerebro_enrichr')

#plot by Tcell_type
plot_dfs <- enriched_by_group$Tfh_type %>% 
  group_by(Tfh_type) %>%
  mutate(number_genes = str_split(Overlap, '/', simplify = T)[,1]) %>%
  mutate(log.Adjusted.P.value = -log10(Adjusted.P.value)) %>%
  group_split()

plot_list <- list()
for (i in seq_along(plot_dfs)){
  
  df <- plot_dfs[[i]]%>%
    slice_max(order_by =log.Adjusted.P.value,  n =40) %>% 
    ungroup() %>%
    arrange(desc(log.Adjusted.P.value)) %>%
    mutate_if(is.factor, as.character) %>%
    mutate(Term = factor(Term, levels = Term))
  
  plot_list[[i]] <- ggplot(df, 
                           aes(log.Adjusted.P.value, Term, 
                               fill= log.Adjusted.P.value)) + 
    geom_bar(stat='identity') +
    scale_fill_viridis(option = 'E', limits = c(0,15)) +
    ggtitle(unique(plot_dfs[[i]]$Tfh_type)) +
    ylab('') +
    xlab('-log10( adjusted P-value )')
  
  
}
wrap_plots(plot_list, ncol = 2)


plot_list[[4]]

# pick pathways ======

write.csv(plot_dfs[[4]], './10x/outs/Set2_Tfh_lineages_enrichr.csv')

plot_df3 <- plot_dfs[[4]]

plot_df3 <- plot_df3 %>%
  slice_max(order_by =log.Adjusted.P.value,  n = 15) %>%
  group_by(Tfh_type) %>%
  arrange(desc(log.Adjusted.P.value), .by_group = T)

plot_df3$Term  <- factor(plot_df3$Term, levels = rev(unique(plot_df3$Term))  )

pan_GO_terms_plot2 <- ggplot(plot_df3,
                             aes(log.Adjusted.P.value, Term)) +
  geom_segment( aes(x = 0, xend=log.Adjusted.P.value,
                    yend=Term, y=Term),
                color = 'grey70', linewidth =1) +
  geom_point( size=3 ) +
  facet_wrap(~Tfh_type, nrow = 1) +
  theme_clean() +
  xlab('-log10( adjusted P-value )')+
  ylab('') +
  NoLegend()

# ggsave('./10x/outs/Set2_Tcells_Tcell_type_Tfh_only_enrichr.pdf',
#        plot = pan_GO_terms_plot2,
#        width = 630/72, height = 315/72, useDingbats =F )
# 



#++++

Tfh_paths <- plot_df3#bind_rows(plot_dfs[[3]], plot_dfs[[4]])
TFH_gene_sets <- list()
for (i in unique(Tfh_paths$Term)){
  
  gs <- Tfh_paths %>%
    filter(Term == i) %>%
    pull(Genes) %>% 
    str_split(.,';') %>% 
    unlist()
  TFH_gene_sets[[i]] <- gs
  
}

# split object by donor
Donors <- c('321-04','321-05')
donor_objs <-list()
for (i in seq(Donors)){
  
  Tfh_d <- subset(Tfh_obj, subset = donor == Donors[i] & tissue == 'FNA' & Tfh_type != 'Treg')
  DefaultAssay(Tfh_d) <- 'RNA'
  Tfh_d<- SeuratReprocess(Tfh_d)
  donor_objs <-list.append( donor_objs , Tfh_d  )
  names(donor_objs)[i] <- Donors[i]
  
}

donor_objs_wTreg <-list()
for (i in seq(Donors)){
  
  Tfh_d <- subset(Tfh_obj, subset = donor == Donors[i] & tissue == 'FNA')
  DefaultAssay(Tfh_d) <- 'RNA'
  Tfh_d<- SeuratReprocess(Tfh_d)
  donor_objs_wTreg <-list.append( donor_objs_wTreg , Tfh_d  )
  names(donor_objs_wTreg)[i] <- Donors[i]
  
}

#make matrix for gsva averaged over cells in the group

groups_tfh <- c('Tfh_type','time_point','ident')
GSVA_wrapper <- function(obj, gene_sets, groups){
  # gene_sets is a names list of gene sets, each entry a vector
  # groups is a vector of grouping factors
  
  avg_matrices <- list()
  GSVA_results <- list()
  
  expressed_genes <- Matrix::rowSums(obj@assays[['RNA']]@data)
  expressed_genes <- which(expressed_genes != 0)
  
  for ( i in seq(groups)){
    
    current_group <- groups[i]
    
    #make matrix for gsva averaged over cells in the group
    
    avg_expression <- AverageExpression(obj, assays = 'RNA', 
                                        features = names(expressed_genes),slot = 'data',
                                        group.by = current_group)
    
    avg_mat <- do.call(cbind, avg_expression) %>% as.matrix() 
    avg_mat2 <- avg_mat[which(rowSums(avg_mat) > 0), ]
    colnames(avg_mat2) <- gsub('RNA.', '', colnames(avg_mat2))
    
    #store it 
    avg_matrices[[i]] <- avg_mat2
    names(avg_matrices)[i] <- current_group
    
    #run gsva and store it
    
    GSVA <- gsva( avg_mat2 , gene_sets )
    GSVA_results[[i]] <- GSVA
    names(GSVA_results)[i] <- current_group
    
    
  }
  
  return(GSVA_results)
  
}
donor_gsva <- map(donor_objs, ~GSVA_wrapper(., TFH_gene_sets, groups_tfh))
donor_Treg_gsva <- map(donor_objs_wTreg, ~GSVA_wrapper(., TFH_gene_sets, groups_tfh))
Tfh_gsva <- GSVA_wrapper(Tfh_obj, TFH_gene_sets, groups_tfh )

#make heatmap of gsva

ann_colors = list(
  day = TimePal,
  year = YearPal
)

#heatmap by donor and across time
gsva_hm_list <- list()
for ( i in seq_along(donor_gsva)){
  
  sample_col <- FetchData(donor_objs[[names(donor_gsva)[i]]], c('time_point','day','year'))
  sample_col <-sample_col%>% distinct_all()
  rownames(sample_col)<-sample_col$time_point
  sample_col <- sample_col %>% select(-time_point)
  gaps_col <- sample_col %>% filter(year ==1) %>% nrow()
  
  gsva_hm_d <- pheatmap(donor_gsva[[i]]$time_point, 
                        color = rev(brewer.rdbu(50)),
                        cluster_cols = FALSE,
                        annotation_col = sample_col, 
                        fontsize_row = 10, border_color= NA,
                        annotation_colors = ann_colors,
                        show_colnames  = FALSE, 
                        main = paste0('Donor ', names(donor_gsva)[i]),
                        gaps_col = gaps_col) %>%
    as.ggplot(.)
  gsva_hm_list[[names(donor_gsva)[i]]] <- gsva_hm_d
  
}

ggsave('./10x/outs/Set2_Tfh_lineages_donor4_time_metabolism_GSVA_HM_noTreg.pdf', 
       gsva_hm_list[[1]], width = 829/72, height =  349/72, useDingbats = F)
ggsave('./10x/outs/Set2_Tfh_lineages_donor5_time_metabolism_GSVA_HM_noTreg.pdf', 
       gsva_hm_list[[2]], width = 829/72, height =  349/72, useDingbats = F)

#heatmap by all
sample_col2 <- FetchData(Tfh_obj, c('time','day','year'))
sample_col2 <-sample_col2 %>% 
  distinct_all() %>%
  filter(!is.na(time))

rownames(sample_col2)<-sample_col2$time


sample_col2 <- sample_col2 %>% select(-time)
gaps_col <- sample_col2 %>% filter(year ==1) %>% nrow()
gsva_hm <- pheatmap(Tfh_gsva$time, 
                    color = rev(brewer.rdbu(50)),
                    cluster_cols = FALSE,
                    annotation_col = sample_col2, 
                    fontsize_row = 10, border_color= NA,
                    annotation_colors = ann_colors,
                    show_colnames  = FALSE, 
                    main = 'All donors', gaps_col = gaps_col) %>%
  as.ggplot(.)


ggsave('./10x/outs/Set2_Tfh_lineages_all_donors_time_metabolism_GSVA_HM_noTreg.pdf', 
       gsva_hm, width = 829/72, height =  349/72, useDingbats = F)


#heatmap by donor and across tfh_type
ann_colors2 = list(
  Tfh_type = TfhPal
)
gsva_hm_list2 <- list()
for ( i in seq_along(donor_gsva)){
  
  sample_col <- FetchData(donor_objs_wTreg[[names(donor_Treg_gsva)[i]]], c('Tfh_type'))
  sample_col <-sample_col%>% distinct_all()
  rownames(sample_col)<-sample_col$Tfh_type
  gsva_hm_d <- pheatmap(donor_Treg_gsva[[i]]$Tfh_type, 
                        color = rev(brewer.rdbu(50)),
                        cluster_cols = FALSE,
                        annotation_col = sample_col, 
                        fontsize_row = 10, border_color= NA,
                        annotation_colors = ann_colors2,
                        show_colnames  = FALSE,
                        main = paste0('Donor ', names(donor_gsva)[i])) %>%
    as.ggplot(.)
  gsva_hm_list2[[names(donor_Treg_gsva)[i]]] <- gsva_hm_d
  
}
ggsave('./10x/outs/Set2_Tfh_lineages_donor4_Tfh_type_metabolism_GSVA_HM.pdf', 
       gsva_hm_list2[[1]], width = 794/72, height =  261/72, useDingbats = F)
ggsave('./10x/outs/Set2_Tfh_lineages_donor5_Tfh_type_metabolism_GSVA_HM.pdf', 
       gsva_hm_list2[[2]], width = 794/72, height =  261/72, useDingbats = F)

#heatmap by all
sample_col3 <- FetchData(Tfh_obj, c('Tfh_type'))
sample_col3 <-sample_col3%>% distinct_all()
rownames(sample_col3)<-sample_col3$Tfh_type

gsva_hm_type <- pheatmap(Tfh_gsva$Tfh_type, 
                         color = rev(brewer.rdbu(50)),
                         cluster_cols = FALSE,
                         annotation_col = sample_col3, 
                         fontsize_row = 10, border_color= 'grey90',
                         annotation_colors = ann_colors2,
                         show_colnames  = FALSE,
                         main = 'All donors') %>%
  as.ggplot(.)

ggsave('./10x/outs/Set2_Tfh_lineages_all_donors_Tfh_type_metabolism_GSVA_HM.pdf', 
       gsva_hm_type, width = 794/72, height =  261/72, useDingbats = F)


