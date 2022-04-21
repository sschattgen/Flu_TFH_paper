source('scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

library(GSVA)
library(cerebroApp)
library(org.Hs.eg.db)

select <- dplyr::select

#import Tfh_object
Tcells <- readRDS(Tcells_path)
Tfh_obj <- readRDS(Tfh_lineages_path)
# look for enriched pathways across factors ====

#groups to scan
groups <- c('Tcell_type')

# find positive markers by group
Tcells <- getMarkerGenes(Tcells, 
                          organism = 'hg',
                          groups = groups, 
                          name = 'cerebro',
                          only_pos = TRUE)

# scan GO biological processes for enriched pathways
Tcells <- getEnrichedPathways(Tcells, 
                               marker_genes_input = 'cerebro',
                               databases = c('GO_Biological_Process_2018')
                              )

# pull results ====
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
enriched_by_group <- pull_cerebro_enrichr_results(Tcells, 'cerebro_enrichr')

#plot by Tcell_type
plot_dfs <- enriched_by_group$Tcell_type %>% 
  group_by(Tcell_type) %>%
  mutate(number_genes = str_split(Overlap, '/', simplify = T)[,1]) %>%
  mutate(log.Adjusted.P.value = -log10(Adjusted.P.value)) %>%
  group_split()
plot_list <- list()
for (i in seq_along(plot_dfs)){
  
  df <- plot_dfs[[i]]%>%
    slice_max(order_by =log.Adjusted.P.value,  n = 15) %>% 
    ungroup() %>%
    arrange(desc(log.Adjusted.P.value)) %>%
    mutate_if(is.factor, as.character) %>%
    mutate(Term = factor(Term, levels = Term))
  plot_list[[i]] <- ggplot(df, 
                           aes(log.Adjusted.P.value, Term, 
                               fill= log.Adjusted.P.value)) + 
    geom_bar(stat='identity') +
    scale_fill_viridis(option = 'E', limits = c(0,15)) +
    ggtitle(unique(plot_dfs[[i]]$Tcell_type)) +
    ylab('') +
    xlab('-log10( adjusted P-value )')
    
  
}
wrap_plots(plot_list, ncol = 1)


#color palette
Tcell_type_pal <- jcolors('pal3')[1:4]
names(Tcell_type_pal) <- c('effector/memory','MAIT/NKT','naive','Tfh')

#plot the effector types
plot_df2 <- do.call(bind_rows, plot_dfs[2:4])
plot_df2 <- plot_df2 %>%
  slice_max(order_by =log.Adjusted.P.value,  n = 15) %>%
  group_by(Tcell_type) %>%
  arrange(desc(log.Adjusted.P.value), .by_group = T)
plot_df2$Term  <- factor(plot_df2$Term, levels = rev(unique(plot_df2$Term))  )

pan_GO_terms_plot <- ggplot(plot_df2, 
       aes(log.Adjusted.P.value, Term)) + 
  geom_segment( aes(x = 0, xend=log.Adjusted.P.value, 
                    yend=Term, y=Term), 
                color = 'grey70', size =1) +
  geom_point(aes(color = Tcell_type), size=3 ) +
  scale_color_manual(values = Tcell_type_pal[c(1,2,4)]) +
  facet_wrap(~Tcell_type, nrow = 1) +
  theme_clean() +
  xlab('-log10( adjusted P-value )')+
  ylab('') + 
  NoLegend()

ggsave('outs/TwoYear_Tcells_Tcell_type_enrichr.pdf', 
       plot = pan_GO_terms_plot,
       width = 1169/72, height = 664/72, useDingbats =F )

#plot the just Tfh
plot_df3 <- plot_dfs[[4]]
plot_df3 <- plot_df3 %>%
  slice_max(order_by =log.Adjusted.P.value,  n = 15) %>%
  group_by(Tcell_type) %>%
  arrange(desc(log.Adjusted.P.value), .by_group = T)
plot_df3$Term  <- factor(plot_df3$Term, levels = rev(unique(plot_df3$Term))  )

pan_GO_terms_plot2 <- ggplot(plot_df3, 
                            aes(log.Adjusted.P.value, Term)) + 
  geom_segment( aes(x = 0, xend=log.Adjusted.P.value, 
                    yend=Term, y=Term), 
                color = 'grey70', size =1) +
  geom_point(aes(color = Tcell_type), size=3 ) +
  scale_color_manual(values = Tcell_type_pal[4]) +
  facet_wrap(~Tcell_type, nrow = 1) +
  theme_clean() +
  xlab('-log10( adjusted P-value )')+
  ylab('') + 
  NoLegend()

ggsave('outs/TwoYear_Tcells_Tcell_type_Tfh_only_enrichr.pdf',
       plot = pan_GO_terms_plot2,
       width = 630/72, height = 315/72, useDingbats =F )

# Focus on some of the top pathways and look for variations in Tfh with time ====
Tfh_paths <- plot_df3
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
  
  Tfh_d <- subset(Tfh_obj, subset = donor == Donors[i] & tissue == 'FNA')
  Tfh_d<- SeuratReprocess(Tfh_d)
  donor_objs <-list.append( donor_objs , Tfh_d  )
  names(donor_objs)[i] <- Donors[i]
  
}

#make matrix for gsva averaged over cells in the group

groups_tfh <- c('Tfh_type','time','ident')
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
                                        features = names(expressed_genes), 
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
Tfh_gsva <- GSVA_wrapper(Tfh_obj, TFH_gene_sets, groups_tfh )

#make heatmap of gsva

ann_colors = list(
  day = TimePal,
  year = YearPal
)

#heatmap by donor and across time
gsva_hm_list <- list()
for ( i in seq_along(donor_gsva)){

  sample_col <- FetchData(donor_objs[[names(donor_gsva)[i]]], c('time','day','year'))
  sample_col <-sample_col%>% distinct_all()
  rownames(sample_col)<-sample_col$time
  sample_col <- sample_col %>% select(-time)
  gaps_col <- sample_col %>% filter(year ==1) %>% nrow()
  
  gsva_hm_d <- pheatmap(donor_gsva[[i]]$time, 
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
ggsave('outs/TwoYear_Tfh_lineages_donor4_time_metabolism_GSVA_HM.pdf', 
       gsva_hm_list[[1]], width = 829/72, height =  349/72, useDingbats = F)
ggsave('outs/TwoYear_Tfh_lineages_donor5_time_metabolism_GSVA_HM.pdf', 
       gsva_hm_list[[2]], width = 829/72, height =  349/72, useDingbats = F)

#heatmap by all
sample_col2 <- FetchData(Tfh_obj, c('time','day','year'))
sample_col2 <-sample_col2%>% distinct_all()
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
                    main = 'Both donors', gaps_col = gaps_col) %>%
  as.ggplot(.)


ggsave('outs/TwoYear_Tfh_lineages_both_donors_time_metabolism_GSVA_HM.pdf', 
       gsva_hm, width = 829/72, height =  349/72, useDingbats = F)


#heatmap by donor and across tfh_type
ann_colors2 = list(
  Tfh_type = TfhPal
)
gsva_hm_list2 <- list()
for ( i in seq_along(donor_gsva)){
  
  sample_col <- FetchData(donor_objs[[names(donor_gsva)[i]]], c('Tfh_type'))
  sample_col <-sample_col%>% distinct_all()
  rownames(sample_col)<-sample_col$Tfh_type
  gsva_hm_d <- pheatmap(donor_gsva[[i]]$Tfh_type, 
                        color = rev(brewer.rdbu(50)),
                        cluster_cols = FALSE,
                        annotation_col = sample_col, 
                        fontsize_row = 10, border_color= NA,
                        annotation_colors = ann_colors2,
                        show_colnames  = FALSE,
                        main = paste0('Donor ', names(donor_gsva)[i])) %>%
    as.ggplot(.)
  gsva_hm_list2[[names(donor_gsva)[i]]] <- gsva_hm_d
  
}
ggsave('outs/TwoYear_Tfh_lineages_donor4_Tfh_type_metabolism_GSVA_HM.pdf', 
       gsva_hm_list2[[1]], width = 794/72, height =  261/72, useDingbats = F)
ggsave('outs/TwoYear_Tfh_lineages_donor5_Tfh_type_metabolism_GSVA_HM.pdf', 
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
                    main = 'Both donors') %>%
  as.ggplot(.)

ggsave('outs/TwoYear_Tfh_lineages_both_donors_Tfh_type_metabolism_GSVA_HM.pdf', 
       gsva_hm_list2[[2]], width = 794/72, height =  261/72, useDingbats = F)

# score some of the interesting genesets and plot on the umap ====

#add oxphos here
#

#oxphos <- read.delim('Z:/ResearchHome/Groups/thomagrp/home/sschattg/TOOLBOX/db/human_oxphos_GO0006119.txt', header = F)
#TFH_gene_sets[['oxidative phosphorylation (GO:0006119)']] <- oxphos$V3#

#Fu_gs_all <- readxl::read_excel('./scripts/Fu_Nature_2021_Tfh_metabolism.xlsx')
#Fu_gs_Tfh <- Fu_gs_all[ which(Fu_gs_all[,4] < 0 & Fu_gs_all[,5] < 0.05  ) ,]
#Tfh_str8_match <- str_to_upper(pull(Fu_gs_Tfh[,2]))#

#TFH_gene_sets[['Fu Tfh lipid metabloism']] <- Tfh_str8_match#

#names(TFH_gene_sets) <- str_replace_all(names(TFH_gene_sets), ' ', '_')
#names(TFH_gene_sets) <- str_replace_all(names(TFH_gene_sets), '[(]', '')
#names(TFH_gene_sets) <- str_replace_all(names(TFH_gene_sets), '[)]', '')
##score gene sets
#for (ii in seq_along(TFH_gene_sets)){
#  Tfh_obj <- AddModuleScore(Tfh_obj, 
#                            features =TFH_gene_sets[ii],
#                            name = names(TFH_gene_sets)[ii] ) 
#}#
#

## overlay on umaps#

#score_df <- FetchData(Tfh_obj, c('ident', 'donor', 'time', 'day', 'year', 
#                                 'Tfh_type', 'UMAP_1','UMAP_2', 
#                     colnames(Tfh_obj@meta.data)[51:ncol(Tfh_obj@meta.data)]))#

#score_umaps <- list()
#keep_columns <- c(9,13,15,18)
#for ( i in colnames(score_df)[keep_columns]){#

#  new_title <- gsub('_', ' ', i)
#  
#  score_umap_p <- ggplot(score_df , aes(UMAP_1, UMAP_2)) +
#    geom_point(aes_string(color = i)) +
#    scale_color_distiller(palette = 'RdBu')+
#    theme_minimal() + 
#    labs(title = new_title, 
#         color= "gene module score") +
#    theme(plot.title = element_text(size = 10),
#          legend.title = element_text(size = 8))
#  
#  score_umaps[[i]] <- score_umap_p
#  
#}
#module_score_umaps <- wrap_plots(score_umaps)
#ggsave('outs/TwoYear_Tfh_lineages_GO_scores_umap.pdf', 
#       module_score_umaps, width = 973/72, height =  657/72, useDingbats = F)#
#

## tfh type umap
#Tfh_type_umap <- ggplot(score_df , aes(UMAP_1, UMAP_2)) +
#  geom_point(aes(color = Tfh_type)) +
#  scale_color_manual(values =TfhPal)+
#  theme_minimal() + 
#  labs(title = 'Tfh types') +
#  theme(plot.title = element_text(size = 10))#
#

## score over time and type====#

#score_dist_tfh_type_donor <- list()
#for (z in unique(score_df$donor)){
#  
#  local_score_df <- score_df %>% filter(donor == z)
#  score_dist_tfh_type <- list()
#  for ( i in colnames(score_df)[keep_columns]){
#    
#    new_title <- gsub('_', ' ', i)
#    
#    score_dist_p <- ggplot(local_score_df , aes_string(x='Tfh_type',y = i, fill = 'Tfh_type')) +
#      geom_boxplot() +
#      scale_fill_manual(values = TfhPal)+
#      theme_minimal() + 
#      labs(title = new_title ) +
#      ylab("gene module score") +
#      xlab('') +
#      theme(plot.title = element_text(size = 10),
#            legend.title = element_text(size = 8)) +
#      NoLegend() +
#      stat_compare_means(method = 'anova',label = 'p.signif')
#    
#    score_dist_tfh_type[[i]] <- score_dist_p
#    
#  }
#  
#  score_dist_tfh_type_donor[[z]] <- score_dist_tfh_type
#  
#}#

#score_dist_tfh_type_d4 <- wrap_plots(score_dist_tfh_type_donor[[1]]) + plot_layout(guides = 'collect')
#score_dist_tfh_type_d5 <- wrap_plots(score_dist_tfh_type_donor[[2]]) + plot_layout(guides = 'collect')#

#ggsave('outs/TwoYear_Tfh_lineages_donor4_Tfh_type_GO_scores.pdf', 
#       score_dist_tfh_type_d4, width = 995/72, height =  606/72, useDingbats = F)
#ggsave('outs/TwoYear_Tfh_lineages_donor5_Tfh_type_GO_scores.pdf', 
#       score_dist_tfh_type_d5, width = 995/72, height =  606/72, useDingbats = F)#
#

#score_dist_time <- list()#

#for ( i in colnames(score_df)[keep_columns]){
#  
#  new_title <- gsub('_', ' ', i)
#  
#  score_dist_p <- ggplot(score_df, aes_string(x='time',y = i, fill = 'time')) +
#    geom_boxplot() +
#    scale_fill_manual(values = TimePal2)+
#    theme_minimal() + 
#    labs(title = new_title ) +
#    ylab("gene module score") +
#    xlab('') +
#    theme(plot.title = element_text(size = 10),
#          legend.title = element_text(size = 8)) +
#    NoLegend()+
#    facet_wrap(~donor) + stat_compare_means(method = 'anova')
#  
#  score_dist_time[[i]] <- score_dist_p
#  
#}
#score_dist_time_wrap <- wrap_plots(score_dist_time) + 
#  plot_layout(guides = 'collect')#
#
#
#
#
#
#

## now let's look at the lipid metabolism ====
#Fu_gs_all <- readxl::read_excel('./scripts/Fu_Nature_2021_Tfh_metabolism.xlsx')
#Fu_gs_Tfh <- Fu_gs_all[ which(Fu_gs_all[,4] < 0 & Fu_gs_all[,5] < 0.05  ) ,]
#Tfh_str8_match <- str_to_upper(pull(Fu_gs_Tfh[,2]))
#donor_lipid <- map(donor_objs, ~GSVA_wrapper(., list(Fu_Tfh_lipid = Tfh_str8_match), groups_tfh))#
#
#
#

#gsva_hm_lipiid <- list()
#for ( i in seq_along(donor_lipid)){
#  
#  sample_col <- FetchData(donor_objs[[names(donor_lipid)[i]]], c('time','day','year'))
#  sample_col <-sample_col%>% distinct_all()
#  rownames(sample_col)<-sample_col$time
#  sample_col <- sample_col %>% select(-time)
#  gaps_col <- sample_col %>% filter(year ==1) %>% nrow()
#  
#  gsva_hm_lipid <- pheatmap(donor_lipid[[i]]$time, 
#                            color = rev(brewer.rdbu(50)),
#                            cluster_cols = FALSE,
#                            cluster_rows = F,
#                            annotation_col = sample_col2, 
#                            fontsize_row = 10, border_color= 'grey90',
#                            annotation_colors = ann_colors,
#                            show_colnames  = FALSE, scale = 'row',
#                            main =names(donor_lipid)[i] ) %>%
#    as.ggplot(.)
#  
#  gsva_hm_lipiid[[names(donor_lipid)[i]]] <- gsva_hm_lipid
#  
#}

# save====
saveRDS(Tfh_obj, Tfh_lineages_path)
saveRDS(Tcells, Tcells_path)



