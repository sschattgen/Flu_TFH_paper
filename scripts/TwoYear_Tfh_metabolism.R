# Take a closer look at the Tfh metabolism genes in Fu, ..., Hongbo Chi Nature 2021
source('scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

if(!require(cerebroApp, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install('romanhaa/cerebroApp')
  require(cerebroApp, quietly = TRUE)
}
if(!require(org.Hs.eg.db, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
  require(org.Hs.eg.db, quietly = TRUE)
}

#import Tfh_object
Tfh_obj <- readRDS(Tfh_lineages_path)

# look for enriched pathways across factors ====

#groups to scan
groups <- c('ident', 'Tfh_type', 'time')

# find positive markers by group
Tfh_obj <- getMarkerGenes(Tfh_obj, 
                          organism = 'hg',
                          groups = groups, 
                          name = 'cerebro',
                          only_pos = TRUE)

# scan GO biological processes for enriched pathways
Tfh_obj <- getEnrichedPathways(Tfh_obj, 
                               marker_genes_input = 'cerebro',
                               databases = 'GO_Biological_Process_2018')

# pull results
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
enriched_by_group <- pull_cerebro_enrichr_results(Tfh_obj, 'cerebro_enrichr')

#plot by Tfh_type
plot_dfs <- enriched_by_group$Tfh_type %>% 
  group_by(Tfh_type) %>%
  slice_head(n = 15) %>%
  mutate(number_genes = str_split(Overlap, '/', simplify = T)[,1]) %>%
  mutate(log.Adjusted.P.value = -log10(Adjusted.P.value)) %>%
  group_split()
plot_list <- list()
for (i in seq_along(plot_dfs)){
  
  plot_dfs[[i]] <- plot_dfs[[i]] %>% 
    ungroup() %>%
    arrange(desc(log.Adjusted.P.value)) %>%
    mutate_if(is.factor, as.character) %>%
    mutate(Term = factor(Term, levels = Term))
  plot_list[[i]] <- ggplot(plot_dfs[[i]], 
                           aes(log.Adjusted.P.value, Term, 
                               fill= log.Adjusted.P.value)) + 
    geom_bar(stat='identity') +
    scale_fill_viridis(option = 'E', limits = c(0,15)) +
    ggtitle(unique(plot_dfs[[i]]$Tfh_type)) +
    ylab('') +
    xlab('-log10( adjusted P-value )')
  
  
}
wrap_plots(plot_list, ncol = 1)

plot_df2 <- bind_rows(plot_dfs[[2]], plot_dfs[[1]], plot_dfs[[3]])
plot_df2$Tfh_type  <- factor(plot_df2$Tfh_type, levels = names(TfhPal)  )
plot_df2 <- plot_df2 %>%
  group_by(Tfh_type) %>%
  arrange(desc(log.Adjusted.P.value), .by_group = T)
plot_df2$Term  <- factor(plot_df2$Term, levels = rev(unique(plot_df2$Term))  )

pan_GO_terms_plot <- ggplot(plot_df2, 
                            aes(log.Adjusted.P.value, Term)) + 
  geom_segment( aes(x = 0, xend=log.Adjusted.P.value, 
                    yend=Term, y=Term), 
                color = 'black') +
  geom_point(aes(color = Tfh_type), size=3 ) +
  scale_color_manual(values = TfhPal) +
  facet_wrap(~Tfh_type) +
  theme_clean() +
  xlab('-log10( adjusted P-value )')+
  ylab('') + 
  NoLegend()

ggsave('outs/TwoYear_Tfh_lineages_Tfh_type_enrichr.pdf',width = 950/72, height = 531/72 )

saveRDS(Tfh_obj, Tfh_lineages_path)

