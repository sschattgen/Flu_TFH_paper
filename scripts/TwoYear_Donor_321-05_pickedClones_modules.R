

# BothDonor clonotype analysis 
# the alluvial plots of clonal lineage abundance for the picked clones is in here, too

source('scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

#import ====
Tfh_obj <- readRDS(Tfh_lineages_path)
Clonedf <- read.delim( clone_df_path , stringsAsFactors = F)

#Define the markers for each subset and score as a module =====
IL10_PM <- read.delim('outs/TwoYear_Tfh_lineages_IL10_TFH_vs_premem_markers.tsv') %>%
  arrange(desc(avg_log2FC))
IL10_GC <- read.delim('outs/TwoYear_Tfh_lineages_IL10_TFH_vs_GC_markers.tsv')%>%
  arrange(desc(avg_log2FC))
EM_GC <- read.delim('outs/TwoYear_Tfh_lineages_premem_vs_GC_markers.tsv')%>%
  arrange(desc(avg_log2FC))

IL10_PM_c <- IL10_PM %>% filter(avg_log2FC > 0 ) %>% pull(gene)
IL10_GC_c <- IL10_GC %>% filter(avg_log2FC > 0 ) %>% pull(gene)
IL10_marks <- intersect(IL10_PM_c,IL10_GC_c)

GC_PM_c <- EM_GC %>% filter(avg_log2FC < 0 ) %>% pull(gene)
GC_IL10_c <- IL10_GC %>% filter(avg_log2FC < 0 ) %>% pull(gene)
GC_marks <- intersect(GC_PM_c,GC_IL10_c)

PM_GC_c <- EM_GC %>% filter(avg_log2FC > 0 ) %>% pull(gene)
PM_IL10_c <- IL10_PM %>% filter(avg_log2FC < 0 ) %>% pull(gene)
PM_marks <- intersect(PM_GC_c,PM_IL10_c)

Tfh_obj <- AddModuleScore(Tfh_obj, list(IL10_marks), name = 'IL10_TFH_module')
Tfh_obj <- AddModuleScore(Tfh_obj, list(GC_marks), name = 'GC_TFH_module')
Tfh_obj <- AddModuleScore(Tfh_obj, list(PM_marks), name = 'premem_TFH_module')
colnames(Tfh_obj@meta.data) <- gsub('module1', 'module', colnames(Tfh_obj@meta.data)  )

module_gene_sets<-data.frame(
  subset = c(rep('IL10 TFH', length(IL10_marks)),
             rep('GC TFH', length(GC_marks)),
             rep('pre/memory TFH', length(PM_marks))),
  gene = c(IL10_marks,GC_marks,PM_marks)
)

write_tsv(module_gene_sets, 'outs/TwoYear_Tfh_subset_module_markers.tsv')


#heatmap for all clones =====

ann_colors = list(
  'Clone ID' = clonePal,
  tissue = TissuePal,
  'Flu specific' = flu_pal,
  day = TimePal,
  year= YearPal,
  Tfh_type =TfhPal
)

# annotation metadata

pcid <- names(clonePal)[c(1,3,11,12)]
Tfh_obj@meta.data$flu_specific <- ifelse(Tfh_obj@meta.data$Tfh_clone_id %in% pcid, 'yes','unknown')

sample_col <- FetchData(Tfh_obj, c('Tfh_clone_id', 'flu_specific' , 'year', 'day', 'tissue', 'time', 'Tfh_type',
                                  'IL10_TFH_module', 'GC_TFH_module','premem_TFH_module', 'barcode' )) %>%
  filter(Tfh_clone_id != 'other') %>%
  arrange( time , Tfh_clone_id)%>%
  select(-time)
colnames(sample_col)[1:2] <- c('Clone ID', 'Flu specific')

# matrix of module scores
mod_mat <- sample_col %>% 
  select( premem_TFH_module, GC_TFH_module,IL10_TFH_module) %>% 
  as.matrix() %>% 
  t()

sample_md <- sample_col %>% 
  select(`Clone ID`, `Flu specific`, Tfh_type, tissue, day, year)

#gaps for year
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

#plot heatmap
pc_hm <- pheatmap(mod_mat, 
                  color = viridis::viridis(50, option = "B"),
                  cluster_cols = FALSE,
                  cluster_rows  = FALSE,
                  annotation_col = sample_md, 
                  fontsize_row = 10, border_color= NA,
                  annotation_colors = ann_colors,
                  show_colnames  = FALSE, gaps_col= gaps_col) %>%
  as.ggplot(.)
ggsave('outs/TwoYear_Donor_321-05_Tfh_pickedClone_modules.pdf', plot = pc_hm,
       height = 8, width = 13, useDingbats = F)


# heatmap for indiv clones ====
lineages <- sample_col %>%
  group_by(`Clone ID`) %>%
  group_split()

lineage_hm <- list()
for ( i in seq_along(lineages)){
  
  df <- lineages[[i]] %>% as.data.frame()
  colnames(df)[1:2] <- c('Clone ID', 'Flu specific')
  rownames(df) <- df$barcode
  md_df <- df %>% select(`Clone ID`, `Flu specific`, Tfh_type, tissue, day, year)

  gaps_col <- df %>% filter(year == 1) %>% nrow()
  
  # matrix of module scores
  mod_mat <- df %>% 
    select( premem_TFH_module, GC_TFH_module,IL10_TFH_module) %>% 
    as.matrix() %>% 
    t()
  
  
  #plot heatmap
  pc_hm <- pheatmap(mod_mat, 
                    color = viridis::viridis(50, option = "B"),
                    cluster_cols = FALSE, cluster_rows = FALSE,
                    annotation_col = md_df, 
                    fontsize_row = 10, border_color= NA,
                    annotation_colors = ann_colors,
                    show_colnames  = FALSE, 
                    annotation_legend = FALSE,
                    main = unique(df$`Clone ID`), gaps_col= gaps_col) %>%
    as.ggplot(.)
  
  lineage_hm <- list.append(lineage_hm , pc_hm)
}
lineage_hm_lay <- wrap_plots(lineage_hm) 
ggsave('outs/TwoYear_Donor_321-05_Tfh_pickedClones_lineages_module_hm.pdf', plot = lineage_hm_lay,
       height = 12.5, width = 17, useDingbats = F)

