

# BothDonor clonotype analysis 
# the alluvial plots of clonal lineage abundance for the picked clones is in here, too

source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

#import ====
Tfh_obj <- readRDS(Set2_integrated_Tfh_path)
Clonedf <- read.delim( clone_df_path , stringsAsFactors = F)

#Define the markers for each subset and score as a module =====


IL10_marks <- read.delim('./10x/outs/Set2_Tfh_lineages_IL10_TFH_vs_premem_markers.tsv') 
GC_marks <- read.delim('./10x/outs/Set2_Tfh_lineages_GC_vs_premem_markers.tsv')
Treg_marks <- read.delim('./10x/outs/Set2_Tfh_lineages_Treg_vs_premem_markers.tsv')

marker_set <-map( list('IL10'= IL10_marks, 'GC' = GC_marks, 'Treg' =Treg_marks) ,
     ~ filter(., gene %notin% IgTcr$genes & grepl('^RP[LS][1-9]',gene) ==F ) %>% 
       arrange(p_val_adj) )



non_PM_pos_marks <- list()
PM_pos_marks <- c()

for(i in seq_along(marker_set)){
  
  
  non_PM_pos_marks[[names(marker_set)[i]]] <- marker_set[[i]] %>% 
    filter(avg_log2FC > 0 ) %>% 
    pull(gene)
  
  PM_tmp <- marker_set[[i]] %>% 
    filter(avg_log2FC < 0 ) %>% 
    pull(gene)
  
  PM_pos_marks <- union(PM_pos_marks, PM_tmp)
  
  
  
  
}




trim <- 300

umark_set <- list(
  'IL10' = setdiff(non_PM_pos_marks$IL10[1:trim] , c(non_PM_pos_marks$GC, non_PM_pos_marks$Treg, PM_pos_marks)),
     'GC' = setdiff(non_PM_pos_marks$GC , c(non_PM_pos_marks$IL10[1:trim], non_PM_pos_marks$Treg,PM_pos_marks)),
     'premem' = setdiff(PM_pos_marks, c(non_PM_pos_marks$GC, non_PM_pos_marks$Treg,non_PM_pos_marks$IL10[1:trim])),
     'Treg' = setdiff(non_PM_pos_marks$Treg , c(non_PM_pos_marks$IL10[1:trim], non_PM_pos_marks$GC,PM_pos_marks))
     )


Tfh_obj@meta.data$IL10_TFH_module <- NULL
Tfh_obj@meta.data$GC_TFH_module <- NULL
Tfh_obj@meta.data$premem_TFH_module <- NULL
Tfh_obj@meta.data$Treg_module <- NULL



Tfh_obj <- AddModuleScore(Tfh_obj, list(umark_set$IL10), nbin =20,  name = 'IL10_TFH_module')
Tfh_obj <- AddModuleScore(Tfh_obj, list(umark_set$GC), nbin =20, name = 'GC_TFH_module')
Tfh_obj <- AddModuleScore(Tfh_obj, list(umark_set$Treg), nbin =20, name = 'Treg_module')
Tfh_obj <- AddModuleScore(Tfh_obj, list(umark_set$premem), nbin =20,name = 'premem_TFH_module')

colnames(Tfh_obj@meta.data) <- gsub('module1', 'module', colnames(Tfh_obj@meta.data)  )

module_gene_sets<-data.frame(
  subset = c(rep('IL10 TFH', length(umark_set$IL10)),
             rep('GC TFH', length(umark_set$GC)),
             rep('pre/memory TFH', length(umark_set$premem)),
             rep('Treg', length(umark_set$Treg))),
  gene = c(umark_set$IL10,umark_set$GC,umark_set$premem, umark_set$Treg)
)

write_tsv(module_gene_sets, './10x/outs/Set2_Tfh_subset_module_markers.tsv')


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

sample_col <- FetchData(Tfh_obj, c('Tfh_clone_id', 'flu_specific' , 'year', 'day', 'tissue', 'time_point', 'Tfh_type',
                                  'IL10_TFH_module', 'GC_TFH_module','premem_TFH_module', 'barcode' )) %>%
  filter(Tfh_clone_id != 'other') %>%
  arrange( time_point , Tfh_clone_id)%>%
  select(-time_point)
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


ggsave('./10x/outs/Set2_Donor_321-05_Tfh_pickedClone_modules.pdf', plot = pc_hm,
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

ggsave('./10x/outs/Set2_Donor_321-05_Tfh_pickedClones_lineages_module_hm.pdf', plot = lineage_hm_lay,
       height = 12.5, width = 17, useDingbats = F)



#save ====
saveRDS(Tfh_obj, Set2_integrated_Tfh_path)


