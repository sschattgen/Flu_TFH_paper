source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

# 
all_Tcells <- readRDS(Set2_integrated_Tcells_path)


# take cells with clone id ====
clone_df <- makeClonedf(all_Tcells, other_metadata = c('donor','dataset')) %>%
  group_by(dataset) %>%
  group_split()



#here we are making the clone_ids between TwoYear and d05 parsimonious ===
# the 

TwoYear_clones <- clone_df[[2]] %>% 
  select(clone_id, match_cdr) %>% 
  distinct_all()

join_clones <- clone_df[[3]] %>%
  select(barcode, match_cdr) %>%
  left_join(.,TwoYear_clones)

matched_clones <- join_clones %>%
  filter(!is.na(clone_id)) %>%
  select(-match_cdr)

new_start <- str_split(TwoYear_clones$clone_id, '_', simplify = T)[,2] %>% 
  as.numeric() %>% 
  max() + 1

new_clones <- join_clones %>%
  filter(is.na(clone_id)) %>%
  select(-clone_id)

new_clone_ids <- new_clones %>% 
  distinct_at(.vars = 'match_cdr') 

new_clone_ids$clone_id <- paste0('clone_',seq(new_start, new_start + nrow(new_clone_ids)-1),'_321-05')


new_clones <- left_join(new_clones, new_clone_ids)%>%
  select(-match_cdr)

fixed_donor5and4 <- do.call(
  bind_rows,
  list(select(clone_df[[2]],c(1,3)), matched_clones, new_clones)
) 


# call donor using clone_id for multiplexed samples ====


donor_dfs <- clone_df[[1]] %>%
  select(clone_id, donor,barcode) %>%
  group_by(donor) %>%
  group_split()


donor_dfs[[3]]$donor <- case_when(
  donor_dfs[[3]]$clone_id %in% donor_dfs[[1]]$clone_id ~ unique(donor_dfs[[1]]$donor),
  donor_dfs[[3]]$clone_id %in% donor_dfs[[2]]$clone_id ~ unique(donor_dfs[[2]]$donor),
  donor_dfs[[3]]$clone_id %in% donor_dfs[[2]]$clone_id & donor_dfs[[3]]$clone_id %in% donor_dfs[[1]]$clone_id ~ 'twodonor',
  .default = NA
)


fixed_new_donor <- do.call(bind_rows, donor_dfs) %>%
  select(-donor)



#dont forget to add the clone_id na barcodes!!!!!! ====

na_clones <- FetchData( all_Tcells, c('barcode','clone_id')) %>%
  filter(is.na(clone_id))



# fixed clones going in ====

fixed_all_cells <- do.call(bind_rows, list(na_clones, fixed_new_donor, fixed_donor5and4))

new_cids <- pull(fixed_all_cells, clone_id)
names(new_cids) <- pull(fixed_all_cells, barcode)

all_Tcells <- AddMetaData(all_Tcells, metadata = new_cids, col.name = 'clone_id')


#update time point too =====

all_Tcells@meta.data$time_point <- factor(paste0('y',all_Tcells@meta.data$year,'_d',all_Tcells@meta.data$day), levels = names(TimePal2))


#save ====
saveRDS(all_Tcells , Set2_integrated_Tcells_path )









