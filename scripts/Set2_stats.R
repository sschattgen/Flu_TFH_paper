source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

Tcells <- readRDS(Set2_integrated_Tcells_path)
Tfh_obj <- readRDS(Set2_integrated_Tfh_path)


Tcells@meta.data %>%
  select(donor, year, day, tissue, clone_id) %>%
  group_by(donor, year, day, tissue) %>%
  add_tally(name = 'number of cells') %>%
  ungroup() %>%
  distinct_all() %>%
  select(-clone_id) %>%
  group_by(donor, year, day, tissue) %>%
  add_tally(name = 'number of clonotypes')%>%
  distinct_all() %>%
  arrange(donor,tissue, year, day) %>%
  mutate(tissue = gsub('FNA','LN',tissue)) %>%
  write_csv(., './10x/outs/Set2_Tcells_stats.csv')


Tfh_obj@meta.data %>%
  select(donor, year, day, tissue, clone_id) %>%
  group_by(donor, year, day, tissue) %>%
  add_tally(name = 'number of cells') %>%
  ungroup() %>%
  distinct_all() %>%
  select(-clone_id) %>%
  group_by(donor, year, day, tissue) %>%
  add_tally(name = 'number of clonotypes')%>%
  distinct_all() %>%
  arrange(donor,tissue, year, day) %>%
  mutate(tissue = gsub('FNA','LN',tissue))%>%
  write_csv(., './10x/outs/Set2_TFH_stats.csv')


Tfh_obj@meta.data %>% group_by(tissue) %>% tally()



tfh_clones <- makeClonedf(Tfh_obj, other_metadata = c('time_point','donor','Tfh_type')) %>% countClones()

treg_clone_id <- tfh_clones %>%
  filter(Tfh_type=='Treg') %>%
  pull(clone_id) %>%
  unique()
tfh_clone_id <- tfh_clones %>%
  filter(Tfh_type!='Treg') %>%
  pull(clone_id) %>%
  unique()



intersect(treg_clone_id, tfh_clone_id) %>% length()


tfh_clones %>%
  select(clone_id, Tfh_type, time_point, donor) %>%
  group_by_all() %>%
  tally() %>%
  ggplot(.,aes(Tfh_type,n)) + 
  geom_jitter() +
  facet_wrap(donor~time_point)