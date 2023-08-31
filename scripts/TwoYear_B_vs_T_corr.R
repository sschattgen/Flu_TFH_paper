source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
rename <- dplyr::rename
Bcells <-readRDS( Bcells_path)
Tcells <- readRDS( Set2_integrated_Tcells_path )
Tfh_only <- readRDS( Set2_integrated_Tfh_path )
Tfh_only@meta.data$Tfh_type2 <- ifelse(Tfh_only@meta.data$Tfh_type !='GC',NA, 'Tfh' )

# estimate frequencies by day and tissue and donor ====

B_table <- FetchData(Bcells, vars = c('donor', 'time', 'tissue', 'B_cell_type'))

B_tally <- B_table %>%
  rename(time_point = time) %>%
  group_by_all() %>% 
  tally() %>% 
  ungroup(B_cell_type) %>%
  mutate( freq = n / sum(n)) 
colnames(B_tally)[4] <- 'subset'


write_csv(B_tally, './10x/outs/TwoYear_LN_Bcell_freq_time.csv')

B_tally <- B_tally %>% 
  filter(tissue == 'FNA') %>%
  mutate(cell_type = 'B_cell')



# focus on the "true" Tfh cells

true_tfh <- Tfh_only@meta.data %>%
  filter(!is.na(Tfh_type2)) %>%
  pull(barcode)


T_table <- FetchData(Tcells, vars = c('donor', 'time_point', 'tissue')) %>%
  mutate(Tfh_type = ifelse(rownames(.) %in% true_tfh, 'Tfh', NA))

T_tally <- T_table %>% 
  group_by_all() %>% 
  tally() %>% 
  ungroup(Tfh_type) %>%
  mutate( freq = n / sum(n)) %>%
  filter(tissue == 'FNA') %>%
  mutate(cell_type = 'T_cell') %>%
  filter(!is.na(Tfh_type))
colnames(T_tally)[4] <- 'subset'


FNA_tally <- bind_rows(T_tally, B_tally)
FNA_tally$subset <- paste(FNA_tally$subset, FNA_tally$cell_type ,sep = '_')
FNA_tally$subset <- gsub('/','_', FNA_tally$subset)

wdf <- pivot_wider(FNA_tally, id_cols = 1:2, 
                   names_from = 'subset', 
                   values_from = freq) 

# generate plots

TfhvsGC <- ggplot(data = wdf, aes(x = Tfh_T_cell, y = GC_B_cell)) +
  geom_point(aes(color = donor, fill = time_point), size = 4, shape =21,stroke =2) +
  scale_fill_manual(values = TimePal2) +
  scale_color_manual(values = DonorPal[1:2]) +
  geom_smooth(method = 'lm', se=FALSE) +
  stat_cor(label.x = 0.06, label.y = 0.05) +
  theme_minimal() +
  xlab('GC Tfh freq') +
  ylab('GC B cell freq') +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        title = element_text(size = 22),
        legend.text = element_text(size = 12))

ggsave('./10x/outs/Set2_GC_vs_Tfh_corr.pdf', plot = TfhvsGC,
       width = 7 , height = 6, useDingbats = F)

#correlation matrix
freq_mat <- wdf %>%
  ungroup()%>%
  select(-donor, -time_point) %>%
  as.matrix(.)
colnames(freq_mat) <- gsub("_._cell","", colnames(freq_mat))
freq_mat[ which( is.na(freq_mat)) ]  <- 0 
T_B_cor <- cor(freq_mat)
T_B_cor_hm <- pheatmap(T_B_cor) %>%
  as.ggplot() +
  labs(title = 'Pearson correlation between GC Tfh and B cell subset frequencies')

ggsave('./10x/outs/Set2_GC_vs_Tfh_corr_heatmap.pdf', 
       plot = T_B_cor_hm,
       width = 6 , height = 6, useDingbats = F)




write_csv(wdf, './10x/outs/TwoYear_LN_Bcell_vs_Tfh_freq_time.csv')
