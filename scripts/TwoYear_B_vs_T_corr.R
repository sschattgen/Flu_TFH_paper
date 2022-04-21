source('scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

#Bcells <-readRDS( Bcells_path)
Tcells <- readRDS( Tcells_path )
Tfh_only <- readRDS( Tfh_lineages_path )

# estimate frequencies by day and tissue and donor ====

# B_table <- FetchData(Bcells, vars = c('donor', 'time', 'tissue', 'B_cell_type'))
# 
# B_tally <- B_table %>% 
#   group_by_all() %>% 
#   tally() %>% 
#   ungroup(B_cell_type) %>%
#   mutate( freq = n / sum(n)) 
# colnames(B_tally)[4] <- 'subset'
# write.csv(B_tally, 'outs/TwoYear_LN_Bcell_freq_time.csv')

B_tally <- read.csv('outs/TwoYear_LN_Bcell_freq_time.csv')

B_tally <- B_tally %>% 
  filter(tissue == 'FNA') %>%
  mutate(cell_type = 'B_cell')

# focus on the "true" Tfh cells
T_table <- FetchData(Tcells, vars = c('donor', 'time', 'tissue', 'Tfh_type'))

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
  geom_point(aes(color = donor, fill = time), size = 4, shape =21,stroke =2) +
  scale_fill_manual(values = TimePal2) +
  scale_color_manual(values = DonorPal[1:2]) +
  geom_smooth(method = 'lm', se=FALSE) +
  stat_cor(label.x = 0.1, label.y = 0.1) +
  theme_minimal() +
  xlab('Tfh freq') +
  ylab('GC B cell freq') +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        title = element_text(size = 22),
        legend.text = element_text(size = 12))

ggsave('outs/TwoYear_GC_vs_Tfh_corr.pdf', plot = TfhvsGC,
       width = 7 , height = 6, useDingbats = F)

#correlation matrix
freq_mat <- wdf %>%
  ungroup()%>%
  select(-donor, -time) %>%
  as.matrix(.)
colnames(freq_mat) <- gsub("_._cell","", colnames(freq_mat))
freq_mat[ which( is.na(freq_mat)) ]  <- 0 
T_B_cor <- cor(freq_mat)
T_B_cor_hm <- pheatmap(T_B_cor) %>%
  as.ggplot() +
  labs(title = 'Pearson correlation between Tfh and B cell subset frequencies')
ggsave('outs/TwoYear_GC_vs_Tfh_corr_heatmap.pdf', 
       plot = T_B_cor_hm,
       width = 6 , height = 6, useDingbats = F)
