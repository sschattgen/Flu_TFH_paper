# BothDonor clonotype analysis 
source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

# import clone_df and Tfh obj ==== 
Clonedf <- read.delim( clone_df_path , stringsAsFactors = F)

# calculate the frequency of Tfh cells by time and tissue ====

typeFreqTT <- Clonedf %>% 
  select(donor, tissue, time, Tfh_type) %>% 
  group_by_all() %>% 
  tally() %>%
  mutate(freq_tissue_time = n / sum(n)) %>%
  select(-n)

# group-split-apply to calc the number of unique clones and per cell number ====

splitClones <- Clonedf %>% 
  filter(tissue =='FNA') %>%
  select(donor, tissue, time, year, day, Tfh_type, clone_id) %>%
  group_by(donor, tissue, time, Tfh_type) %>% 
  group_split()

clonalityCalc <- function( in_list ){
  
  outdf <- data_frame()
  
  for (i in seq_along(in_list)){
    
    tmp <- in_list[[i]][1,1:ncol(in_list[[i]]) -1 ]
    
    cf <- in_list[[i]] %>% 
      add_count(clone_id) %>%
      arrange(desc(n)) %>%
      distinct_all() %>%
      mutate(freq = n /sum(n) )
    
    
    total_N <- sum(cf$n)
    big_N <- floor( total_N / 2)
    little_n = 0 
    for (z in seq(nrow(cf))){
        next_clone <- cf[z,] %>% pull(n)
      if( little_n < big_N ){
        little_n = little_n + next_clone
        nclone = z
      }
    } 
    
    D50 = nclone/nrow(cf) * 100
    
    Top10 = slice_head(cf, n =20) %>% pull(freq) %>% sum(.)
    
    tmp <- tmp %>%
      mutate( unique_clones = n_distinct( in_list[[i]]$clone_id),
              cell_number =  nrow(in_list[[i]]) ) %>%
      mutate( unique_per_cell = unique_clones / cell_number ,
              D50 = D50,
              Top10 = Top10)
      
    outdf <- bind_rows( outdf , tmp )
    
  }
  
  return(outdf)
  
}
clones_freq <- clonalityCalc( splitClones )

# bind together and subset Tfh ====

TypeCloneFreq <- full_join( clones_freq , typeFreqTT , by =c('donor', 'tissue','time', 'Tfh_type') )

Tfh_clones_freq <- TypeCloneFreq %>% 
  filter( !is.na(Tfh_type) ) %>%
  mutate(time = factor(time, levels = time_levels))

# plots ====

uclones <- ggplot( Tfh_clones_freq , aes(x= time, y= unique_per_cell, color = donor)) +
  geom_point(size = 4) +
  geom_line(aes(group= donor), size =2) +
  scale_color_manual(values = DonorPal[1:2], breaks = DonorPal[1:2]) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 30))+
  xlab("Days post vaccination") +
  ylab("Unique clonotypes /\nnumber of cells")

Tfh_clones_freq_yr <- Tfh_clones_freq %>% 
  mutate(across(where(is.factor), as.character)) %>%
  mutate(day = as.numeric(day)) %>%
  mutate(across(where(is.character), as.factor)) %>%
  group_by(year) %>%
  group_split()
  
cfclones_yr1 <- ggplot( Tfh_clones_freq_yr[[1]] , aes(x= day, y= 1/D50, color = donor)) +
  geom_point(size = 4) +
  geom_line(aes(group= donor), size =2) +
  scale_color_manual(values = DonorPal[1:2]) +
  theme_minimal() +
  scale_x_continuous( breaks=c(0,5,12,28,60,90,120,180), 
                      labels=c("0", "5", "12", "28", "60", "90", "120", "180"), 
                      limits=c(-2, 182)) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.title.x = element_blank()) +
  ylab("1 / D50 index") +
  labs(title = 'Year 1')

cfclones_yr2 <- ggplot( Tfh_clones_freq_yr[[2]] , aes(x= day, y= 1/D50, color = donor)) +
  geom_point(size = 4) +
  geom_line(aes(group= donor), size =2) +
  scale_color_manual(values = DonorPal[1:2]) +
  theme_minimal() +
  scale_x_continuous( breaks=c(0,7,14,28,60,90,120),
                      labels=c("0", "7", "14", "28", "60", "90", "120"),
                      limits=c(-2, 182)) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22)) +
  xlab("Time point") + 
  ylab("1 / D50 index") +
  labs(title = 'Year 2')

tfreq_yr1 <- ggplot( Tfh_clones_freq_yr[[1]] , aes(x= day, y= freq_tissue_time, color = donor)) +
  geom_point(size = 4) +
  geom_line(aes(group= donor), size =2) +
  scale_color_manual(values = DonorPal[1:2]) +
  theme_minimal() + 
  scale_x_continuous( breaks=c(0,5,12,28,60,90,120,180), 
                      labels=c("0", "5", "12", "28", "60", "90", "120", "180"), 
                      limits=c(-2, 182)) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.title.x = element_blank())+
  ylab("Tfh freq") +
  labs(title = 'Year 1')

tfreq_yr2 <- ggplot( Tfh_clones_freq_yr[[2]] , aes(x= day, y= freq_tissue_time, color = donor)) +
  geom_point(size = 4) +
  geom_line(aes(group= donor), size =2) +
  scale_color_manual(values = DonorPal[1:2]) +
  theme_minimal() + 
  scale_x_continuous( breaks=c(0,7,14,28,60,90,120),
                      labels=c("0", "7", "14", "28", "60", "90", "120"),
                      limits=c(-2, 182)) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22))+
  xlab("Time point") + 
  ylab("Tfh freq") +
  labs(title = 'Year 2')

freq_plots <- (tfreq_yr1 / tfreq_yr2) + 
  plot_annotation(title =  "Tfh frequency") +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')
cf_plots <- (cfclones_yr1 / cfclones_yr2) + 
  plot_annotation(title =  "1 / D50 index") +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')




freq_plots | cf_plots
# cf vs Tfh frequency ====

# color palette 
#fancy_cols <- c(TimePal, DonorPal)

#Cf_v_freq <- ggplot( Tfh_clones_freq , aes(x= freq_tissue_time, y= cum_freq_expanded)) +
#  geom_smooth(aes(group = Donor), method = 'lm', se=F, linetype = 'dashed', color = 'red') +
#  stat_cor(aes(group = Donor, color = Donor), label.x = 0.1, size =10) +
#  geom_path(aes(color = Donor), arrow = grid::arrow(), size =2) +
#  geom_point(size = 4, aes(color = factor(Day))) +
#  scale_color_manual(values = fancy_cols, breaks = fancy_cols) +
#  xlab("Tfh frequency") + ylab("Cumulative frequency\nof expanded clonotypes") +
#  theme(axis.text = element_text(size = 18))+
#  theme_minimal()


Cf_v_freq <- ggplot( Tfh_clones_freq , aes(x= freq_tissue_time, y= 1/D50)) +
  geom_point(aes(color = donor, fill = time), size = 4, shape =21, stroke =2) +
  scale_fill_manual(values = TimePal2) +
  scale_color_manual(values = DonorPal[1:2]) +
  geom_smooth(method = 'lm', se=FALSE, ) +
  stat_cor(size = 10, label.x.npc = 0.4, label.y.npc = 0) +
  theme_minimal() +
  xlab("Tfh frequency") +
  ylab('1 / D50 index') +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        title = element_text(size = 22),
        legend.text = element_text(size = 12))


#patchwork ====

tfh_plots <-  freq_plots | cf_plots | Cf_v_freq

ggsave("./10x/outs/TwoYear_Tfh_clonality_figure.pdf", 
       plot = tfh_plots,  
       width = 24 , height = 8,
       useDingbats = F)
