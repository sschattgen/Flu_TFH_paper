# BothDonor clonotype analysis 
source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

# import clone_df and Tfh obj ==== 
Clonedf <- read.delim( clone_df_path , stringsAsFactors = F) %>%
  mutate(Tfh_type2 = ifelse(is.na(Tfh_type) | Tfh_type =='Treg',NA, 'Tfh' )) %>%
  filter(donor %in% c('321-04','321-05')) %>%
  mutate(time_point = factor(time_point, level=names(TimePal2)))


# calculate the frequency of Tfh cells by time and tissue ====

typeFreqTT <- Clonedf %>% 
  select(donor, tissue, time_point, Tfh_type2) %>% 
  group_by_all() %>% 
  tally() %>%
  mutate(freq_tissue_time = n / sum(n)) %>%
  select(-n)

# group-split-apply to calc the number of unique clones and per cell number ====

splitClones <- Clonedf %>% 
  #filter(tissue =='FNA') %>%
  select(donor, tissue, time_point, year, day, Tfh_type2, clone_id) %>%
  group_by(donor, tissue, time_point, Tfh_type2) %>% 
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

TypeCloneFreq <- full_join( clones_freq , typeFreqTT , by =c('donor', 'tissue','time_point', 'Tfh_type2') )

time_breaks <- c(0,5,12,28,60,90,120,180,220,227,234,248,280,310,340)
tdf <- data.frame(time_point = names(TimePal2)[c(1:3,5:length(names(TimePal2)))],
                  time_encode = time_breaks)

Tfh_clones_freq <- TypeCloneFreq %>% 
  filter( !is.na(Tfh_type2) & tissue == 'FNA') %>%
  mutate(D50 = 1/D50) %>%
  left_join(., tdf)  %>%
  mutate(time_point = factor(time_point, levels = names(TimePal2))) 

Tfh_clones_freq_yr_split <- Tfh_clones_freq %>%
  group_by(year) %>%
  group_split()


# plots ====

uclones <- ggplot( Tfh_clones_freq , aes(x= time_point, y= unique_per_cell, color = donor)) +
  geom_point(size = 4) +
  geom_line(aes(group= donor), size =2) +
  scale_color_manual(values = DonorPal[1:2], breaks = DonorPal[1:2]) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 30))+
  xlab("Days post vaccination") +
  ylab("Unique clonotypes /\nnumber of cells")

uclones_TwoYr <- ggplot( Tfh_clones_freq , aes(x= time_encode, y= unique_clones, color = donor, group= donor)) +
  geom_point(size = 4) +
  geom_line(size =2) +
  scale_color_manual(values = DonorPal[1:2]) +
  theme_minimal() +
  scale_x_continuous( breaks=time_breaks, 
                      labels=c("0", "5", "12", "28", "60", "90", "120", "180", "0","7","14","28","60","90","120"), 
                      limits=c(-2, 342)) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.title.x = element_blank()) +
  ylab("Unique clones") +
  labs(title = 'TFH clone number')

cfclones_TwoYr <- ggplot( Tfh_clones_freq , aes(x= time_encode, y= D50, color = donor, group= donor)) +
  geom_point(size = 4) +
  geom_line(size =2) +
  scale_color_manual(values = DonorPal[1:2]) +
  theme_minimal() +
  scale_x_continuous( breaks=time_breaks, 
                      labels=c("0", "5", "12", "28", "60", "90", "120", "180", "0","7","14","28","60","90","120"), 
                      limits=c(-2, 342)) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.title.x = element_blank()) +
  ylab("1 / D50 index") +
  labs(title = 'TFH Clonality in LN')

tfreq_TwoYr <- ggplot( Tfh_clones_freq , aes(x= time_encode, y= freq_tissue_time , color = donor, group= donor)) +
  geom_point(size = 4) +
  geom_line(size =2) +
  scale_color_manual(values = DonorPal[1:2]) +
  theme_minimal() +
  scale_x_continuous( breaks=time_breaks, 
                      labels=c("0", "5", "12","28", "60", "90", "120", "180", "0","7","14","28","60","90","120"), 
                      limits=c(-2, 342)) +
  theme(axis.text.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 22),
        axis.title.x = element_blank()) +
  ylab("Freq of T cells") +
  labs(title = 'TFH Frequency in LN')


cfclones_yr1 <- ggplot( Tfh_clones_freq_yr_split[[1]] , aes(x= day, y= D50, color = donor)) +
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

cfclones_yr2 <- ggplot( Tfh_clones_freq_yr_split[[2]] , aes(x= day, y= D50, color = donor)) +
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

tfreq_yr1 <- ggplot( Tfh_clones_freq_yr_split[[1]] , aes(x= day, y= freq_tissue_time, color = donor)) +
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

tfreq_yr2 <- ggplot( Tfh_clones_freq_yr_split[[2]] , aes(x= day, y= freq_tissue_time, color = donor)) +
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


Cf_v_freq <- ggplot( Tfh_clones_freq , aes(x= freq_tissue_time, y= D50)) +
  geom_point(aes(color = donor, fill = time_point), size = 4, shape =21, stroke =2) +
  scale_fill_manual(values = TimePal2) +
  scale_color_manual(values = DonorPal[1:2]) +
  geom_smooth(method = 'lm', se=FALSE, ) +
  stat_cor(size = 10, label.x.npc = 0.05, label.y.npc = 1) +
  theme_minimal() +
  xlab("Tfh frequency") +
  ylab('1 / D50 index') +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        title = element_text(size = 22),
        legend.text = element_text(size = 12))


#patchwork ====

tfh_plots <-  freq_plots | cf_plots | Cf_v_freq

ggsave("./10x/outs/Set2_Tfh_clonality_figure.pdf", 
       plot = tfh_plots,  
       width = 24 , height = 8,
       useDingbats = F)
tfh_plots_new <- (tfreq_TwoYr  / cfclones_TwoYr) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
ggsave("./10x/outs/Set2_Tfh_clonality_figure_v2.pdf", 
       plot = tfh_plots_new,  
       width = 18 , height = 8,
       useDingbats = F)
