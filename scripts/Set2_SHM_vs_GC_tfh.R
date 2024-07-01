source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
rename <- dplyr::rename
library(scales)


shm_files <-c(
  './10x/objects/P04_MutFreq_AgGC_Y1Y2.xlsx',
  './10x/objects/P05_MutFreq_AgGC_Y1Y2.xlsx')


read_files <- function(f){
  
  tdf <- readxl::read_xlsx(f) 
  colnames(tdf) <- gsub('_','_d',colnames(tdf))
  colnames(tdf)[1] <- 'clone_id'
  tdf <- tdf%>%
    pivot_longer(2:ncol(.), names_to = 'time_point', values_to = 'shm') %>%
    filter(!is.na(shm) & grepl('y1',time_point)) %>%
    mutate(time_point = factor(time_point, levels = names(TimePal2))) %>%
    arrange(time_point)
  return(tdf)
}

shm_df <- map(shm_files, ~read_files(.))
all(head(shm_df[[1]]) == head(shm_df[[2]]))


freq_df <- read.csv('./10x/outs/Set2_Tfh_lineages_tfh_type_time_freqs.csv') %>%
  filter(tissue =='FNA' & grepl('y1',time_point)) %>%
  mutate(time_point = factor(time_point, levels = names(TimePal2))) %>%
  rename(shm = GC) %>%
  group_by(donor ) %>%
  group_split()


plist <- list()
for(i in seq_along(shm_df)){
  
  plist[[i]] <- ggplot() +
    geom_point(data =  freq_df[[i]],aes(time_point, shm), color = 'red', size = 2) +
    geom_boxplot(data = shm_df[[i]],aes(time_point, shm), fill = NA) +
    geom_text(
      data =  shm_df[[i]] %>% 
        group_by(time_point) %>% 
        summarise(y = max(shm), 
                  median = round(median(shm), 3), N = n()),
      aes(x = time_point, y = y,label = median), 
      position = position_dodge(width = 1), 
      vjust = -0.5, 
      size = 3, 
      stat = "unique", 
      parse = TRUE
    ) +
    scale_y_continuous(
      "nucleotide mutation frequency\nin IGHV of QIV-specific lineage", 
      sec.axis = sec_axis(~ ., name = "GC TFH frequency"),
    ) +
    theme_bw() +
    ggtitle(unique(freq_df[[i]]$donor)) +
    xlab('Time Point') +
    theme(axis.title.y.right = element_text(color = 'red'))
  
  
}


ggsave('./10x/outs/Set2_SHM_vs_Tfh_freq.png',
       plot = wrap_plots(plist, ncol = 1),
       width =6, height = 5)


ggsave('./10x/outs/Set2_SHM_vs_Tfh_freq.pdf',
       plot = wrap_plots(plist, ncol = 1),
       width =6, height = 5)
