# BothDonor matching to TCRs of known specificity.
# this code is definitely overly complicated and I apologize if someone needs to read it.

library(ggthemes)
library(ggpubr)
library(rlist)
source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)


#import 10x TCR data ====
Clonedf <- read.delim( clone_df_path, stringsAsFactors = F)

# Data from Herati et al https://immunology.sciencemag.org/content/2/8/eaag2152.abstract
# these are all beta chains from adaptive
# this contains multiple subjects but only one year ====
Y1516 <- read.csv('./Wherry_dataset/Samples_1516Data_AsOf09162016.csv', stringsAsFactors = F)

Y1516_fh <- Y1516 %>% filter(SequenceStatus == 'In' & Cell == 'ICOS+CD38+')# remove out of frame or stops and focus on Tfh

dlist <- Y1516_fh %>%  
  select(1,4,6,13) %>%
  group_by(Subject, Day, CDR3.AminoAcid) %>%
  distinct_at(c('Subject', 'Day', 'CDR3.AminoAcid'), .keep_all = T) %>% #let's count by amino acid
  ungroup(CDR3.AminoAcid) %>%
  mutate(rel_freq = Count / sum(Count) ) %>% # recalc freq relative to donor and day
  ungroup(Day) %>%
  mutate(Day = paste0("day_",Day)) %>%
  group_split()

# criteria here is the clone was present at day7 and fold-change from 0 to 7 >= 1
Y1516_fh_spread <- list()
for (i in seq_along(dlist)){
  
  df <- dlist[[i]]
  nm <- unique(df$Subject)
  
  df2 <- df %>% 
    select(-Subject, -Count) %>% 
    pivot_wider( names_from = 'Day', values_from = 'rel_freq') %>% 
    filter(!is.na(day_7)) %>% 
    mutate( day_0 = ifelse( is.na(day_0), 0 , day_0)) %>%
    mutate(fold_change = day_7 / day_0) %>%
    filter(fold_change >= 1) %>%
    arrange(desc(day_7))

  df2$Subject <- nm 
  df2 <- df2 %>% select(Subject, 1,2,3,4 )
  Y1516_fh_spread[[i]] <- df2
  names(Y1516_fh_spread)[i] <- nm
  
}

#plots
fc_plots <- list()
for (z in seq_along(Y1516_fh_spread)){
  fc_plots[[z]] <- ggpaired(Y1516_fh_spread[[z]], "day_0", "day_7", title = paste0( names(Y1516_fh_spread)[z] ) )
}
wrap_plots(fc_plots)  
  
# bind together
Y1516_fh_clones <- do.call( rbind, Y1516_fh_spread)

# this contains six subjects over three years ====
multiYear <- read.csv('./Wherry_dataset/Samples_MultipleYearsData_AsOf09162016.csv', stringsAsFactors = F)

multiYear_fh <- multiYear %>% filter(SequenceStatus == 'In' & Cell == 'ICOS+CD38+')# remove out of frame or stops and focus on Tfh

mylist <- multiYear_fh %>%  
  select(1,2,4,6,13) %>%
  group_by(Subject, Year, Day, CDR3.AminoAcid) %>%
  distinct_at(c('Subject','Year','Day', 'CDR3.AminoAcid'), .keep_all = T) %>% #let's count by amino acid
  ungroup(CDR3.AminoAcid) %>%
  mutate(rel_freq = Count / sum(Count) ) %>% # recalc freq relative to donor and day
  ungroup(Day) %>%
  mutate(Day = paste0("day_",Day)) %>%
  group_split()

mylist2 = subset(mylist, subset = length(unique(Day)) == 1 )
mylist2 = mylist[-c(1,3,12)]

# criteria here is the clone was present at day7 and fold-change from 0 to 7 >= 1
multiYear_fh_spread <- list()
for (i in seq_along(mylist2)){
  
  df <- mylist2[[i]]
  
  nm <- paste0( unique(df$Subject), '_', unique(df$Year) )
  
  df2 <- df %>% 
    select(-Subject,-Year, -Count) %>% 
    pivot_wider( names_from = 'Day', values_from = 'rel_freq') %>% 
    filter(!is.na(day_7)) %>% 
    mutate( day_0 = ifelse( is.na(day_0), 0 , day_0)) %>%
    mutate(fold_change = day_7 / day_0) %>%
    filter(fold_change >= 1) %>%
    arrange(desc(day_7))
  
  df2$Subject <- unique(df$Subject)
  df2$Year <- unique(df$Year)
  df2 <- df2 %>% select(Subject, Year, 1,2,3,4 )
  multiYear_fh_spread[[i]] <- df2
  names(multiYear_fh_spread)[i] <- nm
    
}

#plots
fc_plots2<- list()
for (z in seq_along(multiYear_fh_spread)){
  fc_plots2[[z]] <- ggpaired(multiYear_fh_spread[[z]], "day_0", "day_7", title = paste0( names(multiYear_fh_spread)[z] ) )
}
wrap_plots(fc_plots2)  

# bind together
multiYear_fh_clones <- do.call( rbind, multiYear_fh_spread)
write_csv( multiYear_fh_clones, './Wherry_dataset/multiyear_expanding_clones.csv')



#this is the HA specific cells ====
HAtcr <- read.csv('./Wherry_dataset/Herati_HA_specific.csv', stringsAsFactors = F)
HAtcr_chains <- HAtcr %>% 
  filter(SequenceStatus == 'In') %>% # remove out of frame or stops and focus on Tfh
  pull(CDR3.AminoAcid)



#conveniently bundle all together

HAtcr2 <- HAtcr %>% filter(SequenceStatus == 'In')
all_Tfh <- do.call( rbind, list(HAtcr2 ,Y1516_fh, multiYear_fh))
write_csv( all_Tfh, './Wherry_dataset/all_Tfh_clones.csv')

# import vdjdb-2020-05-20 ====
vdjdb <- read.delim('Z:/ResearchHome/Groups/thomagrp/home/sschattg/TOOLBOX/db/vdjdb-2020-05-20/vdjdb.txt', stringsAsFactors  = F)
vdjdb_II <- vdjdb %>% 
  filter( mhc.class == 'MHCII' & species == 'HomoSapiens' & antigen.species == "InfluenzaA")
vdjdb_II_TRA <- vdjdb_II %>% filter(gene == "TRA")
vdjdb_II_TRB <- vdjdb_II %>% filter(gene == "TRB")
write_csv( vdjdb_II, './repdata/vdjdb_classII_Flu.csv')

# matching chains ====

# subset Clonedf on matches 
Clonedf_match_Y1516 <- Clonedf %>% 
  filter(  cdr3b %in% Y1516_fh_clones$CDR3.AminoAcid )
match_Y1516_ci <- Clonedf_match_Y1516 %>% 
  filter(ident %like% 'Tfh')%>% 
  pull(clone_id) %>% 
  unique()

Clonedf_match_my <- Clonedf %>% 
  filter(cdr3b %in% multiYear_fh_clones$CDR3.AminoAcid )
match_multiYear_ci <- Clonedf_match_my %>% 
  filter(ident %like% 'Tfh') %>% 
  pull(clone_id) %>% 
  unique()

#no matches to HA-specific clones
Clonedf_match_HA <- Clonedf %>% 
  filter(cdr3b %in% HAtcr_chains )
match_HA_ci <- Clonedf_match_HA %>% 
  filter(ident %like% 'Tfh' )%>% 
  pull(clone_id) %>% 
  unique()

#vdj-db class II alpha chains 
Clonedf_match_vdj_TRA <- Clonedf %>%
  filter( cdr3a %in% vdjdb_II_TRA$cdr3 )
match_vdj_TRA_ci <- Clonedf_match_vdj_TRA %>% 
  filter(ident %like% 'Tfh')%>% 
  pull(clone_id) %>% 
  unique()

#vdj-db class II beta chains 
Clonedf_match_vdj_TRB <- Clonedf %>%
  filter(cdr3b %in% vdjdb_II_TRB$cdr3 ) 
match_vdj_TRB_ci <- Clonedf_match_vdj_TRB %>%
  filter(ident %like% 'Tfh')%>% 
  pull(clone_id) %>% 
  unique()

# focus on the Tfh clones. Calc rel freqs by day then calculate the fold-change from day 0 to 12  ====
Tfh_list <- Clonedf %>% 
  #filter(Tissue == 'FNA' & ident %like% 'Tfh') %>% 
  filter(tissue == 'FNA' & !is.na(Tfh_type) ) %>% 
  select(clone_id, donor, time) %>%
  group_by_all() %>%
  add_count(clone_id) %>%
  arrange(desc(n)) %>%
  group_by_all() %>%
  distinct_all( .keep_all = T) %>% #let's count by amino acid
  ungroup(clone_id, n) %>%
  mutate(rel_freq = n / sum(n) ) %>% # recalc freq relative to donor and day
  mutate(cum_freq = cumsum(rel_freq)) %>%
  ungroup(time) %>%
  group_split()

Tfh_list_rel_abundance <- list()
for (i in seq_along(Tfh_list)){
  
    df <- Tfh_list[[i]]
    
    df2 <- df %>% 
      select(-donor,-cum_freq,-n) %>% 
      pivot_wider( names_from = 'time', values_from = 'rel_freq') %>% 
      #filter(!is.na(day_12)) %>% 
      mutate(across(2:ncol(.), ~ nafill(.x, fill = 0))) %>%
      mutate(fold_change_0_to_12 = y1_d12 / y1_d0 ) %>%
      #filter(fold_change_0_to_12 >= 1) %>%
      arrange(desc(y1_d12))
    
    df2$match_Y1516 <- ifelse(df2$clone_id %in% match_Y1516_ci, 'yes', 'no' )
    df2$match_multiYear <- ifelse(df2$clone_id %in% match_multiYear_ci, 'yes', 'no' )
    df2$vdjdb_TRA_hit <- ifelse(df2$clone_id %in% match_vdj_TRA_ci, 'yes', 'no' )
    df2$vdjdb_TRB_hit <- ifelse(df2$clone_id %in% match_vdj_TRB_ci, 'yes', 'no' )
    
    df2 <-  df2 %>% 
      mutate(match_both = ifelse(match_Y1516 == 'yes' &  match_multiYear== 'no', 'yes', 'no' )) %>%
      mutate(match_either = ifelse(match_Y1516 == 'yes' |  match_multiYear== 'yes', 'yes', 'no' ))
      
    df2$donor <- str_split(df2$clone_id, '_', simplify = T)[,3] %>% unique()
    df2 <- df2 %>% 
      select(ncol(.),1:ncol(.)-1)
    Tfh_list_rel_abundance[[i]] <- df2
    names(Tfh_list_rel_abundance)[i] <- unique(df2$donor)
    
}

d4 <- Tfh_list_rel_abundance[[1]]
d5 <- Tfh_list_rel_abundance[[2]]

#plotting the Tfh relative abundances over time

#all days donor 5

time_cols <- which( str_detect(colnames(d5), '^y._') )
non_time_cols <- which( colnames(d5) %like% 'clone_id' |
                    colnames(d5) %like% 'fold_change' |
                      colnames(d5) %like% 'match' |
                      colnames(d5) %like% 'hit' )
d5_long <- d5 %>% 
  select(non_time_cols, time_cols ) %>% 
  gather(c(length(non_time_cols)+1):ncol(.), key = 'time', value = 'Frequency') %>%
  mutate_if(is.character, as.factor) %>%
  mutate(time = factor(time, levels = names(TimePal2)))

in_times <- names(TimePal2)[which(names(TimePal2) %in% d5_long$time)]
time_num_df <- data.frame(time = in_times,
                          time_num = seq(in_times))

d5_long <- left_join(d5_long , time_num_df) %>%
  mutate(time_jit = jitter(time_num, amount = 0.2)) %>%
  arrange(time_num)

d5_long_th <- d5_long  %>% 
  filter(fold_change_0_to_12 >= 1) 

Tfh_clone_id_df <- Clonedf %>% 
  select(clone_id, Tfh_clone_id) %>%
  distinct_all()

d5_long_picked_clones <- left_join(d5_long, Tfh_clone_id_df)
d5_long_picked_clones$Tfh_clone_id <- ifelse(is.na(d5_long_picked_clones$Tfh_clone_id ), 'other',
                                             d5_long_picked_clones$Tfh_clone_id )


# color  donor 5 by matches to the Herati paper, all rel freq
d5_Tfh_clones_overlay <- ggplot( d5_long, aes(y= Frequency)) +
  geom_point(aes(x=time_jit)) +
  geom_point(data = filter( d5_long, match_either == "no"), aes(x=time_jit), color = 'grey80') +
  geom_line(data = filter( d5_long, match_either == "no"),  aes(x=time_jit, group= clone_id), linetype='dashed', color = 'grey80') +
  geom_point(data = filter( d5_long, match_either == "yes"), aes(x=time_jit), color = DonorPal[2] ) + 
  geom_line(data = filter( d5_long, match_either == "yes"),  aes(x=time_jit, group= clone_id), linetype='solid', color = DonorPal[2] ) +
  labs(title = 'All Tfh clones in 321-05', caption = "Blue = Tfh clones from Herati et al. Sci Imm 2017")+
  xlab("Time point") + ylab("Relative frequency\namong Tfh clones") +
  scale_x_continuous(breaks= time_num_df$time_num , 
                     labels= time_num_df$time, 
                     limits=c(min(time_num_df$time_num ) - 1, max(time_num_df$time_num ) + 1 )) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45 , size = 8, hjust =  1 , vjust = 1),
        axis.title = element_text( size = 8))

# color  donor 5 by matches to vdjdb
d5_Tfh_clones_overlay_vdj <- ggplot( d5_long, aes(y= Frequency)) +
  geom_point(aes(x=time_jit)) +
  geom_point(data = filter( d5_long, vdjdb_TRA_hit == "no" & vdjdb_TRB_hit == "no"), aes(x=time_jit), color = 'grey80') +
  geom_line(data = filter( d5_long, vdjdb_TRA_hit == "no" & vdjdb_TRB_hit == "no"),  aes(x=time_jit, group= clone_id), linetype='dashed', color = 'grey80') +
  geom_point(data = filter( d5_long, vdjdb_TRA_hit == "yes"), aes(x=time_jit), color = "red" ) + 
  geom_line(data = filter( d5_long, vdjdb_TRA_hit == "yes"),  aes(x=time_jit, group= clone_id), linetype='solid', color = 'red') +
  geom_point(data = filter( d5_long, vdjdb_TRB_hit == "yes"), aes(x=time_jit), color = "blue" ) + 
  geom_line(data = filter( d5_long, vdjdb_TRB_hit == "yes"),  aes(x=time_jit, group= clone_id), linetype='solid', color = 'blue') +
  labs(title = 'All Tfh clones in 321-05', caption = "Red = Tfh clones with cdr3a match in vdj-db\nBlue = Tfh clones with cdr3b match in vdj-db")+
  xlab("Time point") + ylab("Relative frequency\namong Tfh clones") +
  scale_x_continuous(breaks= time_num_df$time_num , 
                     labels= time_num_df$time, 
                     limits=c(min(time_num_df$time_num ) - 1, max(time_num_df$time_num ) + 1 )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45 , size = 8, hjust =  1 , vjust = 1),
        axis.title = element_text( size = 8))

# color donor 5 by Herati matches, top hits
d5_Tfh_clones_overlay_th <- ggplot( d5_long_th, aes(y= Frequency)) +
  geom_point(aes(x=time_jit)) +
  geom_point(data = filter( d5_long_th, match_either == "no"), aes(x=time_jit), color = 'grey80') +
  geom_line(data = filter( d5_long_th, match_either == "no"),  aes(x=time_jit, group= clone_id), linetype='dashed', color = 'grey80') +
  geom_point(data = filter( d5_long_th, match_either == "yes"), aes(x=time_jit), color = DonorPal[2]  ) + 
  geom_line(data = filter( d5_long_th, match_either == "yes"),  aes(x=time_jit, group= clone_id), linetype='solid', color = DonorPal[2]) +
  labs(title = 'Expanding Tfh clones in 321-05', caption = "Blue = Tfh clones from Herati et al. Sci Imm 2017")+
  xlab("Time point") + ylab("Relative frequency\namong Tfh clones") +
  scale_x_continuous(breaks= time_num_df$time_num , 
                     labels= time_num_df$time, 
                     limits=c(min(time_num_df$time_num ) - 1, max(time_num_df$time_num ) + 1 )) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45 , size = 8, hjust =  1 , vjust = 1),
        axis.title = element_text( size = 8))

# color donor 5 by matches to vdjdb, top hits
d5_Tfh_clones_overlay_vdj_th <- ggplot( d5_long_th, aes(y= Frequency)) +
  geom_point(aes(x=time_jit)) +
  geom_point(data = filter( d5_long_th, vdjdb_TRA_hit == "no" & vdjdb_TRB_hit == "no"), aes(x=time_jit), color = 'grey80') +
  geom_line(data = filter( d5_long_th, vdjdb_TRA_hit == "no" & vdjdb_TRB_hit == "no"),  aes(x=time_jit, group= clone_id), linetype='dashed', color = 'grey80') +
  geom_point(data = filter( d5_long_th, vdjdb_TRA_hit == "yes"), aes(x=time_jit), color = "red" ) + 
  geom_line(data = filter( d5_long_th, vdjdb_TRA_hit == "yes"),  aes(x=time_jit, group= clone_id), linetype='solid', color = 'red') +
  geom_point(data = filter( d5_long_th, vdjdb_TRB_hit == "yes"), aes(x=time_jit), color = "blue" ) + 
  geom_line(data = filter( d5_long_th, vdjdb_TRB_hit == "yes"),  aes(x=time_jit, group= clone_id), linetype='solid', color = 'blue') +
  labs(title = 'Expanding Tfh clones in 321-05', caption = "Red = Tfh clones with cdr3a match in vdj-db\nBlue = Tfh clones with cdr3b match in vdj-db")+
  xlab("Time point") + ylab("Relative frequency\namong Tfh clones") +
  scale_x_continuous(breaks= time_num_df$time_num , 
                     labels= time_num_df$time, 
                     limits=c(min(time_num_df$time_num ) - 1, max(time_num_df$time_num ) + 1 )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45 , size = 8, hjust =  1 , vjust = 1),
        axis.title = element_text( size = 8))


# What about the picked clones?

d5_Tfh_clones_picked_clones <- ggplot( d5_long_picked_clones, aes(y= Frequency, color = Tfh_clone_id)) +
  geom_point(aes(x=time_jit)) +
  geom_point(data = filter( d5_long_picked_clones, Tfh_clone_id == 'other'), aes(x=time_jit)) +
  geom_line(data = filter( d5_long_picked_clones, Tfh_clone_id == 'other'),  aes(x=time_jit, group= clone_id), linetype='dashed') +
  geom_point(data = filter( d5_long_picked_clones, Tfh_clone_id != 'other'), aes(x=time_jit), size = 1.2) +
  geom_line(data = filter( d5_long_picked_clones, Tfh_clone_id != 'other'),  aes(x=time_jit, group= clone_id), size = 1.2) +
  labs(title = 'All Tfh clones in 321-05', caption = "Red = Tfh clones with cdr3a match in vdj-db\nBlue = Tfh clones with cdr3b match in vdj-db")+
  xlab("Time point") + ylab("Relative frequency\namong Tfh clones") +
  scale_x_continuous(breaks= time_num_df$time_num , 
                     labels= time_num_df$time, 
                     limits=c(min(time_num_df$time_num ) - 1, max(time_num_df$time_num ) + 1 )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45 , size = 8, hjust =  1 , vjust = 1),
        axis.title = element_text( size = 8)) + 
  scale_color_manual(values = clonePal )


#all days donor 4

time_cols <- which( str_detect(colnames(d4), '^y._') )
non_time_cols <- which( colnames(d4) %like% 'clone_id' |
                          colnames(d4) %like% 'fold_change' |
                          colnames(d4) %like% 'match' |
                          colnames(d4) %like% 'hit' )
d4_long <- d4 %>% 
  select(non_time_cols, time_cols ) %>% 
  gather(c(length(non_time_cols)+1):ncol(.), key = 'time', value = 'Frequency') %>%
  mutate_if(is.character, as.factor) %>%
  mutate(time = factor(time, levels = names(TimePal2)))

in_times <- names(TimePal2)[which(names(TimePal2) %in% d4_long$time)]
time_num_df <- data.frame(time = in_times,
                          time_num = seq(in_times))

d4_long <- left_join(d4_long , time_num_df) %>%
  mutate(time_jit = jitter(time_num, amount = 0.5)) %>%
  arrange(time_num)

d4_long_th <- d4_long  %>% 
  filter(fold_change_0_to_12 >= 1) 

# color donor 4 by matches to the Herati paper
d4_Tfh_clones_overlay <- ggplot( d4_long, aes(y= Frequency)) +
  geom_point(aes(x=time_jit)) +
  geom_point(data = filter( d4_long, match_either == "no"), aes(x=time_jit), color = 'grey80') +
  geom_line(data = filter( d4_long, match_either == "no"),  aes(x=time_jit, group= clone_id), linetype='dashed', color = 'grey80') +
  geom_point(data = filter( d4_long, match_either == "yes"), aes(x=time_jit), color = DonorPal[1] ) + 
  geom_line(data = filter( d4_long, match_either == "yes"),  aes(x=time_jit, group= clone_id), linetype='solid', color = DonorPal[1]) +
  labs(title = 'All Tfh clones in 321-04', caption = "Green = Tfh clones from Herati et al. Sci Imm 2017")+
  xlab("Time point") + ylab("Relative frequency\namong Tfh clones") +
  scale_x_continuous(breaks= time_num_df$time_num , 
                     labels= time_num_df$time, 
                     limits=c(min(time_num_df$time_num ) - 1, max(time_num_df$time_num ) + 1 )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45 , size = 8, hjust =  1 , vjust = 1),
        axis.title = element_text( size = 8))

d4_Tfh_clones_overlay_vdj <- ggplot( d4_long, aes(y= Frequency)) +
  geom_point(aes(x=time_jit)) +
  geom_point(data = filter( d4_long, vdjdb_TRA_hit == "no" & vdjdb_TRB_hit == "no"), aes(x=time_jit), color = 'grey80') +
  geom_line(data = filter( d4_long, vdjdb_TRA_hit == "no" & vdjdb_TRB_hit == "no"),  aes(x=time_jit, group= clone_id), linetype='dashed', color = 'grey80') +
  geom_point(data = filter( d4_long, vdjdb_TRA_hit == "yes"), aes(x=time_jit), color = "red" ) + 
  geom_line(data = filter( d4_long, vdjdb_TRA_hit == "yes"),  aes(x=time_jit, group= clone_id), linetype='solid', color = 'red') +
  geom_point(data = filter( d4_long, vdjdb_TRB_hit == "yes"), aes(x=time_jit), color = "blue" ) + 
  geom_line(data = filter( d4_long, vdjdb_TRB_hit == "yes"),  aes(x=time_jit, group= clone_id), linetype='solid', color = 'blue') +
  labs(title = 'All Tfh clones in 321-04', caption = "Red = Tfh clones with cdr3a match in vdj-db\nBlue = Tfh clones with cdr3b match in vdj-db")+
  xlab("Time point") + ylab("Relative frequency\namong Tfh clones") +
  scale_x_continuous(breaks= time_num_df$time_num , 
                     labels= time_num_df$time, 
                     limits=c(min(time_num_df$time_num ) - 1, max(time_num_df$time_num ) + 1 )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45 , size = 8, hjust =  1 , vjust = 1),
        axis.title = element_text( size = 8))

# color donor 4 by Herati matches, top hits
d4_Tfh_clones_overlay_th <- ggplot( d4_long_th, aes(y= Frequency)) +
  geom_point(aes(x=time_jit)) +
  geom_point(data = filter( d4_long_th, match_either == "no"), aes(x=time_jit), color = 'grey80') +
  geom_line(data = filter( d4_long_th, match_either == "no"),  aes(x=time_jit, group= clone_id), linetype='dashed', color = 'grey80') +
  geom_point(data = filter( d4_long_th, match_either == "yes"), aes(x=time_jit), color = DonorPal[1]) + 
  geom_line(data = filter( d4_long_th, match_either == "yes"),  aes(x=time_jit, group= clone_id), linetype='solid', color = DonorPal[1]) +
  labs(title = 'Expanding Tfh clones in 321-04', caption = "Red = Tfh clones with cdr3a match in vdj-db\nBlue = Tfh clones with cdr3b match in vdj-db")+
  xlab("Time point") + ylab("Relative frequency\namong Tfh clones") +
  scale_x_continuous(breaks= time_num_df$time_num , 
                     labels= time_num_df$time, 
                     limits=c(min(time_num_df$time_num ) - 1, max(time_num_df$time_num ) + 1 )) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45 , size = 8, hjust =  1 , vjust = 1),
        axis.title = element_text( size = 8))

# color donor 4 by matches to vdjdb, top hits
d4_Tfh_clones_overlay_vdj_th <- ggplot( d4_long_th, aes(y= Frequency)) +
  geom_point(aes(x=time_jit)) +
  geom_point(data = filter( d4_long_th, vdjdb_TRA_hit == "no" & vdjdb_TRB_hit == "no"), aes(x=time_jit), color = 'grey80') +
  geom_line(data = filter( d4_long_th, vdjdb_TRA_hit == "no" & vdjdb_TRB_hit == "no"),  aes(x=time_jit, group= clone_id), linetype='dashed', color = 'grey80') +
  geom_point(data = filter( d4_long_th, vdjdb_TRA_hit == "yes"), aes(x=time_jit), color = "red" ) + 
  geom_line(data = filter( d4_long_th, vdjdb_TRA_hit == "yes"),  aes(x=time_jit, group= clone_id), linetype='solid', color = 'red') +
  geom_point(data = filter( d4_long_th, vdjdb_TRB_hit == "yes"), aes(x=time_jit), color = "blue" ) + 
  geom_line(data = filter( d4_long_th, vdjdb_TRB_hit == "yes"),  aes(x=time_jit, group= clone_id), linetype='solid', color = 'blue') +
  labs(title = 'Expanding Tfh clones in 321-04', caption = "Red = Tfh clones with cdr3a match in vdj-db\nBlue = Tfh clones with cdr3b match in vdj-db")+
  xlab("Time point") + ylab("Relative frequency\namong Tfh clones") +
  scale_x_continuous(breaks= time_num_df$time_num , 
                     labels= time_num_df$time, 
                     limits=c(min(time_num_df$time_num ) - 1, max(time_num_df$time_num ) + 1 )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45 , size = 8, hjust =  1 , vjust = 1),
        axis.title = element_text( size = 8))

Herati_rel_freq_all <- d5_Tfh_clones_overlay / d4_Tfh_clones_overlay
Herati_rel_freq_th <- d5_Tfh_clones_overlay_th / d4_Tfh_clones_overlay_th
vdj_rel_freq_all <- d5_Tfh_clones_overlay_vdj / d4_Tfh_clones_overlay_vdj
vdj_rel_freq_th <- d5_Tfh_clones_overlay_vdj_th / d4_Tfh_clones_overlay_vdj_th

ggsave("./10x/outs/TwoYear_clonal_tracing_Herati_hits_rel_freq_all.pdf", 
       plot = Herati_rel_freq_all, width = 9.25 , height = 7.8 ,useDingbats= F)
ggsave("./10x/outs/TwoYear_clonal_tracing_vdjdb_hits_rel_freq_all.pdf", 
       plot = vdj_rel_freq_all, width = 9.25 , height = 7.8,useDingbats= F)
ggsave("./10x/outs/TwoYear_clonal_tracing_Herati_hits_rel_freq_topHits.pdf", 
       plot = Herati_rel_freq_th, width = 9.25 , height = 7.8,useDingbats= F)
ggsave("./10x/outs/TwoYear_clonal_tracing_vdjdb_hits_rel_freq_topHits.pdf", 
       plot = vdj_rel_freq_th, width = 9.25 , height = 7.8,useDingbats= F)

ggsave("./10x/outs/TwoYear_clonal_tracing_pickedClones.pdf", 
       plot = d5_Tfh_clones_picked_clones, width = 9.25 , height = 4,useDingbats= F)


# cf summed vs time ====

Tfh_combo <- bind_rows(Tfh_list[[1]], Tfh_list[[2]])

Tfh_cf <- Tfh_combo %>%
  select(donor , time, n) %>% 
  group_by_all() %>% 
  tally(name = 'freq') %>% 
  mutate( naive_cumfreq = freq / sum(freq))  %>%
  filter(n == 1) %>%
  mutate( expanded_cumfreq = 1 - naive_cumfreq ) %>%
  select(-n, -freq) %>%
  gather(3:4, key = 'type', value = 'freq' ) %>%
  ungroup() %>%
  group_by(type) %>%
  group_split()
  
#expanded fraction
cf_expanded <- ggplot( Tfh_cf[[1]], aes(y= 100 * freq, color = Donor)) +
  geom_point(data = Tfh_cf[[1]], aes(x=Day), size = 3) + 
  geom_line(data = Tfh_cf[[1]],  aes(x=Day, group= Donor), linetype='solid', size = 1) +
  scale_x_continuous(breaks=c(0,5,12,28,60,90,120,180), labels=c("Day 0", "Day 5", "Day 12", "Day 28", "Day 60", "Day 90", "Day 120", "Day 180"), limits=c(-2, 182)) +
  theme_minimal() +
  scale_color_jcolors("pal3") +
  xlab("Days post immunization") + ylab("Cumulative frequency (%)") +
  labs(title = 'Cumulative frequency of expanded Tfh clones')+
  theme(axis.text.x = element_text(angle = 45 , size = 8, hjust =  1 , vjust = 1),
        axis.title = element_text( size = 8))

cf_singleton <- ggplot( Tfh_cf[[2]], aes(y= 100 * freq, color = Donor)) +
  geom_point(data = Tfh_cf[[2]], aes(x=Day), size = 3) + 
  geom_line(data = Tfh_cf[[2]],  aes(x=Day, group= Donor), linetype='dashed', size = 1) +
  scale_x_continuous(breaks=c(0,5,12,28,60,90,120,180), labels=c("Day 0", "Day 5", "Day 12", "Day 28", "Day 60", "Day 90", "Day 120", "Day 180"), limits=c(-2, 182)) +
  theme_minimal() +
  scale_color_jcolors("pal3") +
  xlab("Days post immunization") + ylab("Cumulative frequency (%)") +
  labs(title = 'Cumulative frequency of singleton Tfh clones')+
  theme(axis.text.x = element_text(angle = 45 , size = 8, hjust =  1 , vjust = 1),
        axis.title = element_text( size = 8))

cf_plots <- cf_expanded / cf_singleton

ggsave("./10x/outs/BothDonors_Tfh_expanded_singleton_cumfreq.png", plot = cf_plots, width = 9.25 , height = 6)

