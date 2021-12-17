# reformat and plot the tfh in vitro screening data

source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

# function for cleaning the input data files 
clean_files <- function( files ){
  
  # returns a combined dataframe formatted for plotting
  
  files_list <- map(files, ~ read_delim(., delim = '\t'))
  
  short_names <- str_remove(files, "./MG_data/") %>%
    str_remove(., '.txt')
  names(files_list) <- short_names
  
  clean_df_list <- list()
  for ( i in seq_along(files_list)){
    
    df <- files_list[[i]]
    assay_details <- str_split(names(files_list)[i], '_', simplify = T)
    df$readout <- assay_details[,1]
    df$clone <- assay_details[,2]
    df$assay <- assay_details[,3]
    
    colnames(df)[which(colnames(df) %like% "TFH")] <- 'IAV'
    skip_cols <- (ncol(df)-3)
    df_long <- pivot_longer(df, 1:(ncol(df)-3), values_to = 'frequency', names_to = 'stim') %>%
      filter(!is.na(frequency))
    
    df_long$readout <- gsub('TCR','CD3', df_long$readout)
    df_long$stim <- gsub('Mock','Control', df_long$stim )
    
    clean_df_list[[i]] <- df_long
    names(clean_df_list)[i] <- names(files_list)[i]
  }
  combo_df <- do.call(rbind, clean_df_list)
  combo_df$clone <- factor(combo_df$clone, levels = c('Tfh1','Tfh3','Tfh11','Tfh12'))
  
  
  return(combo_df)
  
  
}

# read in the files
IAV_files <- list.files("./MG_data/", "IAV", full.names = T)
trans_files <- list.files("./MG_data/", "transfection", full.names = T)
trans_files_1 <- trans_files[which(str_detect(trans_files,'Tfh11', negate = T))]
peptide_files <- list.files("./MG_data/", "peptide", full.names = T)
Tfh_11_files <- trans_files[which(str_detect(trans_files,'Tfh11'))]

# make df for each
IAV_df <- clean_files(IAV_files)

tranf_df <- clean_files(trans_files_1)
tranf_df$stim <- factor(tranf_df$stim, 
                        levels = c("Control","M","NS","PA","PB1","PB2","NP"))

peptide_df <- clean_files(peptide_files)
pep_levels <- c("Control",
                "NIIGILHLILWILDRLF","LHLILWILDRLFFKCIY","WILDRLFFKCIYRRFKY","LFFKCIYRRFKYGLKGG", 
                "IYRRFKYGLKGGPSTEG", "LKTRPILSPLTKGILGF","ILSPLTKGILGFVFTLT","TKGILGFVFTLTVPSER",
                "GFVFTLTVPSERGLQRR","LTVPSERGLQRRRFVQN", "EQMETDGERQNATEIRA","DGERQNATEIRASVGKM",
                "NATEIRASVGKMIGGIG","RASVGKMIGGIGRFYIQ", "KMIGGIGRFYIQMCTEL","IGRFYIQMCTELKLSDY",
                "IQMCTELKLSDYEGRLI"
)
peptide_df$stim <- factor(peptide_df$stim, levels = pep_levels)



Tfh_11_df <- clean_files(Tfh_11_files)
Tfh_11_df$stim <- factor(Tfh_11_df$stim, levels = 
                           c('unstim',	'HA_60',	'HA_160',	
                             'HA_260',	'HA_360',	'HA_460',	
                             'full_HA'))

# remake the bar plots ====
#color palette
tfh_clone_names <- c('Tfh1','Tfh3','Tfh11', 'Tfh12')
clonePal2 <- c('#943CB4',"#C6CF6E", "#D590DA", "#E48F1B")
names(clonePal2) <- tfh_clone_names

#plotting function
stderror <- function(x) sd(x)/sqrt(length(x))
plotting_func <- function(df, readout_param){
  
  #returns the plotting df, summary stat df, and plot in list structure
  return_list <- list()
  for (i in seq(tfh_clone_names)) {
    
    df2 <- df %>%
      filter(readout == readout_param & 
               clone == tfh_clone_names[i])
    
    df_summary <- df2 %>%
      group_by(clone, stim) %>%
      summarise(mean_freq = mean(frequency),
                sem = stderror(frequency))
    
    
    p <- ggplot() + 
      geom_col(data = df_summary, 
               aes(stim, mean_freq, fill= clone)) +
      geom_errorbar(data = df_summary, 
                    aes(stim, mean_freq, ymin=mean_freq-sem, ymax=mean_freq+sem), 
                    width=.2) +
      geom_point(data = df2, 
                 aes(stim, frequency)) +
      geom_violin() +
      scale_fill_manual(values = clonePal2) +
      theme_minimal() +
      theme(legend.position = 'none',
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
      ylab(paste0("% ", readout_param, "+")) +
      labs(subtitle = tfh_clone_names[i]) 
    
    #pack it up into a list
    return_list[[i]] <- list(df2, df_summary, p)
    names(return_list[[i]]) <- c('readout_df', 'readout_summary', 'plot')
    
  }
  names(return_list) <- tfh_clone_names
  return(return_list)
  
}
add_ttest <- function( plotting_func_results_list, sided = "two.sided" , ref_group = 'Control'){
  
  
  
  for (i in seq_along(plotting_func_results_list)){
    
    plotting_func_results_list[[i]]$plot <- plotting_func_results_list[[i]]$plot + 
      stat_compare_means(data =plotting_func_results_list[[i]]$readout_df, 
                         aes( stim, frequency, label = ..p.signif.. ),
                         ref.group = ref_group, method = 't.test',  hide.ns = TRUE, 
                         method.args = list( alternative =sided))
    
  }
  
  return(plotting_func_results_list)
  
  
}

#make plots

#IAV
IAV_CD69 <- plotting_func(IAV_df, 'CD69')
IAV_CD69 <- add_ttest(IAV_CD69, sided = 'greater')
IAV_CD3 <- plotting_func(IAV_df, 'CD3')
IAV_CD3 <- add_ttest(IAV_CD3, sided = 'less')

IAV_plots <- list(IAV_CD69[[1]]$plot, IAV_CD69[[2]]$plot, 
                  IAV_CD69[[3]]$plot, IAV_CD3[[1]]$plot, 
                  IAV_CD3[[2]]$plot, IAV_CD3[[3]]$plot)%>%
  wrap_plots()
#tranfections
transf_CD3 <- plotting_func(tranf_df, 'CD3')
transf_CD3 <- add_ttest(transf_CD3, sided = 'less')
transf_CD69 <- plotting_func(tranf_df, 'CD69')
transf_CD69 <- add_ttest(transf_CD69, sided = 'greater')

transftransf_plots_p <- list(transf_CD3[[1]]$plot, transf_CD3[[2]]$plot, 
                 transf_CD3[[3]]$plot, transf_CD69[[1]]$plot, 
                 transf_CD69[[2]]$plot, transf_CD69[[3]]$plot) %>%
  wrap_plots()

#peptide
peptide_CD69 <- plotting_func(peptide_df, 'CD69')
peptide_CD69 <- add_ttest(peptide_CD69, sided = 'greater')

peptide_plots <- list(peptide_CD69[[1]]$plot, 
                  peptide_CD69[[2]]$plot, 
                  peptide_CD69[[3]]$plot) %>% 
  wrap_plots()

# Tfh11
Tfh11_CD69 <- plotting_func(Tfh_11_df, 'CD69')
Tfh11_IFNg <- plotting_func(Tfh_11_df, 'IFNg')
Tfh11_plots <- list(Tfh11_CD69$Tfh11$plot, Tfh11_IFNg$Tfh11$plot) %>% 
  wrap_plots()

# save plots 
ggsave('./MG_data/replot/pickedClones_IAV_infection.pdf', plot = IAV_plots,
       width = 7, height = 8, useDingbats = F)
H
ggsave('./MG_data/replot/pickedClones_seg_transfection.pdf', plot = transf_plots,
       width = 8, height = 8, useDingbats = F)

ggsave('./MG_data/replot/pickedClones_peptide_stims.pdf', plot = peptide_plots,
       width = 8, height = 4, useDingbats = F)

ggsave('./MG_data/replot/pickedClones_Tfh11.pdf', plot = Tfh11_plots,
       width = 6, height = 3, useDingbats = F)
