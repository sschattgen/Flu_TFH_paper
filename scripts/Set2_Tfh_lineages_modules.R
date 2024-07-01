#module score heatmap
source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
library(ggpubr)
#
Tfh_obj <- readRDS(Set2_integrated_Tfh_path)
Tfh_obj@meta.data$day <- factor(Tfh_obj@meta.data$day, levels = names(TimePal))
Tfh_obj@meta.data$time_point <- factor(Tfh_obj@meta.data$time, levels = names(TimePal2))

Clonedf <- read.delim( clone_df_path , stringsAsFactors = F)

flu_pal <- c('grey90','red')
names(flu_pal) <- c('unknown','yes')

for (donor1 in c('321-05','321-04')){
  
  Tfh_obj2 <- subset(Tfh_obj, 
                     Tfh_type %notin% c('naive','Treg') & 
                       donor == donor1 & 
                       tissue =='FNA')
  
  ann_colors = list(
    'Clone ID' = clonePal,
    tissue = TissuePal,
    'Flu specificity' = flu_pal,
    'Time point' = TimePal2,
    Tfh_type =TfhPal
  )
  
  # annotation metadata
  sample_col <- FetchData(Tfh_obj2, c('Tfh_clone_id', 'flu_specific' ,'time_point', 'tissue','Tfh_type'))%>%
    arrange( time_point , Tfh_clone_id)
  
  colnames(sample_col)[1:3] <- c('Clone ID', 'Flu specificity','Time point')
  sample_md <- sample_col %>% 
    select(`Clone ID`, `Flu specificity`, `Time point`, Tfh_type, tissue)
  
  
  mod_df <- FetchData(Tfh_obj2, c('premem_TFH_module','GC_TFH_module','IL10_TFH_module', 'time_point')) %>% 
    arrange(time_point) %>%
    select(-time_point) %>%
    as.matrix() %>%
    t()
  
  
  pc_hm2 <- pheatmap(mod_df, 
                     color = viridis::viridis(50, option = "B"),
                     cluster_cols = FALSE,
                     cluster_rows = FALSE,
                     annotation_col = sample_md, 
                     fontsize_row = 10, border_color= NA,
                     annotation_colors = ann_colors,
                     show_colnames  = FALSE) %>%
    as.ggplot(.)
  
  ggsave(paste0('./10x/outs/Set2_Tfh_lineages_module_hm_',donor1,'.pdf'), plot = pc_hm2,
         height = 8, width = 17)
  
  mod_df2 <- FetchData(Tfh_obj2, c('premem_TFH_module','GC_TFH_module','IL10_TFH_module', 'time_point')) %>%
    pivot_longer(1:3, names_to = 'module', values_to = 'score') %>%
    mutate(module = gsub('_TFH_module','', module)) %>%
    mutate(module = gsub('premem','pre/memory', module)) %>%
    mutate(module = gsub('IL10','IL10 TFH', module))
  
  p1 <- ggplot(mod_df2, aes(time_point,score, fill = module))+
    geom_boxplot(outlier.shape = NA, alpha = 0.7)+
    geom_jitter(alpha = 0.1) +
    facet_wrap(~module, ncol = 1, scales = 'free_y') +
    scale_fill_manual(values = TfhPal) + 
    theme_bw() +
    ylab('Module Score') +
    xlab('Time Point') +
    labs(title = donor1)
  
  width = (length(unique(mod_df2$time_point))*0.667) +1
  
  ggsave(paste0('./10x/outs/Set2_Tfh_lineages_module_boxplots_',donor1,'.pdf'), plot = p1,
         height = 5, width = width, useDingbats = F)
  
  
  
}

Tfh_obj3 <- subset(Tfh_obj, tissue !='PBMC')
mod_df3 <- FetchData(Tfh_obj3, c('premem_TFH_module','GC_TFH_module','IL10_TFH_module', 'time_point','donor')) %>%
  pivot_longer(1:3, names_to = 'module', values_to = 'score') %>%
  mutate(module = gsub('_TFH_module','', module)) %>%
  mutate(module = gsub('premem','pre/memory', module)) %>%
  mutate(module = gsub('IL10','IL10 TFH', module))

width = (length(unique(mod_df3$time_point))*0.667) +1

p2 <- ggplot(mod_df3, aes(time_point,score, fill = module, color = donor))+
  geom_boxplot(outlier.shape = NA)+
  #geom_jitter(alpha = 0.1) +
  facet_wrap(~module, ncol = 1, scales = 'free_y') +
  scale_fill_manual(values = TfhPal) + 
  scale_color_manual(values = DonorPal) +
  theme_bw() +
  ylab('Module Score') +
  xlab('Time Point') 

ggsave('./10x/outs/Set2_Tfh_lineages_module_boxplots_allDonors.pdf', plot = p2,
       height = 5, width = width, useDingbats = F)





mod_df4 <- mod_df3 %>%
  filter(time_point %in% c('y1_d0','y1_d5','y2_d0','y2_d7')) %>%
  mutate(time_point = factor(time_point, levels = names(TimePal2))) %>%
  mutate(year = gsub('y','',str_split(time_point,pattern = '_',simplify = T)[,1])) %>%
  mutate(time = ifelse(grepl('d0',time_point),'pre','post')) 


p3 <- ggplot(mod_df4, aes(time_point, score, color = year, fill = time_point))+
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means(aes(label = after_stat(p.signif)),label.y.npc = 0,
                     comparisons = list(c('y1_d0','y1_d5'),c('y2_d0','y2_d7'), c('y1_d0','y2_d0'), c('y1_d5','y2_d7')),
                     tip.length = 0.01,vjust = 0.4)  +
  scale_fill_manual(values = TimePal2) + 
  scale_color_manual(values = YearPal) + 
  scale_y_continuous(expand = c(0,0))+
  theme_bw() +
  ylab('Module Score') +
  xlab('Year') +
  labs(fill = 'Time Point', color = 'Year')+
  #annotate("rect", xmin = 2.5, xmax = 5, ymin =-0.5, ymax = 1,alpha = .1,fill = YearPal[2])+ 
  #annotate("rect", xmin = 0, xmax = 2.5, ymin =-0.5, ymax = 1,alpha = .1,fill = YearPal[1])+
  facet_wrap(~module, nrow = 3) +
  scale_y_continuous(limits = c(-.5, 1.6))

ggsave('./10x/outs/Set2_Tfh_lineages_module_boxplots_321-05_early_times.pdf', plot = p3,
       height = 6, width = 3.5, useDingbats = F)






library("rjson")


geneset_files <- c('https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/glycolysis+and+gluconeogenesis/KEGG+Pathways',
                   'https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/oxidative+phosphorylation/KEGG+Pathways',
                   'https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/Cellular+response+to+hypoxia/Reactome+Pathways')


geneset_df <- map(geneset_files, ~fromJSON(paste(readLines(.), collapse="")))

genes_lists <- list(map_chr(geneset_df[[1]]$associations, ~.$gene$symbol),
                    map_chr(geneset_df[[2]]$associations, ~.$gene$symbol),
                    map_chr(geneset_df[[3]]$associations, ~.$gene$symbol))

double_genes <- intersect(genes_lists[[1]], genes_lists[[2]])
genes_lists <-map(genes_lists, ~.[which(. %notin%double_genes & nchar(.) != 0)])
names(genes_lists) <- c('glycolysis','oxphos','hypoxia')



# 
# geneset_files <- c('./scripts/Set2/KEGG_GLYCOLYSIS_GLUCONEOGENESIS.v2023.2.Hs.tsv',
#                    './scripts/Set2/KEGG_OXIDATIVE_PHOSPHORYLATION.v2023.2.Hs.tsv')
# geneset_df <- map(geneset_files, ~read.delim(.))
# genes_lists <-map(geneset_df, ~str_split(.[17,2], ',') %>% unlist() %>% unique())
# double_genes <- intersect(genes_lists[[1]], genes_lists[[2]])
# genes_lists <-map(genes_lists, ~.[which(. %notin%double_genes & nchar(.) != 0)])
# names(genes_lists) <- c('glycolysis','oxphos')
# 

# Fu_gs_all <- readxl::read_excel('./scripts/Fu_Nature_2021_Tfh_metabolism.xlsx')
# Fu_gs_Tfh <- Fu_gs_all[ which(Fu_gs_all[,4] < 0 & Fu_gs_all[,5] < 0.05  ) ,]
# Tfh_str8_match <- str_to_upper(pull(Fu_gs_Tfh[,2]))
# Tfh_str8_match <- Tfh_str8_match[which(Tfh_str8_match %in% rownames(Tfh_obj))]
# 
# genes_lists[['FAO']] <- Tfh_str8_match

Tfh_obj <- AddModuleScore(Tfh_obj, features = genes_lists ,name =names(genes_lists),  nbin =16,assay = 'RNA')
colnames(Tfh_obj@meta.data)[which(colnames(Tfh_obj@meta.data) == 'glycolysis1')] <- 'glycolysis_module'
colnames(Tfh_obj@meta.data)[which(colnames(Tfh_obj@meta.data) == 'oxphos2')] <- 'oxphos_module'
colnames(Tfh_obj@meta.data)[which(colnames(Tfh_obj@meta.data) == 'hypoxia3')] <- 'hypoxia_module'


df <- FetchData(Tfh_obj, c('donor','tissue','time_point','Tfh_type','glycolysis_module','oxphos_module','hypoxia_module'))


df2 <- df %>%
  filter(donor %in% c('321-04','321-05') &
           tissue =='FNA' & 
           Tfh_type %notin% c('naive','Treg')
         ) %>%
  pivot_longer(glycolysis_module:hypoxia_module, names_to = 'module',values_to = 'score')


median_df <- df2 %>%
  group_by(donor, time_point, module) %>%
  summarise(median = median(score))


time_breaks <- c(0,5,12,28,60,90,120,180,220,227,234,248,280,310,340)
tdf <- data.frame(time_point = factor(names(TimePal2)[c(1:3,5:length(names(TimePal2)))], levels = names(TimePal2)[c(1:3,5:length(names(TimePal2)))]),
                  time_encode = as.integer(time_breaks))

median_df <- left_join(median_df, tdf) 
df2 <- left_join(df2, tdf) 

donor_median <- median_df %>% 
  group_by(donor) %>% 
  group_split()
donor_all <- df2 %>% 
  group_by(donor) %>% 
  group_split()

d4_ <- ggplot()+
  geom_boxplot(data = donor_all[[1]], aes(time_point, score, color = module)) +
  geom_point(data = donor_median[[1]], aes(time_point, median, color = module), size = 2) +
  geom_path(data = donor_median[[1]], aes(time_point, median, color = module, group = module), linewidth = 1) +
  #scale_x_continuous( breaks=time_breaks[which(time_breaks %in% donor_median[[1]]$time_encode)], 
  #                    labels=c("0", "5", "12", "120", "0","14","90","120"), 
  #                    limits=c(-2, 342)) +
  theme_minimal() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.title.x = element_blank()) +
  xlab("Days post vaccination") +
  ylab("Gene set module score") +
  scale_color_brewer(palette = "Set1") 

d5_ <- ggplot()+
  geom_boxplot(data = donor_all[[2]], aes(time_point, score, color = module)) +
  geom_point(data = donor_median[[2]], aes(time_point, median, color = module), size = 2) +
  geom_path(data = donor_median[[2]], aes(time_point, median, color = module, group = module), linewidth = 1) +
  #scale_x_continuous( breaks=time_breaks[which(time_breaks %in% donor_median[[2]]$time_encode)], 
  #                    labels=c("0", "5", "12","28","60","90","180", "0","7","28","60","120"), 
  #                    limits=c(-2, 342)) +
  theme_minimal() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.title.x = element_blank()) +
  xlab("Days post vaccination") +
  ylab("Gene set module score") +
  scale_color_brewer(palette = "Set1") 

mod_df4 <- donor_all[[2]] %>%
  filter(time_point %in% c('y1_d0','y1_d5','y2_d0','y2_d7')) %>%
  mutate(time_point = factor(time_point, levels = names(TimePal2))) %>%
  mutate(year = gsub('y','',str_split(time_point,pattern = '_',simplify = T)[,1])) %>%
  mutate(time = ifelse(grepl('d0',time_point),'pre','post'))


ggplot(mod_df4, aes(year, score, color = year, fill = time_point))+
  geom_boxplot(outlier.shape = NA)+
  #geom_jitter(alpha = 0.1) +
  facet_wrap(~module, ncol = 1, scales = 'free_y') +
  scale_fill_manual(values = TimePal2) + 
  scale_color_manual(values = YearPal) + 
  theme_bw() +
  ylab('Module Score') +
  xlab('Year') +
  labs(fill = 'Time Point', color = 'Year')

ggplot()+
  geom_boxplot(data = donor_all[[2]], aes(time_point, score, color = module)) +
  theme_minimal() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.title.x = element_blank()) +
  xlab("Days post vaccination") +
  ylab("Gene set module score") +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~Tfh_type)

Tfh_obj2 <- subset(Tfh_obj, 
                   Tfh_type %notin% c('naive','Treg') & 
                     donor %in% c('321-04','321-05') &
                     tissue =='FNA')
mat <- AverageExpression(Tfh_obj,assays = 'SCT', features = genes_lists$hypoxia, group.by = 'time_point')$SCT 
mat <- mat[rowSums(mat) != 0 , ]

pheatmap(mat,cluster_cols = F)



GOI_list <- list(
  oxphos = c('MT-ATP6','MT-ND3','MT-ND4L'), 
  glycolysis = c('LDHB','GAPDH','HK1'))
plist <- list()
for (i in seq_along(GOI_list)){
  
  
  in_name <- names(GOI_list)[i]
  
  oxp_df <- FetchData(Tfh_obj2, c(GOI_list[[i]],'time_point','donor')) %>%
    pivot_longer(1:length(GOI_list[[i]]), names_to = 'gene', values_to = 'expression')
  
  median_df <- oxp_df %>%
    group_by(donor, time_point, gene) %>%
    summarise(median = median(expression))
  
  plist[[in_name]] <- ggplot(oxp_df, aes(time_point, expression, fill = gene)) +
    geom_boxplot() +
    geom_point(data = median_df, aes(time_point, median, color = gene), size = 2) +
    geom_path(data = median_df, aes(time_point, median, color = gene, group = gene), linewidth = 1) +
    facet_wrap(~donor, ncol = 1) +
    theme_bw() +
    ggtitle(in_name) +
    theme(legend.position = 'bottom')
  
}

wrap_plots(plist, nrow = 1)
