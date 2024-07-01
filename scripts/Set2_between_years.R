# TwoYear Tfh phenotypes over time

source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

rename <- dplyr::rename
# import T cells Seurat object and df ----
Tfh_obj <- readRDS(Set2_integrated_Tfh_path)
Clonedf <- read.delim( clone_df_path , stringsAsFactors = F) 

GC_TFH <- subset(Tfh_obj, Tfh_type == 'GC' & donor == '321-05' & tissue =='FNA')


earlyMarks <- FindMarkers(GC_TFH, ident.1 = 'y1_d5', ident.2 = 'y2_d7', group.by = 'time_point')

comp1_mark <- earlyMarks %>%
  filter(rownames(.) %notin% IgTcr$genes) %>%
  arrange(desc(avg_log2FC )) %>%
  mutate( across(c('p_val','p_val_adj'), ~log10(.))) %>%
  mutate( across(c('p_val','p_val_adj'), abs)) %>%
  mutate( time_point = ifelse(avg_log2FC  > 0 , 'y1_d5', "y2_d7")) %>%
  mutate( gene = rownames(.))


comp1_vp <- ggplot(comp1_mark, aes(avg_log2FC, p_val_adj)) + 
  geom_point(aes(color = time_point)) + 
  theme_minimal() +
  geom_text_repel( 
    data=  filter(comp1_mark, p_val_adj > 10 & abs(avg_log2FC) > 0.5 & grepl('^RP[LS][0-9]', gene) == F), 
    aes(label=gene),
    size = 3,
    nudge_x = .15,
    box.padding = 0.5,
    nudge_y = 2,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20,
    max.overlaps = 15
  ) + 
  ylab("-log10( adj p-value )") + 
  xlab("Year 1 day 5 vs. Year 2 day 7\n( avg log2 fold-change )") +
  coord_cartesian(clip = 'off') +
  scale_color_manual(values = TimePal2)+
  labs(color = 'Time Point')





yearMarks <- FindMarkers(GC_TFH, ident.1 = 1, ident.2 = 2, group.by = 'year')

comp2_mark <- yearMarks %>%
  filter(rownames(.) %notin% IgTcr$genes) %>%
  arrange(desc(avg_log2FC )) %>%
  mutate( across(c('p_val','p_val_adj'), ~log10(.))) %>%
  mutate( across(c('p_val','p_val_adj'), abs)) %>%
  mutate( year = ifelse(avg_log2FC  > 0 , 1, 2)) %>%
  mutate( gene = rownames(.)) 

comp2_vp <- ggplot(comp2_mark, aes(avg_log2FC, p_val_adj)) + 
  geom_point(aes(color = factor(year))) + 
  theme_minimal() +
  geom_text_repel( 
    data=  filter(comp2_mark, p_val_adj > 10 & abs(avg_log2FC) > 0.4 & grepl('^RP[LS][0-9]', gene) == F), 
    aes(label=gene),
    size = 3,
    nudge_x = .15,
    box.padding = 0.5,
    nudge_y = 2,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20,
    max.overlaps = 15
  ) + 
  ylab("-log10( adj p-value )") + 
  xlab("Year 1 vs. Year 2\n( avg log2 fold-change )") +
  coord_cartesian(clip = 'off') +
  scale_color_manual(values = YearPal)+
  labs(color = 'Year')




monthMarks <- FindMarkers(GC_TFH, ident.1 = 'y1_d28', ident.2 = 'y2_d28', group.by = 'time_point')

comp3_mark <- monthMarks %>%
  filter(rownames(.) %notin% IgTcr$genes) %>%
  arrange(desc(avg_log2FC )) %>%
  mutate( across(c('p_val','p_val_adj'), ~log10(.))) %>%
  mutate( across(c('p_val','p_val_adj'), abs)) %>%
  mutate( time_point = ifelse(avg_log2FC  > 0 , 'y1_d28', "y2_d28")) %>%
  mutate( gene = rownames(.))


comp3_vp <- ggplot(comp3_mark, aes(avg_log2FC, p_val_adj)) + 
  geom_point(aes(color = time_point)) + 
  theme_minimal() +
  geom_text_repel( 
    data=  filter(comp3_mark, p_val_adj > 10 & abs(avg_log2FC) > 0.5 & grepl('^RP[LS][0-9]', gene) == F), 
    aes(label=gene),
    size = 3,
    nudge_x = .15,
    box.padding = 0.5,
    nudge_y = 2,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20,
    max.overlaps = 15
  ) + 
  ylab("-log10( adj p-value )") + 
  xlab("Year 1 day 28 vs. Year 2 day 28\n( avg log2 fold-change )") +
  coord_cartesian(clip = 'off') +
  scale_color_manual(values = TimePal2)+
  labs(color = 'Time Point')



wk1 <- head(comp1_mark,n = 10)$gene
wk2 <- tail(comp1_mark,n = 10)$gene
mo1 <- head(comp3_mark,n = 10)$gene
mo2 <- tail(comp3_mark,n = 10)$gene
yr1 <- head(comp2_mark,n = 10)$gene
yr2 <- tail(comp2_mark,n = 10)$gene
#marks <- c(yr1,yr2)[which(grepl('^RP[LS][0-9]', c(yr1,yr2)) == F)]
marks <-c(wk1, wk2, mo1, mo2, yr1,yr2)
lin <- c('CXCR5','PDCD1','TOX2','IL21','CXCL13')
marks <-c(marks, lin)



firstdays <- subset(GC_TFH, time_point %in% c('y1_d0','y1_d5','y1_d28','y2_d0','y2_d7','y2_d28'))
firstdays@meta.data$time_point <-  factor(firstdays@meta.data$time_point,levels =c('y1_d0','y1_d5','y1_d28','y2_d0','y2_d7','y2_d28'))

hm <- DoHeatmap(firstdays, marks, group.by = 'time_point')  + 
  scale_fill_viridis(option = 'A')

widedf <- hm$data %>%
  filter(!is.na(Expression)) %>%
  pivot_wider( names_from = Feature, values_from = Expression) %>%
  arrange(Identity)

mat <- as.matrix(widedf[,3:ncol(widedf)])
rownames(mat) <- widedf$Cell 
mat<- t(mat)
annmeta <- data.frame(time_point = widedf$Identity, row.names = widedf$Cell)

annmeta$year <- as.factor(ifelse(grepl('y1', annmeta$time_point), 1,2))

breakhm <- nrow(filter(annmeta, year ==1))
TimePal2c <- TimePal2[which(names(TimePal2) %in%annmeta$time_point )]
color_list = list(time_point = TimePal2c, year = YearPal)

hm2 <- pheatmap(mat, scale = 'none',
         show_colnames = F, 
         annotation_col = annmeta,
         annotation_colors =  color_list,
         fontsize_row = 8,
         cluster_cols = F,
         gaps_col = breakhm ,color = viridisLite::inferno(50)) %>% 
  as.ggplot()



ggsave('./10x/outs/Set2_GC_crossyear_markers_hm.png', 
       plot = hm2, width = 6.5, height = 8)
ggsave('./10x/outs/Set2_GC_crossyear_markers_hm.pdf', 
       plot = hm2, width = 6.5, height = 8)
#Violin of IL21
vln <- scCustomize::Stacked_VlnPlot(GC_TFH, 
                                    c('IL21','CXCL13','CXCR5'),
                                    group.by = 'time_point',
                                    x_lab_rotate = T,pt.size = 0.2)+
  plot_annotation('Donor 321-05')

ggsave('./10x/outs/Set2_d05_GC_crossyear_cyto_vln.png', 
       plot = vln, width = 6.5, height = 3)

ggsave('./10x/outs/Set2_d05_GC_crossyear_cyto_vln.pdf', 
       plot = vln, width = 6.5, height = 3)

panel1 <- (comp1_vp / comp3_vp/ comp2_vp) 

ggsave('./10x/outs/Set2_GC_crossyear_volcano.png', 
       plot = panel1,  width = 5, height = 8)

ggsave('./10x/outs/Set2_GC_crossyear_volcano.pdf', 
       plot = panel1,  width = 5, height = 8)





GC_TFH2 <- subset(Tfh_obj, Tfh_type == 'GC' & 
                    tissue =='FNA' & 
                    donor %in% c('321-04'))

vln2 <- scCustomize::Stacked_VlnPlot(GC_TFH2, 
                                    c('IL21','CXCL13','CXCR5'),
                                    group.by = 'time_point',
                                    x_lab_rotate = T,pt.size = 0.2) +
  plot_annotation('Donor 321-04')
ggsave('./10x/outs/Set2_d04_GC_crossyear_cyto_vln.png', 
       plot = vln2, width = 6.5, height = 3)
ggsave('./10x/outs/Set2_d04_GC_crossyear_cyto_vln.pdf', 
       plot = vln2, width = 6.5, height = 3)








df <- FetchData(GC_TFH, c('time_point','IL21','CXCL13')) %>%
  filter(time_point %in% c('y1_d5','y1_d28','y2_d7','y2_d28')) %>%
  pivot_longer(IL21:CXCL13, values_to = 'exp',names_to = 'gene')


ggplot(df, aes(time_point, exp))+
  geom_boxplot()+
  facet_wrap(~gene)+
  stat_compare_means(comparisons = list(c('y1_d5','y2_d7'),c('y1_d28','y2_d28')))
