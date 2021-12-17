# TwoYear Tfh phenotypes

source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

# import T cells Seurat object and df ----
Tfh_obj <- readRDS(Tfh_lineages_path)
Clonedf <- read.delim(clone_df_path, stringsAsFactors = F)

# Find markers between Tfh types ====

Tfh_type <- FetchData( Tfh_obj, vars = c('Tfh_type'))
Tfr_cells <- Tfh_type %>% filter(Tfh_type == 'IL10 TFH') %>% rownames(.)
GC_cells <- Tfh_type %>% filter(Tfh_type == 'GC') %>% rownames(.)
EM_cells <- Tfh_type %>% filter(Tfh_type == 'pre/memory') %>% rownames(.)

Tfr_vs_GC <- FindMarkers(Tfh_obj, ident.1 =Tfr_cells, ident.2 = GC_cells) 
Tfr_vs_GC$gene <- rownames(Tfr_vs_GC)
write_tsv( Tfr_vs_GC, './10x/outs/TwoYear_Tfh_lineages_IL10_TFH_vs_GC_markers.tsv')
Tfr_vs_EM <- FindMarkers(Tfh_obj, ident.1 = Tfr_cells, ident.2 = EM_cells)
Tfr_vs_EM$gene <- rownames(Tfr_vs_EM)
write_tsv( Tfr_vs_EM, './10x/outs/TwoYear_Tfh_lineages_IL10_TFH_vs_premem_markers.tsv')
EM_vs_GC <- FindMarkers(Tfh_obj, ident.1 =EM_cells, ident.2 = GC_cells) 
EM_vs_GC$gene <- rownames(EM_vs_GC)
write_tsv( EM_vs_GC, './10x/outs/TwoYear_Tfh_lineages_premem_vs_GC_markers.tsv')

Tfh_obj <- SetIdent(Tfh_obj, value = 'Tfh_type')
Tfh_type_markers <- FindAllMarkers(Tfh_obj)
Tfh_obj <- SetIdent(Tfh_obj, value = 'ident')
write_tsv( Tfh_type_markers, './10x/outs/TwoYear_Tfh_lineages_Tfh_type_markers.tsv')

Tfh_type_markers_up <- Tfh_type_markers %>%
  filter(avg_log2FC > 0 ) %>%
  group_by(cluster) %>% 
  arrange(cluster, desc(avg_log2FC))
write_tsv( Tfh_type_markers_up, './10x/outs/TwoYear_Tfh_lineages_Tfh_type_up_markers.tsv')


Tfr_vs_GC2 <- Tfr_vs_GC %>%
  arrange(desc(avg_log2FC )) %>%
  mutate( across(c('p_val','p_val_adj'), log10)) %>%
  mutate( across(c('p_val','p_val_adj'), abs)) %>%
  mutate( Tfh_type = ifelse(avg_log2FC  > 0 , 'IL10 TFH', "GC")) %>%
  mutate( gene = rownames(.)) 

Tfr_vs_EM2 <- Tfr_vs_EM %>%
  arrange(desc(avg_log2FC )) %>%
  mutate( across(c('p_val','p_val_adj'), log10)) %>%
  mutate( across(c('p_val','p_val_adj'), abs)) %>%
  mutate( Tfh_type = ifelse(avg_log2FC  > 0 , 'IL10 TFH', "pre/memory")) %>%
  mutate( gene = rownames(.)) 

EM_vs_GC2 <- EM_vs_GC %>%
  arrange(desc(avg_log2FC )) %>%
  mutate( across(c('p_val','p_val_adj'), log10)) %>%
  mutate( across(c('p_val','p_val_adj'), abs)) %>%
  mutate( Tfh_type = ifelse(avg_log2FC  > 0 , 'pre/memory', "GC")) %>%
  mutate( gene = rownames(.)) 


# volcano plot
Tfr_vs_GC_vp <- ggplot(Tfr_vs_GC2, aes(avg_log2FC, p_val_adj)) + 
  geom_point(aes(color = Tfh_type )) + 
  theme_minimal() +
  scale_color_manual(values = TfhPal) +
  geom_text_repel( 
    data=  filter(Tfr_vs_GC2, p_val_adj > 50 | abs(avg_log2FC) > 1.5 ), 
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
  xlab("IL10 Tfh vs. GC Tfh expression ( avg log2 fold-change )") +
  theme(axis.text = element_text(size = 10)) +
  NoLegend()


Tfr_vs_EM_vp <- ggplot(Tfr_vs_EM2, aes(avg_log2FC, p_val_adj)) + 
  geom_point(aes(color = Tfh_type )) + 
  theme_minimal() +
  scale_color_manual(values = TfhPal) +
  geom_text_repel( 
    data=  filter(Tfr_vs_EM2, p_val_adj > 50 | abs(avg_log2FC) > 1.5 ), 
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
  xlab("IL10 Tfh vs. pre/memory Tfh expression ( avg log2 fold-change )") +
  theme(axis.text = element_text(size = 10)) +
  NoLegend()

EM_vs_GC_vp <- ggplot(EM_vs_GC2, aes(avg_log2FC, p_val_adj)) + 
  geom_point(aes(color = Tfh_type )) + 
  theme_minimal() +
  scale_color_manual(values = TfhPal) +
  geom_text_repel( 
    data=  filter(EM_vs_GC2, p_val_adj > 50 | abs(avg_log2FC) > 1.5 ), 
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
  xlab("pre/memory Tfh vs. GC Tfh expression ( avg log2 fold-change )") +
  theme(axis.text = element_text(size = 10)) +
  NoLegend()


vol_plots <- EM_vs_GC_vp | Tfr_vs_EM_vp | Tfr_vs_GC_vp

ggsave('./10x/outs/TwoYear_Tfh_lineages_Tfh_type_volcano.pdf', vol_plots, 
       height = 5, width = 15, useDingbats = F)

# Tfh_only umap scatter plots =====
cpal <- cols25(12)
names(cpal) <- seq(12) - 1

ident_umap_Tfh <- DimPlot(Tfh_obj, pt.size = 2, ) + 
  scale_color_manual(values = cpal) + 
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        plot.subtitle = element_text(size = 20)) + 
  labs( subtitle = 'GEX cluster' ) +
  ggtitle(label = NULL)

tissue_umap_Tfh <- DimPlot(Tfh_obj, group.by = "tissue", order = T, pt.size = 2) + 
  scale_color_manual(values = TissuePal) + 
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        plot.subtitle = element_text(size = 20))+ 
  labs( subtitle = 'Tissue origin' ) +
  ggtitle(label = NULL)

type_umap_Tfh <- DimPlot(Tfh_obj, group.by = 'Tfh_type', pt.size = 2) + 
  scale_color_manual(values = TfhPal) + 
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        plot.subtitle = element_text(size = 20))+ 
  labs( subtitle = 'T cell type') +
  ggtitle(label = NULL)

donor_umap_Tfh <- DimPlot(Tfh_obj, group.by = "donor", pt.size = 2) + 
  scale_color_manual(values = DonorPal[1:2]) + 
  theme( axis.line = element_blank(),
         axis.title = element_blank(),
         plot.subtitle = element_text(size = 20)) +
  labs( subtitle = 'Donor') +
  ggtitle(label = NULL)

picked_clone_umap_Tfh <- DimPlot(Tfh_obj, group.by = "Tfh_clone_id", pt.size = 2, order = T) + 
  scale_color_manual(values = clonePal) + 
  theme( axis.line = element_blank(),
         axis.title = element_blank(),
         plot.subtitle = element_text(size = 20)) +
  labs( subtitle = 'Clonotypes picked for screening') +
  ggtitle(label = NULL)

year_umap_Tfh <- DimPlot(Tfh_obj, group.by = "year", pt.size = 2, order = T) + 
  scale_color_manual(values = YearPal) + 
  theme( axis.line = element_blank(),
         axis.title = element_blank(),
         plot.subtitle = element_text(size = 20)) +
  labs( subtitle = 'Study Year') +
  ggtitle(label = NULL)

umaps_Tfh <- (ident_umap_Tfh / type_umap_Tfh) 

# Tfh_only supp umaps
TimeUMAP_Tfh <- DimPlot(Tfh_obj, group.by = "time", pt.size = 2) + 
  scale_color_manual(values = TimePal2) + 
  theme( axis.line = element_blank(),
         axis.title = element_blank(),
         plot.subtitle = element_text(size = 20)) + 
  labs(subtitle = "Days post-immunization")+
  theme(plot.subtitle = element_text(size = 20)) +
  ggtitle(label = NULL)

TimeUMAP_Tfh2 <- DimPlot(Tfh_obj, group.by = "day", pt.size = 2, split.by = 'year') + 
  scale_color_manual(values = TimePal) + 
  theme( axis.line = element_blank(),
         axis.title = element_blank(),
         plot.subtitle = element_text(size = 20)) + 
  labs(subtitle = "Days post-immunization")+
  theme(plot.subtitle = element_text(size = 20)) +
  ggtitle(label = NULL)

CloneUMAP_Tfh <- FeaturePlot(Tfh_obj, features = 'clone_size_log', pt.size = 2 ) +
  theme( axis.line = element_blank(),
         axis.title = element_blank(),
         plot.subtitle = element_text(size = 20,hjust = 0)) +
  labs(subtitle = 'Clone size') +
  scale_color_viridis(option = 'C') +
  ggtitle(label = NULL)

supp_umaps_Tfh1 <- (tissue_umap_Tfh | donor_umap_Tfh) / (TimeUMAP_Tfh | year_umap_Tfh) / (picked_clone_umap_Tfh | CloneUMAP_Tfh)

ggsave('./10x/outs/TwoYear_Tfh_lineages_supp_umaps_v1.png', supp_umaps_Tfh1, height = 15, width = 15)

supp_umaps_Tfh2 <- (tissue_umap_Tfh | donor_umap_Tfh) / (TimeUMAP_Tfh2 | year_umap_Tfh) / (picked_clone_umap_Tfh | CloneUMAP_Tfh)

ggsave('./10x/outs/TwoYear_Tfh_lineages_supp_umaps_v1.png', supp_umaps_Tfh2, height = 15, width = 15)

#figure umap
figs_umaps_Tfh <- ident_umap_Tfh | donor_umap_Tfh | tissue_umap_Tfh | TimeUMAP_Tfh 
figs_umaps_Tfh_v2 <- ident_umap_Tfh | donor_umap_Tfh | tissue_umap_Tfh | TimeUMAP_Tfh2 

ggsave('./10x/outs/TwoYear_Tfh_lineages_figure_umaps.png', figs_umaps_Tfh, height = 5, width = 20)
ggsave('./10x/outs/TwoYear_Tfh_lineages_figure_umaps_v2.png', figs_umaps_Tfh_v2, height = 5, width = 20)

type_umap_Tfh2 <- type_umap_Tfh + theme(legend.position = 'top')
ggsave('./10x/outs/TwoYear_Tfh_lineages_Tfh_type.png',type_umap_Tfh2 , height = 5, width = 6)

# splitting by donor to make some fancier heatmaps and violins ====

tfh_gene_set_c <- c('CXCR5', 'PDCD1','CCR7','SELL',
                   'VIM', 'TMSB10', 'EMP3','TCF7',
                    'TOX2','GNG4', 'CXCL13', 'IL21',
                    'FOXB1','DUSP4', 'DDIT4', 'IL10')
ViolinPlots <- list()
for (i in tfh_gene_set_c){
  vp <- VlnPlot(Tfh_obj, features = i, group.by = 'Tfh_type', split.by = 'donor', combine = T) +
    theme(axis.title = element_blank()) + 
    scale_fill_manual(values = DonorPal[1:2])
  ViolinPlots <- list.append(ViolinPlots, vp)
}

bigVP <- wrap_plots(ViolinPlots, ncol = 4) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')

ggsave('./10x/outs/TwoYear_Tfh_lineages_Tfh_type_markers_violin.pdf', 
       bigVP, width = 1143/72, height = 861/72, useDingbats = F)

# DotPlot 


plot_feats <- list('Adhesion / Migration' = c('CXCR4','CXCR5', 'CXCR6','CCR7','SELL', 'VIM', 'TMSB10', 'EMP3', 'CADM1'),
                   'Cytokines' = c('IL21','CXCL13', 'IL10'),
                   'Co-stimulation' = c('PDCD1','ICOS','LAG3','TNFRSF18'),
                   'Transcription Factors' = c('TCF7','BCL6','TOX2','NR4A2','FOXB1','PRDM1','IKZF3','NFATC1'),
                   'Signaling' = c('DDIT4','RGS1','MAP3K8','DUSP4')
)

pfeat2 <- unlist(plot_feats) 
names(pfeat2) <-NULL
p <- DotPlot(Tfh_obj, features = pfeat2, 
             group.by = 'Tfh_type',
             split.by = 'donor', dot.scale = 10 , cols = 'RdYlBu') +
  theme(axis.text.x = element_text(angle =90, vjust = 0.5, hjust = 1)) +
  xlab('') + ylab('') 
yorder <- p$data$id %>% unique() %>% as.character()
yorder <- yorder[c(4,3,2,1,6,5)] %>%rev()
p <- p + scale_y_discrete(limits = yorder)

ggsave('./10x/outs/TwoYear_Tfh_lineages_Tfh_type_markers_DotPlot.pdf', 
       p, width = 12, height = 4.1, useDingbats = F)


#save ====
saveRDS(Tfh_obj, Tfh_lineages_path )




