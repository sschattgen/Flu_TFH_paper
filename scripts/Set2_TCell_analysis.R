
source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
library(scCustomize)
# imports ====

Tcells <- readRDS(Set2_integrated_Tcells_path)
TMarkers <- read.csv( './10x/outs/integrated_markers.csv' , stringsAsFactors = F)
Clonedf <- read.delim(clone_df_path, stringsAsFactors = F)
Tcells@meta.data$time_point <- factor(Tcells@meta.data$time_point, levels = names(TimePal2)) 
top_DEGs <- Extract_Top_Markers(TMarkers, num_genes = 5) %>% unique()

# DotPlots of clusters ====

basic <-c("CD4","CD8A","CD8B","CD6")
#KIRs <- c("IKZF2", "KLRC3", "KIR2DL3")
Act <- c( "NKG7",  "GZMH", "GZMB", "HLA-DRB1", "CCR7","CCR5","CCR6", "SELL")
Mem <- c("CD27", "CD69", "KLRG1", "KLRD1", "KLRB1","CX3CR1", "ITGAE", "IL7R","IL2RA","IL2RB")
TF <- c( "TCF7","STAT3","STAT4","PRDM1","ID2","ID3","ZNF683","BCL6","EOMES","TBX21","RORC","FOXP3", "GATA3", "IKZF3", "TOX2")
cytokines <- c("IFNG","TNF","IL2","IL21","IL10","AREG","IL21R","IL23R","IL12RB2","IL18R1")

genesA <- c(basic, TF, cytokines, Act, Mem)

# DotPlot(Tcells, features =  genesA , dot.scale = 10) + 
#   scale_y_discrete(limits = levels(Tcells@active.ident)) +
#   coord_flip() + scale_color_viridis() + 
#   theme_minimal(base_size = 12) + 
#   theme(axis.text.x = element_text(angle = 90))
dp <- Clustered_DotPlot(Tcells, top_DEGs)

pdf(file = "./10x/outs/Set2_dotplot_markers.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 13)
dp[[2]]
dev.off()
png(file = "./10x/outs/Set2_dotplot_markers.png",   # The directory you want to save the file in
    width = 10,
    height = 13,
    units = 'in',
    res = 300)
dp[[2]]
dev.off()

# umap by T cell phenoytpe ====

TtypeUMAP <- DimPlot(Tcells, group.by = "Tcell_type") + 
  scale_color_jcolors("pal3") + 
  theme( axis.line = element_blank() ) + 
  labs(title = "T cell phenotype")

# What's moving with time? ----

TimeUMAP <- DimPlot(Tcells, group.by = "day", split.by = 'year') + 
  scale_color_manual(values = TimePal) + 
  theme( axis.line = element_blank() ) + 
  labs(title = "Days post-immunization")

TimeUMAPsplit <- DimPlot(Tcells, group.by = "day", split.by = "ident") + 
  scale_color_manual(values = TimePal) + 
  theme( axis.line = element_blank() )

TimeUMAPsplitDay <- DimPlot(Tcells, group.by = "ident", split.by = "day") + 
  theme( axis.line = element_blank() )

TimeUMAP2 <- DimPlot(Tcells, group.by = "time_point") + 
  scale_color_manual(values = TimePal2) + 
  theme( axis.line = element_blank() ) + 
  labs(title = "Time point")

YearUMAP <- DimPlot(Tcells, group.by = "year") + 
  scale_color_manual(values = YearPal) + 
  theme( axis.line = element_blank() ) + 
  labs(title = "Study year")

# What's changing with the location? ----

# mark LN vs blood samples
TissueUMAP <- DimPlot(Tcells, group.by = "tissue") + 
  scale_color_manual(values = TissuePal) + 
  theme( axis.line = element_blank() ) + 
  labs(title = "Tissue source", subtitle = "PBMC = blood \nFNA = lymph node")
TissueUMAPsplitID <- DimPlot(Tcells, group.by = "tissue", split.by = "ident") + 
  scale_color_manual(values = TissuePal) + 
  theme( axis.line = element_blank() )
TissueUMAPsplitTissue <- DimPlot(Tcells, group.by = "ident", split.by = "tissue") + 
  scale_color_manual(values = TissuePal) + 
  theme( axis.line = element_blank() )

#What's changing with the Donor? ----
# mark LN vs blood samples
DonorUMAP <- DimPlot(Tcells, group.by = "donor") + 
  scale_color_manual(values = DonorPal) +
  theme( axis.line = element_blank() ) + 
  labs(title = "Donor source")
DonorUMAPsplitID <- DimPlot(Tcells, group.by = "donor", split.by = "ident") + 
  scale_color_manual(values = DonorPal) +
  theme( axis.line = element_blank() )
DonorUMAPsplitDonor <- DimPlot(Tcells, group.by = "ident", split.by = "donor") + 
  theme( axis.line = element_blank() )

# The original tissue tends to drive the differences between clusters. You can see nice changes over time for a "super cluster "
TimeTissueDonor <- FetchData(Tcells, vars = c("ident", "day", "tissue", "donor" ,"Tcell_type", 'year', 'time_point') )
TimeTissueDonor <- setDT(TimeTissueDonor, keep.rownames = T) %>% as.data.frame()
colnames(TimeTissueDonor)[1] <- "cell"

DayGrp <- TimeTissueDonor %>% 
  group_by(time_point, ident) %>% 
  tally() %>% 
  mutate(freq = n / sum(n))
TissueGrp <- TimeTissueDonor %>% 
  group_by(ident, tissue) %>% 
  tally() %>% 
  mutate(freq = n / sum(n)) %>% 
  mutate_if(is.character, as.factor)

DonorGrp <- TimeTissueDonor %>% 
  group_by(ident, donor) %>% 
  tally() %>% 
  mutate(freq = n / sum(n)) %>% 
  mutate_if(is.character, as.factor)

YearGrp <- TimeTissueDonor %>% 
  group_by(year, ident) %>% 
  tally() %>% 
  mutate(freq = n / sum(n)) %>% 
  mutate_if(is.character, as.factor)

day_grp_plt <- ggplot(DayGrp, aes(ident, freq, fill =time_point)) + 
  geom_col(position = position_fill(reverse = TRUE)) + 
  scale_fill_manual(values = TimePal2) + 
  theme_minimal() +
  labs(title = "Freq of cells at each time per cluster") + 
  theme(axis.text.x = element_text(angle = 90))

tissue_grp_plt <- ggplot(TissueGrp, aes(ident, freq, fill = tissue)) + 
  geom_col(position = position_fill(reverse = TRUE)) + 
  scale_fill_manual(values = TissuePal) + 
  theme_minimal() +
  labs(title = "Freq of cells from tissue per cluster") + 
  theme(axis.text.x = element_text(angle = 90))

donor_grp_plt <- ggplot(DonorGrp, aes(ident, freq, fill = donor)) + 
  geom_col(position = position_fill(reverse = TRUE)) + 
  scale_fill_manual(values = DonorPal) + 
  theme_minimal() +
  labs(title = "Freq of cells from donor per cluster") + 
  theme(axis.text.x = element_text(angle = 90))

year_grp_plt <- ggplot(YearGrp, aes(ident, freq, fill = year)) + 
  geom_col(position = position_fill(reverse = TRUE)) + 
  scale_fill_manual(values = YearPal) + 
  theme_minimal() +
  labs(title = "Freq of cells from by year per cluster") + 
  theme(axis.text.x = element_text(angle = 90))

dist_big <- ( TissueUMAP| DonorUMAP | YearUMAP |  TimeUMAP2 ) / ( tissue_grp_plt | donor_grp_plt | year_grp_plt | day_grp_plt )  
ggsave('./10x/outs/Set2_Tcell_factor_cluster_distributions_w_UMAPs.png', plot= dist_big, height = 12, width = 16)

dist <- ( tissue_grp_plt | donor_grp_plt) / ( year_grp_plt | day_grp_plt) 
ggsave('./10x/outs/Set2_Tcell_factor_cluster_distributions.png', plot= dist, height = 10, width = 8)
dist2 <- ( tissue_grp_plt / donor_grp_plt / year_grp_plt / day_grp_plt) 
ggsave('./10x/outs/Set2_Tcell_factor_cluster_distributions.pdf', plot= dist2, height = 16, width = 8)

# Where are the clonal expansions ====

CloneUMAP <- FeaturePlot(Tcells, features = 'clone_size_log' ) + 
  theme( axis.line = element_blank(),
         plot.title = element_text(hjust = 0),
         plot.subtitle = element_text(hjust = 0)) + 
  labs(title = 'Clone size', subtitle = 'log( clone size + 1)') +
  scale_color_viridis(option = 'C')


#Make pretty umaps  ====
cpal <- cols25(length(levels(Tcells@meta.data$ident)))
names(cpal) <- levels(Tcells@meta.data$ident)

UMAP <- DimPlot(Tcells, group.by = "ident", label =T) + 
  theme( axis.line = element_blank() ) + 
  labs(title = "Cluster identity") +
  scale_color_manual(values = cpal)

UMAPdonor <- DimPlot(Tcells, group.by = "donor") + 
  scale_color_manual(values = DonorPal) + 
  theme( axis.line = element_blank() ) + 
  labs(title = "Donor") 

uplots <- ( UMAP | UMAPdonor | TtypeUMAP ) / ( TissueUMAP| YearUMAP | TimeUMAP2 )

ggsave('./10x/outs/Set2_Tcell_umaps.pdf', plot = uplots, height = 12, width = 21, useDingbats =F)
ggsave('./10x/outs/Set2_Tcell_umaps.png', plot = uplots, height = 12, width = 21)


# Pretty feature plots ====

Feats <- c('CD8B', 'CD4', 'CCR7', 'SELL', "KLRG1", 'KLRB1', 
           "FOXP3", 'ICOS','TOX2','CXCR5', "PDCD1", 'IL21')

FPs_panel <- FeaturePlot_scCustom(Tcells, Feats, num_columns = 6) & 
  NoLegend() & 
  xlab('UMAP 1') &
  xlab('UMAP 2') &
  theme(plot.title = element_text(face = 'italic', hjust = 0.5, size = 7) , 
        axis.text = element_text(size = 7), 
        axis.title = element_text(size = 7)) 

ggsave('./10x/outs/Set2_Tcell_feature_plots.png', plot = FPs_panel, height = 5, width = 12)
ggsave('./10x/outs/Set2_Tcell_feature_plots.pdf', plot = FPs_panel, height = 5, width = 12)

# one page plot ====
adj <- theme(axis.text = element_text(size = 10) , axis.title = element_text(size = 12))

big_plot_top <-  ( UMAP | UMAPdonor| TtypeUMAP)  / ( TissueUMAP | YearUMAP | TimeUMAP2 ) & adj
                  
big_plot <- big_plot_top + FPs_panel + plot_layout( heights = c(1,1,2)) & theme(text = element_text(size = 16) ) 

ggsave('./10x/outs/Set2_Tcell_figure.png', plot = big_plot, height = 14, width = 16)

#save updated T cell object ----

saveRDS(Tcells , Set2_integrated_Tcells_path )


