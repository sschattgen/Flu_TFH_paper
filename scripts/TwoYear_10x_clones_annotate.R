# this script generates table of clonotype data for all cells in the BothDonors_Tcells dataset 
# match to sc and bulk rep data and associated plots

source('scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

# import T cell seurat object, and the scPCR and bulk data ----

Tcells <- readRDS(Tcells_path)

# paired scPCR and bulk data from sorted Tfh cells
pairedSC <- read.delim('data/Ali_paired_sc_clones_JCC283.tsv', 
                       stringsAsFactors = F)

bulk_clones.vdj <- read.delim("data/Ali_bulk_allDonor_clones.tsv", 
                              stringsAsFactors = F)

# generate a clone_df from the Tcells metadata ====

Clonedf <- makeClonedf(Tcells, 
                       stripAlleles = TRUE, 
                       other_metadata =  c( "donor", "time", "year", 
                                            "day", "tissue", "Tcell_type", 
                                            "strict.match.cdr")
                       )

# adding a less strict concatenated cdr3 sequences for matching to the scPCR data
Clonedf <- Clonedf %>%
  mutate(match.cdr = paste(cdr3a, cdr3b, sep = ";"))
# count clones within cluster and in total. 
Clonedf <- countClones(Clonedf)

# are any clones shared between the donors?
donorShared_clones <- Clonedf %>% 
  select(strict.match.cdr, donor) %>%
  group_by_all()%>% 
  distinct(.keep_all = T) %>%
  ungroup(donor) %>%
  add_count(strict.match.cdr, name = "num_libraries") %>% 
  select(-donor) %>%
  distinct(.keep_all = T) %>%
  filter(num_libraries >= 2 )

# one clone matched across the donor, and it's one cell each, but it is a Tfh!
multidonor_clone <- Clonedf %>% 
  filter(strict.match.cdr %in% donorShared_clones$strict.match.cdr ) %>%
  pull(barcode)
Clonedf <- Clonedf %>%
  mutate(clone_in_both_donors = ifelse(barcode %in% multidonor_clone,
                                       'yes', 'no')
         )

#mark MAITs and NKTs
Clonedf <- detectUnconventionals(Clonedf, species = "Human")

# plotting of some basic info about the clones ====

# compare CDR3 lengths 
dClonesLen <- Clonedf %>% 
  select(ident, clone_id, cdr3a_length, cdr3b_length) %>% 
  group_by(ident, clone_id) %>% distinct(clone_id, .keep_all = T) %>% 
  ungroup(clone_id)

#whats the global mean?
uC <- Clonedf %>% 
  select(clone_id, cdr3a_length, cdr3b_length) %>%
  distinct(clone_id, .keep_all = T) 

alenm = median(uC$cdr3a_length)
blenm = median(uC$cdr3b_length) 

cdr3a_len <- ggplot(dClonesLen, aes(x = factor(ident), y =cdr3a_length, group = ident)) + 
  geom_boxplot() + 
  theme_minimal() +
  labs(title = "CDR3 length, alpha chain") + 
  geom_hline(yintercept = alenm) + 
  xlab("ident") +
  theme(axis.text.x = element_text(angle = 90)) 

cdr3b_len <- ggplot(dClonesLen, aes(x = factor(ident), y =cdr3b_length, group = ident)) + 
  geom_boxplot() + 
  theme_minimal() + 
  labs(title = "CDR3 length, beta chain") + 
  geom_hline(yintercept = blenm) + 
  xlab("ident") + 
  theme(axis.text.x = element_text(angle = 90)) 

# where are the expanded clones? 
countident <- Clonedf %>% 
  group_by(ident, clone_id) %>% 
  distinct(clone_id, .keep_all = T) %>% 
  select(clone_id, ident, count_cluster)

ctClusterUMAP <- ggplot(Clonedf , aes(UMAP_1, UMAP_2, color = log10(1 + count_total) )) + 
  geom_point(size = 0.5 ) +  scale_color_viridis(option = "C") + 
  theme_minimal() + 
  labs(title = "Clone size total dataset")

ctClusterJitter <- ggplot(countident, aes(ident, log10(1 +count_cluster))) + geom_jitter(aes(color = log10(1 +count_cluster))) + 
  scale_color_viridis() + 
  theme_minimal() + 
  labs(title = "Clone size within cluster")  + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) 

csize_plot <- (ctClusterUMAP / ctClusterJitter) + plot_layout(heights = c(2,1))
ggsave('outs/TwoYear_Tcell_cloneSizes.png', plot = csize_plot, height = 10, width = 7)


# plot the top100clones and their cluster freq
TopClones <- countTopClones(Clonedf , clone_counts = c(100))
tC <- TopClones[[1]]$clone_id

TopIdent <- Clonedf %>% 
  filter(clone_id %in% tC) %>% 
  mutate(clone_id = factor(clone_id, levels = tC)) %>%
  group_by(clone_id, ident) %>% 
  select(clone_id, ident) %>% 
  tally() %>% 
  mutate(freq = n / sum(n), clone_size = sum(n)) 

TopSplitCluster <- ggplot(TopIdent , aes(clone_id, freq, fill = ident)) +  
  geom_col(position = "fill") +
  theme_minimal() + 
  labs(title = "Cluster distribution of Top100 clones") + 
  geom_text(data=TopClones[[1]], aes(clone_id, y =1.1, label= n), inherit.aes = F, angle=90,hjust=1, vjust= 0.5) +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom" ) 

TopSplitCluster
ggsave('outs/TwoYear_TopClones_clusters.png', height = 5, width = 20)


# match to paired single cell sorted data ----

#split by donor
d4_sc <- pairedSC %>% filter( subject %like% "321_04")
d5_sc <- pairedSC %>% filter( subject %like% "321_05")

#match pairs
Clonedf$donor_321_04_paired_scPCR_match <- ifelse(Clonedf$match.cdr %in% d4_sc$match.cdr, "yes", "no" )
Clonedf$donor_321_05_paired_scPCR_match <- ifelse(Clonedf$match.cdr %in% d5_sc$match.cdr, "yes", "no" )
Clonedf$donor_321_05_paired_scPCR_match <- ifelse(Clonedf$match.cdr %in% d5_sc$match.cdr, "yes", "no" )
#match single chains
Clonedf$donor_321_04_alpha_scPCR_match <- ifelse(Clonedf$match.cdr %in% d4_sc$cdr3a, "yes", "no" )
Clonedf$donor_321_05_alpha_scPCR_match <- ifelse(Clonedf$match.cdr %in% d5_sc$cdr3a, "yes", "no" )
Clonedf$donor_321_04_beta_scPCR_match <- ifelse(Clonedf$match.cdr %in% d4_sc$cdr3b, "yes", "no" )
Clonedf$donor_321_05_beta_scPCR_match <- ifelse(Clonedf$match.cdr %in% d5_sc$cdr3b, "yes", "no" )

# plotting

cols <- c("grey90", DonorPal[1])
names(cols) <- c("no", "yes")
cols2 <- c("grey90", DonorPal[2])
names(cols2) <- c("no", "yes")

d4_match_scPCR_UMAP <- ggplot(Clonedf) + 
  geom_point(data = subset(Clonedf, donor_321_04_paired_scPCR_match == 'no'),
             aes(UMAP_1, UMAP_2, colour = donor_321_04_paired_scPCR_match)) +
  geom_point(data = subset(Clonedf, donor_321_04_paired_scPCR_match != 'no'),
             aes(UMAP_1, UMAP_2, colour = donor_321_04_paired_scPCR_match)) +
  scale_colour_manual(values = cols , breaks = cols ) +
  theme_minimal() +
  labs(title = "Clone match to paired scPCR from Donor 321-04")

d5_match_scPCR_UMAP <- ggplot(Clonedf) + 
  geom_point(data = subset(Clonedf, donor_321_05_paired_scPCR_match == 'no'),
             aes(UMAP_1, UMAP_2, colour = donor_321_05_paired_scPCR_match)) +
  geom_point(data = subset(Clonedf, donor_321_05_paired_scPCR_match != 'no'),
             aes(UMAP_1, UMAP_2, colour = donor_321_05_paired_scPCR_match)) +
  scale_colour_manual(values = cols2 , breaks = cols2 ) +
  theme_minimal() +
  labs(title = "Clone match to paired scPCR from Donor 321-05")

d4_match_scPCR_UMAP | d5_match_scPCR_UMAP
ggsave('outs/TwoYear_scPCR_pairedMatch.png', height = 6.5, width = 12)

# bulk TCRs -----
#match the bulk back to seurat object

d5_bulk <- bulk_clones.vdj %>% filter(donor == "321-05" )
d4_bulk <- bulk_clones.vdj %>% filter(donor == "321-04" )

Clonedf$matched_321_04_alpha_bulk <- ifelse(Clonedf$cdr3a %in% d4_bulk$cdr3aa, "yes", "no" )
Clonedf$matched_321_05_alpha_bulk <- ifelse(Clonedf$cdr3a %in% d5_bulk$cdr3aa, "yes", "no" )
Clonedf$matched_321_04_beta_bulk <- ifelse(Clonedf$cdr3b %in% d4_bulk$cdr3aa, "yes", "no" )
Clonedf$matched_321_05_beta_bulk <- ifelse(Clonedf$cdr3b %in% d5_bulk$cdr3aa, "yes", "no" )


d4_match_alpha_UMAP <- ggplot(Clonedf) + 
  geom_point(data = subset(Clonedf, matched_321_04_alpha_bulk  == 'no'),
             aes(UMAP_1, UMAP_2, colour = matched_321_04_alpha_bulk )) +
  geom_point(data = subset(Clonedf, matched_321_04_alpha_bulk  != 'no'),
             aes(UMAP_1, UMAP_2, colour = matched_321_04_alpha_bulk )) +
  scale_colour_manual(values = cols , breaks = cols ) +
  theme_minimal() +
  labs(title = "CDR3a matched bulk from Donor 321-04")

d4_match_beta_UMAP <- ggplot(Clonedf) + 
  geom_point(data = subset(Clonedf, matched_321_04_beta_bulk  == 'no'),
             aes(UMAP_1, UMAP_2, colour = matched_321_04_beta_bulk )) +
  geom_point(data = subset(Clonedf, matched_321_04_beta_bulk  != 'no'),
             aes(UMAP_1, UMAP_2, colour = matched_321_04_beta_bulk )) +
  scale_colour_manual(values = cols , breaks = cols ) +
  theme_minimal() +
  labs(title = "CDR3b matched bulk from Donor 321-04")

d5_match_alpha_UMAP <- ggplot(Clonedf) + 
  geom_point(data = subset(Clonedf, matched_321_05_alpha_bulk  == 'no'),
             aes(UMAP_1, UMAP_2, colour = matched_321_05_alpha_bulk )) +
  geom_point(data = subset(Clonedf, matched_321_05_alpha_bulk  != 'no'),
             aes(UMAP_1, UMAP_2, colour = matched_321_05_alpha_bulk )) +
  scale_colour_manual(values = cols2 , breaks = cols2 ) +
  theme_minimal() +
  labs(title = "CDR3a matched bulk from Donor 321-05")

d5_match_beta_UMAP <- ggplot(Clonedf) + 
  geom_point(data = subset(Clonedf, matched_321_05_beta_bulk  == 'no'),
             aes(UMAP_1, UMAP_2, colour = matched_321_05_beta_bulk )) +
  geom_point(data = subset(Clonedf, matched_321_05_beta_bulk  != 'no'),
             aes(UMAP_1, UMAP_2, colour = matched_321_05_beta_bulk )) +
  scale_colour_manual(values = cols2 , breaks = cols2 ) +
  theme_minimal() +
  labs(title = "CDR3b matched bulk from Donor 321-05")

d4_match_alpha_col <- ggplot( Clonedf ,aes(x = ident, y = freq_cluster, fill = matched_321_04_alpha_bulk ) ) + 
  geom_col(position = "fill") + 
  scale_fill_manual(values = cols , breaks = cols ) + 
  theme_minimal() +
  labs(title = "CDR3a matched bulk from Donor 321-04") + 
  theme(axis.text.x = element_text(angle = 90)) 

d4_match_beta_col <- ggplot( Clonedf ,aes(x = ident, y = freq_cluster, fill = matched_321_04_beta_bulk ) ) + 
  geom_col(position = "fill") + 
  scale_fill_manual(values = cols , breaks = cols ) + 
  theme_minimal() +
  labs(title = "CDR3b matched bulk from Donor 321-04") + 
  theme(axis.text.x = element_text(angle = 90)) 

d5_match_alpha_col <- ggplot( Clonedf ,aes(x = ident, y = freq_cluster, fill = matched_321_05_alpha_bulk ) ) + 
  geom_col(position = "fill") + 
  scale_fill_manual(values = cols2 , breaks = cols2 ) + 
  theme_minimal() +
  labs(title = "CDR3a matched bulk from Donor 321-05") + 
  theme(axis.text.x = element_text(angle = 90)) 

d5_match_beta_col <- ggplot( Clonedf ,aes(x = ident, y = freq_cluster, fill = matched_321_05_beta_bulk ) ) + 
  geom_col(position = "fill") + 
  scale_fill_manual(values = cols2 , breaks = cols2 ) + 
  theme_minimal() +
  labs(title = "CDR3b matched bulk from Donor 321-05") + 
  theme(axis.text.x = element_text(angle = 90)) 

(d4_match_alpha_UMAP | d4_match_beta_UMAP | d5_match_alpha_UMAP | d5_match_beta_UMAP) / 
  (d4_match_alpha_col | d4_match_beta_col | d5_match_alpha_col | d5_match_beta_col)

ggsave('outs/TwoYear_bulk_matched_chains.png', height = 10, width = 22)


# add the picked clones from the original analysis ====

# the one year picked clones sheet
pickedClones_yr1 <- read.csv('data/Donor_321-05_pickedClones_yr1.csv', stringsAsFactors = F)

pickedClones_yr1_short <- pickedClones_yr1 %>%
  mutate(strict.match.cdr = paste(va_gene, cdr3a, ja_gene,
                                  vb_gene, cdr3b, jb_gene,
                                  sep = ";"),
         Tfh_clone_id = clone_id) %>%
  select(Tfh_clone_id, strict.match.cdr)

pickedClones_yr1_short$strict.match.cdr %in% Clonedf$strict.match.cdr

Clonedf <- left_join(Clonedf, pickedClones_yr1_short)


# add clone_size to Tcells
counts <- Clonedf$count_total
names(counts) <- Clonedf$barcode
Tcells <- AddMetaData( Tcells, metadata = counts, col.name = 'clone_size')

counts_log <- log(Clonedf$count_total + 1)
names(counts_log) <- Clonedf$barcode
Tcells <- AddMetaData( Tcells, metadata = counts_log, col.name = 'clone_size_log')

# save the Clonedf and Tcells ====

write_tsv( Clonedf , clone_df_path )
saveRDS(Tcells , Tcells_path )
