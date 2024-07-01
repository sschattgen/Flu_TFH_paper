
# Picking out potentially flu-specific clonotypes for donor 321-05 for cloning
# NOTE: the clones were originally picked when we had data just for 321-05 in year 1. 
# Adding the second year made this hard to reproduce since the sizes of the orignal clones that were 
# originally picked  are different now. So this uses the pickedClone table from year 1 to guide the setup for both years.

source('Z:/ResearchHome/ResearchHomeDirs/thomagrp/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

# imports ====
Tfh_obj <- readRDS(Set2_integrated_Tfh_path)

# read in the 10X clones file 
Clonedf <- read.delim(clone_df_path, stringsAsFactors = F)

# tally by day and tissue
Tfh_clones_tally <- Clonedf %>%
  filter( !is.na(Tfh_clone_id)) %>%
  group_by(time_point, tissue, Tfh_clone_id, clone_id) %>%
  tally() %>%
  ungroup() %>%
  group_by(clone_id) %>%
  mutate(Total = sum(n)) %>%
  pivot_wider(values_from = n, names_from = time_point, values_fill = 0) %>%
  arrange(desc(Total)) %>%
  select(2:3,1,4,5:ncol(.)) 

#add tcr info
Clonedf_meta <- Clonedf %>%
  filter( !is.na(Tfh_clone_id)) %>%
  select(Tfh_clone_id, 
         va_gene, ja_gene, cdr3a, cdr3a_nucseq, 
         vb_gene, jb_gene, cdr3b, cdr3b_nucseq, 
         donor_321_05_paired_scPCR_match ,matched_321_05_alpha_bulk , matched_321_05_beta_bulk ) %>%
  distinct_all()

Tfh_clones <- left_join(Tfh_clones_tally, Clonedf_meta)

# mark PB lineage clones
tdata <- data.frame(Tfh_clone_id = unique(Tfh_clones$Tfh_clone_id))
pb_lin <- c()
for (i in tdata$Tfh_clone){
  print(i)
  test <- Tfh_clones_tally %>%
    filter(Tfh_clone_id == i) 
  pb_lin_in <- ifelse(nrow(test) == 2, 'yes', 'no')
  pb_lin <- append(pb_lin, pb_lin_in)
  
}
tdata$PB_lineage <- pb_lin
Tfh_clones <- left_join(Tfh_clones, tdata) %>% 
  arrange(Tfh_clone_id)

write_csv( Tfh_clones, "./10x/outs/Set2_Donor_321-05_pickedClones.csv")


# make compatible for stitchr ====

pickedClone2 <- Tfh_clones %>% 
  mutate_if(is.factor, as.character )

stitchrdf <- data.frame( TCR_name = pickedClone2$Tfh_clone_id,
                         TRAV = pickedClone2$va_gene,
                         TRAJ = pickedClone2$ja_gene,
                         TRA_CDR3 = pickedClone2$cdr3a,
                         TRBV = pickedClone2$vb_gene,
                         TRBJ = pickedClone2$jb_gene,
                         TRB_CDR3 = pickedClone2$cdr3b,
                         TRAC = "",
                         TRBC = "",
                         TRA_leader = "",
                         TRB_leader = "",
                         Linker = "P2A"
)
stitchrdf <- stitchrdf %>% distinct_all()
write_tsv( stitchrdf, './10x/outs/Set2_Donor_321-05_pickedClones_stitchr.tsv')

# any hits in Herati or vdjdb ----

#this is the HA specific cells
HAtcr <- read.csv('./Wherry_dataset/Herati_HA_specific.csv', stringsAsFactors = F)
HAtcr_chains <- HAtcr %>% 
  filter(SequenceStatus == 'In') %>% # remove out of frame or stops and focus on Tfh
  pull(CDR3.AminoAcid)

# import vdjdb-2020-05-20 
vdjdb <- read.delim('Z:/ResearchHome/Groups/thomagrp/home/sschattg/TOOLBOX/db/vdjdb-2020-05-20/vdjdb.txt', stringsAsFactors  = F)
vdjdb_II <- vdjdb %>% filter( mhc.class == 'MHCII' & species == 'HomoSapiens') 
vdjdb_II_TRA <- vdjdb_II %>% filter(gene == "TRA")
vdjdb_II_TRB <- vdjdb_II %>% filter(gene == "TRB")

Tfh_clones$Herati_match_CDR3b <- ifelse(Tfh_clones$cdr3b %in% HAtcr_chains , 'yes', 'no')
Tfh_clones$vdjdb_MHCII_match_CDR3a <- ifelse(Tfh_clones$cdr3a %in% vdjdb_II_TRA$cdr3 , 'yes', 'no')
Tfh_clones$vdjdb_MHCII_match_CDR3b <- ifelse(Tfh_clones$cdr3b %in% vdjdb_II_TRB$cdr3 , 'yes', 'no')


#nice pickedclones table ====
ag_specific <- c('Tfh_321-05_clone_1', 
                 'Tfh_321-05_clone_3', 
                 'Tfh_321-05_clone_6', 
                 'Tfh_321-05_clone_11', 
                 'Tfh_321-05_clone_12')

Tfh_clones$confirmed_flu_specific <- ifelse(Tfh_clones$Tfh_clone_id %in% ag_specific, 'yes', 'no')


write_csv( Tfh_clones, "./10x/outs/Set2_Donor_321-05_pickedClones_annotated.csv")

# save changes to clonedf ====
write_tsv(Clonedf, clone_df_path)


