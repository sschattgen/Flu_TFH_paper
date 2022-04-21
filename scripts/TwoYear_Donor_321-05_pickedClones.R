
# Picking out potentially flu-specific clonotypes for donor 321-05 for cloning
# NOTE: the clones were originally picked when we had data just for 321-05 in year 1. 
# Adding the second year made this hard to reproduce since the sizes of the orignal clones that were 
# originally picked  are different now. So this uses the pickedClone table from year 1 to guide the setup for both years.

source('scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)


# imports ====
Tfh_obj <- readRDS(Tfh_lineages_path)

# read in the 10X clones file 
Clonedf <- read.delim(clone_df_path, stringsAsFactors = F)

# tally by day and tissue
Tfh_clones_tally <- Clonedf %>%
  filter( !is.na(Tfh_clone_id)) %>%
  group_by(time, tissue, Tfh_clone_id, clone_id) %>%
  tally() %>%
  ungroup() %>%
  group_by(clone_id) %>%
  mutate(Total = sum(n)) %>%
  pivot_wider(values_from = n, names_from = time, values_fill = 0) %>%
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

write_csv( Tfh_clones, "outs/TwoYear_Donor_321-05_pickedClones.csv")


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
write_tsv( stitchrdf, 'outs/TwoYear_Donor_321-05_pickedClones_stitchr.tsv')

