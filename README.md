## Data and scripts for study by Schattgen, Turner, Ghonim et al. 2023 bioRxiv 2023.08.29.555186; doi: https://doi.org/10.1101/2023.08.29.555186

This R code will run the entire workflow beginning after preprocessing and subsetting on only T cells. The B cell data is being used for a different publication and needs to be withheld for now. `TwoYear_VDJ_parsing.py` and `TwoYear_gen_seurat.R` show the code used for preprocessing prior to subsetting on T cells. 

## How to run full workflow (assumes you have R installed)
1. Clone this repository
2. Download the the dataset from Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6476022.svg)](https://doi.org/10.5281/zenodo.6476022)
3. Unpack with `tar -xzvf flu_vaxx_v2_dataset.tar.gz` and move the `data` folder to the cloned repo from step 1.
4. Open `scripts/tfh_pkg_paths_vars.R` and replace `path/to/repository` on line 2 with the path to the repo and data 

### Script Order

```
# These scripts were used for preprocessing and included to show the steps taken

TwoYear_gen_seurat.R
d05_add_PBMC_gen_seurat.R
Set2_gen_seurat.R
harmony_integration.R
integrate_clonotypes.R

# Starting with intergrated_Tcells_harmony_bydonor.rds the following scripts will run the analyses

setup.R
Set2_10x_clones_annotate.R
Set2_Tcell_annotation.R
Set2_TCell_analysis.R
Set2_Tcell_metabolism.R
Set2_CD4_subset.R
Set2_Tfh_subset.R
Set2_Tfh_phenotypes.R
Set2_Tfh_time.R
Set2_Tfh_clonotypes.R
Set2_Donor_321-05_pickedClones.R
Set2_Donor_321-05_pickedClones_gex_variation.R
Set2_Donor_321-05_pickedClones_modules.R
Set2_monocle.R
Set2_Tfh_lineage_network_graph.R
```