## Data and scripts for study by Schattgen, Turner, Ghonim et al. 2022

This R code will run the entire workflow beginning after preprocessing and subsetting on only T cells with TCR information. The B cell data is being used for a different publication and needs to be withheld for now. `TwoYear_VDJ_parsing.py` and `TwoYear_gen_seurat.R` show the code used for preprocessing prior to subsetting on T cells. 

## How to run full workflow (assumes you have R installed)
1. Clone this repository
2. Download the the dataset from Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6476022.svg)](https://doi.org/10.5281/zenodo.6476022)
3. Unpack with `tar -xzvf flu_vaxx_dataset.tar.gz` and move the `data` folder to the cloned repo from step 1.
4. Open `scripts/tfh_pkg_paths_vars.R` and replace `path/to/repository` on line 2 with the path to the repo and data 
5. Open command prompt, move to repo folder, run `bash flu_tfh_pipeline.sh`

### Alternatively, each script can be run interactively in order.

The workflow will take several hours to run. The output will be generated in the `outs` and `data` directories. 
