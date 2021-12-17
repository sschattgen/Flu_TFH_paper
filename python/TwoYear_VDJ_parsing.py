# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import glob
import re
import pandas as pd
import numpy as np
os.chdir("/Users/sschattg/conga")
import conga
from conga.tcrdist.make_10x_clones_file import make_10x_clones_file_batch 
os.chdir("/Users/sschattg/VDJoutputs")

vdj_files = glob.glob("./*_filtered_contig_annotations.csv")

for i in vdj_files:
    
    df = pd.read_csv(i)
    
    bcr = df.loc[(df['chain'] == "IGH" ) | (df['chain'] ==  "IGL") | (df['chain'] ==  "IGK")]
    tcr = df.loc[(df['chain'] == "TRA" ) | (df['chain'] == "TRB")]
    
    if bcr.shape[0] != 0:
        bcr.to_csv("./bcr_only/" + i , index = False )
    if tcr.shape[0] != 0:
        tcr.to_csv("./tcr_only/" + i , index = False )

vdj_meta = pd.read_csv('vdj_metadata.csv')

tcr_files = glob.glob("./tcr_only/*_filtered_contig_annotations.csv")
tcr_metadata = vdj_meta.copy()
tcr_metadata['file'] = './tcr_only/' + tcr_metadata['file'] 
tcr_metadata = tcr_metadata[ tcr_metadata['file'].isin(tcr_files) ]
tcr_metadata.to_csv('tcr_metadata.csv', index = False)

bcr_files = glob.glob("./bcr_only/*_filtered_contig_annotations.csv")
bcr_metadata = vdj_meta.copy()
bcr_metadata['file'] = './bcr_only/' + bcr_metadata['file'] 
bcr_metadata = bcr_metadata[ bcr_metadata['file'].isin(bcr_files) ]
bcr_metadata.to_csv('bcr_metadata.csv', index = False)

make_10x_clones_file_batch(metadata_file = 'tcr_metadata.csv', 
                           organism = 'human',
                           clones_file = 'TwoYear_TCR_clones.tsv', 
                           stringent = True,
                           prefix_clone_ids_with_tcr_type = True)


make_10x_clones_file_batch(metadata_file = 'bcr_metadata.csv', 
                           organism = 'human_ig',
                           clones_file = 'TwoYear_BCR_clones.tsv', 
                           stringent = True,
                           prefix_clone_ids_with_tcr_type = True)


