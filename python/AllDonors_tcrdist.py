# -*- coding: utf-8 -*-
"""
Created on Mon May 24 09:29:30 2021

@author: sschattg
"""

import conga
import numpy as np
import pandas as pd
from conga.tcrdist.tcr_distances import TcrDistCalculator
import seaborn as sns

from sklearn.manifold import TSNE
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial import distance
import matplotlib.pyplot as plt
import os 
os.chdir('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/')


# import scPCR clones from all three donors
sc_pcr_df = pd.read_csv('./repdata/sc/clone_tables/Ali_paired_sc_clones_JCC283.tsv', sep = '\t')
sc_pcr_df['match.cdr'] = sc_pcr_df['cdr3a'] + ";" + sc_pcr_df['cdr3b']
sc_pcr_df['epitope'] = sc_pcr_df['epitope'].str.replace('_CTF', '')
sc_pcr_df['epitope'] = sc_pcr_df['epitope'].str.replace('_', '-')
sc_pcr_df['clone_id'] = sc_pcr_df['clone_id'].str.replace('.clone', '') + "_" + sc_pcr_df['epitope'] 
sc_pcr_df = sc_pcr_df.rename(  columns =  {'epitope': 'Donor'})

for i in ['va_gene','ja_gene','vb_gene','jb_gene']:
    sc_pcr_df[i] = sc_pcr_df[i].str.replace('[*]..', '', regex= True)
    
# import the 10x clones file and subset on the Tfh only clones
x_clone = pd.read_csv('./10x/outs/TwoYear_10x_clones.tsv', sep = '\t')
tfh_x_clone = x_clone[ x_clone['Tcell_type'] == 'Tfh' ].copy()
tfh_cols = tfh_x_clone.columns.tolist()

# reformat the scpcr to merge with tfh_x_clone
sc_pcr_df2 = sc_pcr_df[sc_pcr_df.columns.intersection(tfh_cols)]
sc_pcr_df2 = sc_pcr_df2[ ~ sc_pcr_df2['match.cdr'].isin(tfh_x_clone['match.cdr']) ]
sc_pcr_df2 = sc_pcr_df2.drop_duplicates( subset=['match.cdr'])

# match the cols
new_cols = sc_pcr_df2.columns.tolist()
tfh_x_clone = tfh_x_clone[tfh_x_clone.columns.intersection(new_cols)]
tfh_x_clone = tfh_x_clone.drop_duplicates( subset=['match.cdr'])

# merge df
all_clones_df = pd.concat([tfh_x_clone, sc_pcr_df2])
all_clones_df = all_clones_df.reset_index(drop = True)
all_clones_df = all_clones_df.rename( columns = {'va_gene':'va',
                                                 'ja_gene':'ja',
                                                 'vb_gene':'vb',
                                                 'jb_gene':'jb'})
all_clones_df['donor'] = all_clones_df['clone_id'].str.split('_', expand=True)[[2]]

# adding alleles
for i in ['va','ja','vb','jb']:
    all_clones_df[i] = all_clones_df[i] + '*01'
    
# picked clones seqs 
pickedClones = pd.read_csv('./10x/outs/TwoYear_Donor_321-05_pickedClones_annotated.csv')
match_seq_df = pickedClones[ (pickedClones.confirmed_flu_specific == 'yes') & 
             (pickedClones.clone_id != 'Tfh_321-05_clone_6')]    
match_seq = match_seq_df.cdr3a + ";" + match_seq_df.cdr3b

# note the picked clones
all_clones_df['picked'] = 'no'
all_clones_df.loc[np.where( all_clones_df['match.cdr'].isin(match_seq))[0].tolist(), 
                  'picked'] = 'yes'

#save input df
all_clones_df.to_csv('./10x/outs/AllDonors_Tfh_clones_tcrdist.csv', index = False )

#run TCRdist and save
tcrs = conga.tcr_clumping.tcrs_from_dataframe_helper(all_clones_df, add_j_and_nucseq=True )
tdist = TcrDistCalculator( 'human' )
tcrdist_mat = np.array([tdist(x,y) for x in tcrs for y in tcrs]).reshape((len(tcrs), len(tcrs)))
np.savetxt('./10x/outs/AllDonors_Tfh_tcrdist_matrix.csv', tcrdist_mat, delimiter=",")

#kPC reduction and save
from sklearn.decomposition import KernelPCA
transformer = KernelPCA(n_components=50, kernel='linear')
tcrdist_mat_kpc = transformer.fit_transform(tcrdist_mat)
kpc_df = pd.DataFrame(tcrdist_mat_kpc, index = all_clones_df.clone_id)
kpc_df.to_csv('./10x/outs/AllDonors_Tfh_tcrdist_kpc.csv')





"""
from sklearn.cluster import KMeans
kmeans = KMeans(n_clusters=15, random_state=0).fit(tcrdist_mat_kpc)
assert len(kmeans.labels_) == all_clones_df.shape[0]
all_clones_df['kmeans'] = kmeans.labels_

X_embedded = TSNE(n_components=2).fit_transform(tcrdist_mat)



### prep for plotting
tsne_df1 = pd.DataFrame(X_embedded, columns = ["tSNE1", "tSNE2"], index = all_clones_df.index )

# add kmeans cluster and tSNE coordinates to df
all_clones_df = all_clones_df.join(tsne_df1)


# plot results
#heatmap of raw pw tcrdist
tcrdist_df = pd.DataFrame(tcrdist_mat, 
                          columns = all_clones_df.donor , 
                          index = all_clones_df.donor)

lut = dict(zip(tcrdist_df.index.unique(), ["#29BF12", "#00A5CF", "#DE1A1A"]))
row_colors = tcrdist_df.index.map(lut)
col_colors = tcrdist_df.columns.map(lut)
tcrdist_hm = sns.clustermap(tcrdist_df, 
                            row_colors = row_colors, 
                            col_colors = col_colors )

plt.savefig('./10x/outs/AllDonors_Tfh_tcrdist_hm.png', dpi=400)
"""
