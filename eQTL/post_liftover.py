# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 14:39:49 2019

@author: bingkun

@project hfb -- H-MAGMA/eQTL

Script for processing post-liftover interations
S2 of pipeline
* collect boths ends of interactions after liftover
* data cleaning
* output file for intersection
"""

import pandas as pd

cellType='RGC'

#%%
################################## PART ONE #####################################

# read in interactions before liftover, with ID
raw_all = pd.read_csv(r'C:\Users\libin\UCSF\MAGMA\{}.5k.downsample.bedpe.withID'.format(cellType), sep="\t")
print ("number of interactions: {}".format(raw_all.shape[0]))
raw_all_left = raw_all[['ID', 'chr1', 'start1', 'end1']]
raw_all_right = raw_all[['ID', 'chr2', 'start2', 'end2']]

# read in both ends of interactions after lifting over
end1 = pd.read_csv(r'C:\Users\libin\UCSF\MAGMA\{}.5k.downsample.bedpe.hg19.end1'.format(cellType), sep="\t")
end2 = pd.read_csv(r'C:\Users\libin\UCSF\MAGMA\{}.5k.downsample.bedpe.hg19.end2'.format(cellType), sep="\t")
end_both = pd.merge(end1, end2, how="inner", on='ID')
end_both = end_both[['chr1', 'start1', 'chr2', 'start2', 'ID']]
# remove interactions whose ends got mapped to diff chr after liftover
misplaced_index = end_both[end_both["chr1"] != end_both["chr2"]].index
end_both = end_both.drop(misplaced_index)

interactions_all_liftover = pd.merge(end_both, raw_all, how="inner", on=["chr1", "chr2", "ID"])
# compute distance of interactions before/after liftover
interactions_all_liftover["distance_new"] = abs(interactions_all_liftover["start2_x"] - interactions_all_liftover["start1_x"])
interactions_all_liftover["distance_old"] = abs(interactions_all_liftover["start2_y"] - interactions_all_liftover["start1_y"])
interactions_all_liftover["distance_diff"] = abs(interactions_all_liftover["distance_old"] - interactions_all_liftover["distance_new"])

distance_changed = interactions_all_liftover[interactions_all_liftover["distance_diff"] != 0]
interactions_all_liftover["end1_new"] = interactions_all_liftover["start1_x"] + 5000
interactions_all_liftover["end2_new"] = interactions_all_liftover["start2_x"] + 5000

interactions_all_liftover = interactions_all_liftover[[
 'chr1',
 'start1_x',
 'end1_new',
 'chr2',
 'start2_x',
 'end2_new',
 'count',
 'expected',
 'fdr',
 'ClusterLabel',
 'ClusterSize',
 'ClusterType',
 'ClusterNegLog10P',
 'ClusterSummit',
 'ID',
 'distance_new',
 'distance_old',
 'distance_diff']].rename(columns={"start1_x": "start1", "start2_x": "start2", "end1_new": "end1", "end2_new": "end2"})
interactions_all_liftover["end1"] = interactions_all_liftover["start1"] + 5000
interactions_all_liftover["end2"] = interactions_all_liftover["start2"] + 5000
print("number of interactons after liftover: {}".format(interactions_all_liftover.shape[0]))
#interactions_all_liftover_lh = interactions_all_liftover[['chr1', 'start1', 'end1', 'ID']]
#interactions_all_liftover_rh = interactions_all_liftover[['chr2', 'start2', 'end2', 'ID']]
#%%

############################ PART TWO -- LABEL PROMOTERS ############################
anchor_bins = pd.read_csv(r'C:\Users\libin\UCSF\MAGMA\MACS2_peak_5KbBin_Union4.txt', sep="\t")
anchor_bins["anchor_ID"] = ["anchor_{}".format(i) for i in range(anchor_bins.shape[0])]

raw_all_left_anchor_label = pd.merge(raw_all_left, anchor_bins, left_on=['chr1', 'start1', 'end1'],
                                     right_on = ['chr', 'start', 'end'], how="inner")
left_anchors_list = raw_all_left_anchor_label[['ID', 'anchor_ID']]

raw_all_right_anchor_label = pd.merge(raw_all_right, anchor_bins, left_on=['chr2', 'start2', 'end2'],
                                     right_on = ['chr', 'start', 'end'], how="inner")
right_anchors_list = raw_all_right_anchor_label[['ID', 'anchor_ID']]


anchor_anchor_list = pd.merge(left_anchors_list, right_anchors_list, on="ID", how="inner")
#### remove a-a interactions first
interactions_all_liftover_XOR = interactions_all_liftover[~interactions_all_liftover["ID"].isin(anchor_anchor_list["ID"])]
print ("number of XOR interactions after liftover: {}".format(interactions_all_liftover_XOR.shape[0]))

interactions_all_liftover_XOR_left_anchor = pd.merge(interactions_all_liftover_XOR, left_anchors_list, on="ID", how="inner")
print ("number of XOR interactions with bin1 as anchor: {}".format(interactions_all_liftover_XOR_left_anchor.shape[0]))
interactions_all_liftover_XOR_right_anchor = pd.merge(interactions_all_liftover_XOR, right_anchors_list, on="ID", how="inner")
print ("number of XOR interactions with bin2 as anchor: {}".format(interactions_all_liftover_XOR_right_anchor.shape[0]))

left_anchor_interaction_anchors = \
interactions_all_liftover_XOR_left_anchor[["chr1", "start1", "end1", "ID", "anchor_ID"]] \
.rename(columns={"chr1":"chr_ac", "start1":"start_ac", "end1":"end_ac"})

left_anchor_interaction_targets = \
interactions_all_liftover_XOR_left_anchor[["chr2", "start2", "end2", "ID", "anchor_ID"]] \
.rename(columns={"chr2":"chr_tg", "start2":"start_tg", "end2":"end_tg"})

right_anchor_interaction_anchors = \
interactions_all_liftover_XOR_right_anchor[["chr2", "start2", "end2", "ID", "anchor_ID"]] \
.rename(columns={"chr2":"chr_ac", "start2":"start_ac", "end2":"end_ac"})

right_anchor_interaction_targets = \
interactions_all_liftover_XOR_right_anchor[["chr1", "start1", "end1", "ID", "anchor_ID"]] \
.rename(columns={"chr1":"chr_tg", "start1":"start_tg", "end1":"end_tg"})


interaction_all_anchors = pd.concat([left_anchor_interaction_anchors, right_anchor_interaction_anchors])
interaction_all_targets = pd.concat([left_anchor_interaction_targets, right_anchor_interaction_targets])

interaction_all_anchors.to_csv(r'C:\Users\libin\UCSF\MAGMA\{}_interaction_anchorBins.hg19.bed'.format(cellType), sep="\t", header=False, index=False)
interaction_all_targets.to_csv(r'C:\Users\libin\UCSF\MAGMA\{}_interaction_targetBins.hg19.bed'.format(cellType), sep="\t", header=False, index=False)






