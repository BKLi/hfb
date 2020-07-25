# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 11:16:03 2020

@author: bingkun
@project: hfb -- eQTL

sample nonsig interactions from MAPS pospoisson file
only sample non-sig interactions with the same anchors as sig-interactions
including both AND and XOR interactions

* S2 of alternative alternative pipeline
"""

#%%
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['figure.dpi']= 200


####### read in complete MAPS interaction sets
cellType='neuron'
#sample_power
power = 1
iteration = 5
date = 20200212
#%%

inter_all = pd.read_csv(r'C:\Users\libin\UCSF\hfb\eqtl\hfb\MAPS_allPeaks\{}.all.MAPS2_pospoisson.reform.withID'.format(cellType), sep="\t",
                      names=['ID','bin1_mid','bin2_mid','count','X1D_peak_bin1','X1D_peak_bin2','effective_length1',
                             'gc1','mappability1','short_count1','effective_length2','gc2','mappability2',
                             'short_count2','dist','logl','loggc','logm','logdist','logShortCount','chr','expected',
                             'p_val','expected2','ratio2','p_val_reg2','p_bonferroni','fdr'])
inter_all = inter_all[['ID', 'chr', 'bin1_mid', 'bin2_mid', 'X1D_peak_bin1', 'X1D_peak_bin2', 'fdr', "dist"]]
inter_all["dist"] = inter_all["dist"]*5000
print ("number of all interactions: {}".format(inter_all.shape[0]))

#%%
inter_all_AND = inter_all[(inter_all["X1D_peak_bin1"] == 1) & (inter_all["X1D_peak_bin2"] == 1)]
inter_all_bin1_anchor = inter_all[(inter_all["X1D_peak_bin1"] == 1) & (inter_all["X1D_peak_bin2"] == 0)]
inter_all_bin2_anchor = inter_all[(inter_all["X1D_peak_bin1"] == 0) & (inter_all["X1D_peak_bin2"] == 1)]

# rename so bin1 is always anchor bin for XOR interactions
inter_all_bin2_anchor = inter_all_bin2_anchor.rename(columns={"bin1_mid":"bin2_mid", "bin2_mid":"bin1_mid"})
inter_all_bin2_anchor = inter_all_bin2_anchor[['ID', 'chr', 'bin1_mid', 'bin2_mid', 'X1D_peak_bin1', 'X1D_peak_bin2', 'fdr', "dist"]]

inter_all_reformed = pd.concat([inter_all_AND, inter_all_bin1_anchor, inter_all_bin2_anchor])
inter_all_reformed['bin1_mid'] = inter_all_reformed['bin1_mid'].astype('int')
inter_all_reformed['bin2_mid'] = inter_all_reformed['bin2_mid'].astype('int')

inter_all_reformed['bin1_end'] = inter_all_reformed['bin1_mid'] + 5000
inter_all_reformed['bin2_end'] = inter_all_reformed['bin2_mid'] + 5000

inter_all_reformed["chr1"] = inter_all_reformed["chr"]
inter_all_reformed["chr2"] = inter_all_reformed["chr"]

inter_all_reformed = inter_all_reformed.rename(columns={"bin1_mid":"start1", "bin2_mid":"start2", "bin1_end":"end1", "bin2_end":"end2"})
inter_all_reformed = inter_all_reformed[['ID','chr','start1','start2','X1D_peak_bin1','X1D_peak_bin2','fdr','chr1','chr2','end1','end2',"dist"]]

### ---------------------------------------------------------------------------------------------------------------


#%%
####### read in & process significant MAPS interaction sets
inter_sig = pd.read_csv(r'C:\Users\libin\UCSF\hfb\interactions_final\{}.5k.downsample.bedpe'.format(cellType), sep="\t")
inter_sig["sig_ID"] = ["sig_{}_{}".format(cellType, i) for i in range(inter_sig.shape[0])]
print ("number of signif interactions: {}".format(inter_sig.shape[0]))

anchor_bins = pd.read_csv(r'C:\Users\libin\UCSF\hfb\interactions_final\MACS2_peak_5KbBin_Union4.txt', sep="\t")
anchor_bins["anchor_ID"] = ["anchor_{}".format(i) for i in range(anchor_bins.shape[0])]

inter_sig_left_anchor = pd.merge(inter_sig, anchor_bins, left_on=['chr1', 'start1', 'end1'],
                                    right_on = ['chr', 'start', 'end'], how="inner")
left_anchors_list = inter_sig_left_anchor[['sig_ID', 'anchor_ID']]

inter_sig_right_anchor = pd.merge(inter_sig, anchor_bins, left_on=['chr2', 'start2', 'end2'],
                                     right_on = ['chr', 'start', 'end'], how="inner")
right_anchors_list = inter_sig_right_anchor[['sig_ID', 'anchor_ID']]

anchor_anchor_list = pd.merge(left_anchors_list, right_anchors_list, on="sig_ID", how="inner")

inter_sig_XOR = inter_sig[~inter_sig["sig_ID"].isin(anchor_anchor_list["sig_ID"])]
print ("number of signif XOR interactions: {}".format(inter_sig_XOR.shape[0]))

inter_sig_AND = inter_sig[inter_sig["sig_ID"].isin(anchor_anchor_list["sig_ID"])]
inter_sig_AND.loc[:, "anchor"] = "both"

print ("number of signif AND interactions: {}".format(inter_sig_AND.shape[0]))


inter_sig_left_anchor_XOR = pd.merge(inter_sig_XOR, anchor_bins, left_on=['chr1', 'start1', 'end1'],
                                    right_on = ['chr', 'start', 'end'], how="inner")
inter_sig_left_anchor_XOR = inter_sig_left_anchor_XOR[['chr1','start1','end1','chr2','start2','end2','count','expected','fdr','ClusterLabel','ClusterSize',
                                                         'ClusterType','ClusterNegLog10P','ClusterSummit','sig_ID']]
inter_sig_left_anchor_XOR.loc[:, "anchor"] = "left"
print ("number of signif interactions with bin1 as anchor bin: {}".format(inter_sig_left_anchor_XOR.shape[0]))

inter_sig_right_anchor_XOR = pd.merge(inter_sig_XOR, anchor_bins, left_on=['chr2', 'start2', 'end2'],
                                     right_on = ['chr', 'start', 'end'], how="inner")
print ("number of signif interactions with bin2 as anchor bin: {}".format(inter_sig_right_anchor_XOR.shape[0]))
# rename so bin1 is always anchor bin
inter_sig_right_anchor_XOR = inter_sig_right_anchor_XOR.rename(columns={"chr1":"chr2", "start1":"start2", "end1":"end2", "chr2":"chr1", "start2":"start1", "end2":"end1"})
inter_sig_right_anchor_XOR = inter_sig_right_anchor_XOR[['chr1','start1','end1','chr2','start2','end2','count','expected','fdr','ClusterLabel','ClusterSize',
                                                         'ClusterType','ClusterNegLog10P','ClusterSummit','sig_ID']]
inter_sig_right_anchor_XOR.loc[:, "anchor"] = "right"

inter_sig_reformed = pd.concat([inter_sig_AND, inter_sig_left_anchor_XOR, inter_sig_right_anchor_XOR])
inter_sig_reformed = inter_sig_reformed[['chr1','start1','end1','chr2','start2','end2','count','expected','fdr','ClusterLabel','ClusterSize',
                                                         'ClusterType','ClusterNegLog10P','ClusterSummit','sig_ID', "anchor"]]

inter_sig_reformed["dist_sig"] = abs(inter_sig_reformed["start2"] - inter_sig_reformed["start1"])
inter_sig_anchor_list = inter_sig_reformed[['chr1','start1','end1']].drop_duplicates()
#inter_sig_XOR_reformed_anchor = inter_sig_XOR_reformed[['chr1','start1','end1','sig_ID']]
#inter_sig_XOR_reformed_anchor.to_csv(r'C:\Users\libin\UCSF\hfb\eqtl\20200115\interaction_sig_anchor_{}'.format(cellType),sep="\t", index=False, header=False)
#inter_sig_XOR_reformed_target = inter_sig_XOR_reformed[['chr2','start2','end2','sig_ID']]
#inter_sig_XOR_reformed_target.to_csv(r'C:\Users\libin\UCSF\hfb\eqtl\20200115\interaction_sig_target_{}'.format(cellType),sep="\t", index=False, header=False)

#%%
##### set up parameters for distance-matchng
print (inter_sig_reformed["dist_sig"].describe())

distance_bins = list(np.concatenate(
        (np.arange(0, 20000, 10000), np.arange(20000, 80000, 30000), np.arange(80000, 200000, 20000), 
         np.arange(200000, 300000, 25000), np.arange(300000, 500000, 50000), 
         np.arange(500000, 700000, 100000),np.arange(700000, 1000001, 150000)
         )))

inter_sig_reformed_binned = inter_sig_reformed.groupby(pd.cut(inter_sig_reformed["dist_sig"].abs(), distance_bins))
print ("distance distribution (sig interaction): ")
print(inter_sig_reformed_binned.count()["dist_sig"])

sample_size_mdist = [i*int(power) for i in inter_sig_reformed_binned.count()["dist_sig"].tolist()]
#%%


### ---------------------------------------------------------------------------------------------------------------
##### intersect sig set with complete set
inter_all_reformed = inter_all_reformed[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'fdr', 'ID', 'X1D_peak_bin1', 'X1D_peak_bin2', 'dist']]

inter_all_reformed_sig = pd.merge(inter_all_reformed, inter_sig_reformed, how="inner",
                                  on=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'])
# check for duplicates -- should be 0
print ("duplicates:", inter_all_reformed_sig[["sig_ID"]].shape[0] - inter_all_reformed_sig[["sig_ID"]].drop_duplicates().shape[0])
print ("duplicates:", inter_all_reformed_sig[["ID"]].shape[0] - inter_all_reformed_sig[["ID"]].drop_duplicates().shape[0])

# filter for non-sig interactions
inter_all_reformed_sig_ID = inter_all_reformed_sig["ID"]
inter_all_reformed_nonsig = inter_all_reformed[~inter_all_reformed["ID"].isin(inter_all_reformed_sig_ID)]
# take ones with same anchors as in sig 
inter_all_reformed_nonsig_sanchor = pd.merge(inter_sig_anchor_list, inter_all_reformed_nonsig, on=['chr1', 'start1', 'end1'], how="inner")

print ("duplicates:", inter_all_reformed_nonsig_sanchor[["ID"]].shape[0] - inter_all_reformed_nonsig_sanchor[["ID"]].drop_duplicates().shape[0])

### ---------------------------------------------------------------------------------------------------------------

##### sample non-sig interactions
##### distance matched
# sample_size = inter_sig_XOR_reformed.shape[0]

#%%
## check distance distribution first               
inter_all_reformed_nonsig_sanchor_binned = inter_all_reformed_nonsig_sanchor.groupby(pd.cut(inter_all_reformed_nonsig_sanchor["dist"], distance_bins))
print ("distance distribution before sampling (all pairs): ")
print (inter_all_reformed_nonsig_sanchor_binned.count()["dist"])

#%%
inter_all_reformed_nonsig_sanchor_binned_df_list = [group.reset_index().drop(["index"], axis=1) for _, group in inter_all_reformed_nonsig_sanchor_binned]    
inter_all_reformed_nonsig_sanchor_sampled_list = [inter_all_reformed_nonsig_sanchor_binned_df_list[i].sample(n=sample_size_mdist[i]) for i in range(len(sample_size_mdist))]                    
inter_all_reformed_nonsig_sampled_sanchor_concat = pd.concat(inter_all_reformed_nonsig_sanchor_sampled_list).reset_index().drop(["index"], axis=1)
inter_all_reformed_nonsig_sampled_sanchor_binned = inter_all_reformed_nonsig_sampled_sanchor_concat.groupby(pd.cut(inter_all_reformed_nonsig_sampled_sanchor_concat["dist"], distance_bins))
print ("distance distribution after sampling (all pairs): ")
print(inter_all_reformed_nonsig_sampled_sanchor_binned.count()["dist"])

#compare_dist = pd.DataFrame()
#compare_dist["sig"] = inter_sig_XOR_reformed["dist_sig"]
#compare_dist["nonsig"] = inter_all_reformed_nonsig_sampled_concat[["dist"]].reset_index()["dist"]
#compare_dist_melt = pd.melt(compare_dist, value_name="dist", var_name="")

plt.figure(figsize=(10,7))
sns.distplot(inter_sig_reformed["dist_sig"], hist=False, color="blue", label="sig")
sns.distplot(inter_all_reformed_nonsig_sampled_sanchor_concat["dist"], hist=False, color="orange", label="nonsig")
plt.xlabel('distance distribution ({})'.format(cellType), fontsize=15)
plt.savefig(r'C:\Users\libin\UCSF\hfb\eqtl\20200115\{}_dist_compare_{}.pdf'.format(cellType, iteration), transparent=True)


#%%

#sampled_nonsig_anchor = inter_all_reformed_nonsig_sampled_sanchor_concat[['chr1','start1','end1','fdr','ID', "dist"]]
#sampled_nonsig_anchor.to_csv(r'C:\Users\libin\UCSF\hfb\eqtl\20200115\interaction_sampled_sanchor_nonsig_anchor_{}_{}'.format(cellType, iteration), sep="\t", index=False, header=False)
#sampled_nonsig_target = inter_all_reformed_nonsig_sampled_sanchor_concat[['chr2','start2','end2','fdr','ID', "dist"]]
#sampled_nonsig_target.to_csv(r'C:\Users\libin\UCSF\hfb\eqtl\20200115\interaction_sampled_sanchor_nonsig_target_{}_{}'.format(cellType, iteration), sep="\t", index=False, header=False)

inter_all_reformed_nonsig_sampled_sanchor_concat.to_csv(r'C:\Users\libin\UCSF\hfb\eqtl\20200115\interaction_sampled_sanchor_nonsig_{}_{}_{}'.format(cellType, iteration, date), sep="\t", index=False, header=True)












