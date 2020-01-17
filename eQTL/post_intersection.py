# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 02:02:34 2020

@author: bingkun

@ project hfb -- eQTL

* Data processing after intersecting interactions with eQTL-gene pairs using bedtools
* S4 of pipeline
* previous step : intersecting with bedtools
"""


#%%
#import sys
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['figure.dpi']= 200
#%%
# def post_intersection(raw_inter_path, dir_path, cellType, date, run):
raw_inter_path = r'C:\Users\libin\UCSF\MAGMA'
dir_path = r'C:\Users\libin\UCSF\hfb\eqtl\20191227'
criteria = "count"
self_ctrl_prefix = "1X_mdist_SigGeneControl"
cellType= "neuron"
control_cellType = "Lung"
neg_control_prefix = "nonsig.mdist.1X"
date = "0109"
run = "6"

# read in raw interactions
raw_inter_all = pd.read_csv(r'{}\{}.5k.downsample.bedpe.withID'.format(raw_inter_path, cellType), sep="\t")

# read in signif intersections
anchor_gene_intersection_sig = pd.read_csv(r'{}\{}_anchor_gene.intersect'.format(dir_path ,cellType), sep="\t",
                                       names=["chr_ac", "start_ac", "end_ac", "inter_id", "anchor_id", "chr_gene", "start_gene", "end_gene", "gene_id", "npval", "pair_id"])
target_snp_intersection_sig = pd.read_csv(r'{}\{}_target_snp.intersect'.format(dir_path, cellType), sep="\t",
                                      names=["chr_tg", "start_tg", "end_tg", "inter_id", "anchor_id", "chr_snp", "start_snp", "end_snp", "snp_id", "npval", "pair_id"])

interaction_eqtlPair_intersection_sig = pd.merge(anchor_gene_intersection_sig, target_snp_intersection_sig, on=["inter_id", "anchor_id", "npval", "pair_id"], how="inner")
print ("number of total interactions involved in sig intersection: {}".format(interaction_eqtlPair_intersection_sig.shape[0]))
# take all interaction IDs appeared in sig intersection
inter_ID_sig = interaction_eqtlPair_intersection_sig[["inter_id"]].drop_duplicates()
inter_sig = pd.merge(inter_ID_sig, raw_inter_all, left_on="inter_id", right_on="ID")
inter_sig["fdr"] = 0-inter_sig["fdr"].apply(np.log2)
inter_sig["norm_count"] = inter_sig["count"]/inter_sig["expected"]
if criteria == "fdr":
        print (inter_sig["{}".format(criteria)].describe())
else:
        print (inter_sig["norm_count"].describe())
        
# take IDs of eqtl-gene pairs that intersected with interactions
inter_pairID = interaction_eqtlPair_intersection_sig["pair_id"]

####------------- self control

anchor_gene_intersection_control = pd.read_csv(r'{}\{}_anchor_gene.{}.intersect'.format(dir_path, cellType, self_ctrl_prefix), sep="\t",
                                       names=["chr_ac", "start_ac", "end_ac", "inter_id", "anchor_id", "chr_gene", "start_gene", "end_gene", "gene_id", "npval", "pair_id"])
target_snp_intersection_control = pd.read_csv(r'{}\{}_target_snp.{}.intersect'.format(dir_path, cellType, self_ctrl_prefix), sep="\t",
                                      names=["chr_tg", "start_tg", "end_tg", "inter_id", "anchor_id", "chr_snp", "start_snp", "end_snp", "snp_id", "npval", "pair_id"])

interaction_eqtlPair_intersection_control = pd.merge(anchor_gene_intersection_control, target_snp_intersection_control, on=["inter_id", "anchor_id", "npval", "pair_id"], how="inner")
print ("number of total interactions involved in control intersection: {}".format(interaction_eqtlPair_intersection_control.shape[0]))

# take all interaction IDs appeared in control intersection
inter_ID_control = interaction_eqtlPair_intersection_control[["inter_id"]].drop_duplicates()
inter_control = pd.merge(inter_ID_control, raw_inter_all, left_on="inter_id", right_on="ID")
inter_control["fdr"] = 0-inter_control["fdr"].apply(np.log2)
inter_control["norm_count"] = inter_control["count"]/inter_control["expected"]
if criteria == "fdr":
        print (inter_control["{}".format(criteria)].describe())
else:
        print (inter_control["norm_count"].describe())

# take IDs of eqtl-gene pairs that intersected with interactions
inter_control_pairID = interaction_eqtlPair_intersection_control["pair_id"]


#### ------------------ negative control -- sig

anchor_gene_intersection_NegControl_sig = pd.read_csv(r'{}\{}_anchor_{}_gene.sig.intersect'.format(dir_path, cellType, control_cellType), sep="\t",
                                       names=["chr_ac", "start_ac", "end_ac", "inter_id", "anchor_id", "chr_gene", "start_gene", "end_gene", "pair_id", "gene_id", "npval"])
target_snp_intersection_NegControl_sig = pd.read_csv(r'{}\{}_target_{}_snp.sig.intersect'.format(dir_path, cellType, control_cellType), sep="\t",
                                      names=["chr_tg", "start_tg", "end_tg", "inter_id", "anchor_id", "chr_snp", "start_snp", "end_snp", "snp_id", "pair_id", "npval"])

interaction_eqtlPair_intersection_NegControl_sig = pd.merge(anchor_gene_intersection_NegControl_sig, target_snp_intersection_NegControl_sig, on=["inter_id", "anchor_id", "npval", "pair_id"], how="inner")
print ("number of total interactions involved in sig negative control intersection: {}".format(interaction_eqtlPair_intersection_NegControl_sig.shape[0]))

# take all interaction IDs appeared in sig neg control intersection
inter_ID_NegControl_sig = interaction_eqtlPair_intersection_NegControl_sig[["inter_id"]].drop_duplicates()
inter_NegControl_sig = pd.merge(inter_ID_NegControl_sig, raw_inter_all, left_on="inter_id", right_on="ID")
inter_NegControl_sig["fdr"] = 0-inter_NegControl_sig["fdr"].apply(np.log2)
inter_NegControl_sig["norm_count"] = inter_NegControl_sig["count"]/inter_NegControl_sig["expected"]
if criteria == "fdr":
        print (inter_NegControl_sig["{}".format(criteria)].describe())
else:
        print (inter_NegControl_sig["norm_count"].describe())
        
# take IDs of eqtl-gene pairs that intersected with interactions
inter_NegControl_sig_pairID = interaction_eqtlPair_intersection_NegControl_sig["pair_id"]

#### ----------------- negative control -- nonsig

anchor_gene_intersection_NegControl_nonsig = pd.read_csv(r'{}\{}_anchor_{}_gene.{}.intersect'.format(dir_path, cellType, control_cellType, neg_control_prefix), sep="\t",
                                       names=["chr_ac", "start_ac", "end_ac", "inter_id", "anchor_id", "chr_gene", "start_gene", "end_gene", "pair_id", "gene_id", "npval"])
target_snp_intersection_NegControl_nonsig = pd.read_csv(r'{}\{}_target_{}_snp.{}.intersect'.format(dir_path, cellType, control_cellType, neg_control_prefix), sep="\t",
                                      names=["chr_tg", "start_tg", "end_tg", "inter_id", "anchor_id", "chr_snp", "start_snp", "end_snp", "snp_id", "pair_id", "npval"])

interaction_eqtlPair_intersection_NegControl_nonsig = pd.merge(anchor_gene_intersection_NegControl_nonsig, target_snp_intersection_NegControl_nonsig, on=["inter_id", "anchor_id", "npval", "pair_id"], how="inner")
print ("number of total interactions involved in negative control intersection: {}".format(interaction_eqtlPair_intersection_NegControl_nonsig.shape[0]))
# take all interaction IDs appeared in nonsig neg control intersection
inter_ID_NegControl_nonsig = interaction_eqtlPair_intersection_NegControl_nonsig[["inter_id"]].drop_duplicates()
inter_NegControl_nonsig = pd.merge(inter_ID_NegControl_nonsig, raw_inter_all, left_on="inter_id", right_on="ID")
inter_NegControl_nonsig["fdr"] = 0-inter_NegControl_nonsig["fdr"].apply(np.log2)
inter_NegControl_nonsig["norm_count"] = inter_NegControl_nonsig["count"]/inter_NegControl_nonsig["expected"]
if criteria == "fdr":
        print (inter_NegControl_nonsig["{}".format(criteria)].describe())
else:
        print (inter_NegControl_nonsig["norm_count"].describe())
        
# take IDs of eqtl-gene pairs that intersected with interactions
inter_NegControl_nonsig_pairID = interaction_eqtlPair_intersection_NegControl_nonsig["pair_id"]


### make plot
y_ax = criteria
compare_sig_control = pd.DataFrame()

if criteria == "fdr":
        compare_sig_control["{}_sig".format(cellType)] = inter_sig[criteria]
        compare_sig_control["{}_control".format(cellType)] = inter_control[criteria]
        compare_sig_control["{}_sig".format(control_cellType)] = inter_NegControl_sig[criteria]
        compare_sig_control["{}_nonsig".format(control_cellType)] = inter_NegControl_nonsig[criteria]

else:
        compare_sig_control["{}_sig".format(cellType)] = inter_sig["norm_count"]
        compare_sig_control["{}_control".format(cellType)] = inter_control["norm_count"]
        compare_sig_control["{}_sig".format(control_cellType)] = inter_NegControl_sig["norm_count"]
        compare_sig_control["{}_nonsig".format(control_cellType)] = inter_NegControl_nonsig["norm_count"]

        
compare_sig_control_melt = pd.melt(compare_sig_control, value_name=y_ax, var_name="")

print("sig vs inter_control: ", stats.ks_2samp(inter_sig[criteria], inter_control[criteria]))
print("sig vs neg_control_sig: ", stats.ks_2samp(inter_sig[criteria], inter_NegControl_sig[criteria]))
print("sig vs neg_control_nonsig: ", stats.ks_2samp(inter_sig[criteria], inter_NegControl_nonsig[criteria]))

# compare_sig_control_melt.to_csv(r'{}\{}_{}_{}_run{}.melt'.format(dir_path, cellType, control_cellType, date, run), sep="\t", index=False, header=True)

plt.figure(figsize=(10,7))
sns.violinplot(x="", y=y_ax, data=compare_sig_control_melt, palette="Pastel1", cut=0)
plt.savefig(r'{}\{}_{}_{}_run{}_{}.pdf'.format(dir_path, cellType, control_cellType, date, run, criteria), transparent=True)

### ---------------------------------------------------------------------------

        #%%
plot_distance = True
if plot_distance:
        ### eQTL distance distributions
        distance_distribution = pd.read_csv(r'C:\Users\libin\UCSF\hfb\eqtl\20191227\Lung_sig_gene_shuf_dist_matched_for_plot', sep="\t")
        #distance_distribution = distance_distribution[['rand_gene_all', 'rand_gene_nonsig', 'sig']]

        distance_distribution_melt = pd.melt(distance_distribution, value_name="dist", var_name="")
        
        plt.figure(figsize=(10,7))
        sns.violinplot(x="", y="dist", data=distance_distribution_melt, palette="Pastel1", cut=0)
        plt.savefig(r'C:\Users\libin\UCSF\hfb\eqtl\20191227\Lung_sig_gene_shuf_dist_matched_for_plot.pdf', transparent=True)

        print (distance_distribution["sig"].describe())
        print (distance_distribution["rand_gene_all"].describe())
        print (distance_distribution["rand_gene_nonsig"].describe())
        
        distance_bins = list(np.concatenate(
                (np.arange(10000, 80000, 10000), np.arange(80000, 200000, 20000), np.arange(200000, 300000, 50000), np.arange(300000, 600000, 150000), np.arange(600000, 1000001, 400000)
                 )))
        
        sig_eqtl_binned = distance_distribution[["sig"]].groupby(pd.cut(distance_distribution["sig"], distance_bins))
        print(sig_eqtl_binned.count()["sig"])
        
        sig_eqtl_binned = distance_distribution[["rand_gene_all"]].groupby(pd.cut(distance_distribution["rand_gene_all"], distance_bins))
        print(sig_eqtl_binned.count()["rand_gene_all"])
        
        sig_eqtl_binned = distance_distribution[["rand_gene_nonsig"]].groupby(pd.cut(distance_distribution["rand_gene_nonsig"], distance_bins))
        print(sig_eqtl_binned.count()["rand_gene_nonsig"])

# post_intersection(raw_inter_path=sys.argv[1], dir_path=sys.argv[2], cellType=sys.argv[3], date=sys.argv[4], run=sys.argv[5])