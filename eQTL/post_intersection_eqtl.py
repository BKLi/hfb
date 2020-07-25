# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 02:15:22 2020

@author: bingkun
@project: hfb -- eqtl

* process of eqtl tss/snp region after intersecting with interactions
"""

import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['figure.dpi']= 200


# def post_intersect_eqtl(sig_tss, sig_snp, nonsig_tss, nonsig_snp, cellType):
date = "0115"        
cellType = "RGC"

sig_tss = r'C:\Users\libin\UCSF\hfb\eqtl\20200115\{}_sig_intersect_etql.tss'.format(cellType)
sig_snp = r'C:\Users\libin\UCSF\hfb\eqtl\20200115\{}_sig_intersect_etql.snp'.format(cellType)
nonsig_tss = r'C:\Users\libin\UCSF\hfb\eqtl\20200115\{}_nonsig_intersect_etql.tss'.format(cellType)
nonsig_snp = r'C:\Users\libin\UCSF\hfb\eqtl\20200115\{}_nonsig_intersect_etql.snp'.format(cellType)

sig_tss_df = pd.read_csv(sig_tss, sep="\t", names=["chr", "start_tss", "end_tss", "pair_id", "gene_id", "npval", "dist"])
sig_snp_df = pd.read_csv(sig_snp, sep="\t", names=["chr", "start_snp", "end_snp", "pair_id", "gene_id", "npval", "dist"])
nonsig_tss_df = pd.read_csv(nonsig_tss, sep="\t", names=["chr", "start_tss", "end_tss", "pair_id", "gene_id", "npval", "dist"])
print ("nonsig tss duplicates: ", nonsig_tss_df[["pair_id"]].shape[0] - nonsig_tss_df[["pair_id"]].drop_duplicates().shape[0])
nonsig_snp_df = pd.read_csv(nonsig_snp, sep="\t", names=["chr", "start_snp", "end_snp", "pair_id", "gene_id", "npval", "dist"])
print ("nonsig snp duplicates: ", nonsig_snp_df[["pair_id"]].shape[0] - nonsig_snp_df[["pair_id"]].drop_duplicates().shape[0])
#nonsig_snp_df_dup_bool = nonsig_snp_df.duplicated(subset=['pair_id'], keep=False)
#nonsig_snp_df_dup = nonsig_snp_df.loc[nonsig_snp_df_dup_bool == True]

sig_merged = pd.merge(sig_tss_df, sig_snp_df, on = ["chr", "pair_id", "gene_id", "npval", "dist"])
# check dups -- should be 0
print ("sig duplicates: ", sig_merged[["pair_id"]].shape[0] - sig_merged[["pair_id"]].drop_duplicates().shape[0])
nonsig_merged = pd.merge(nonsig_tss_df, nonsig_snp_df, on = ["chr", "pair_id", "gene_id", "npval", "dist"])
#nonsig_merged_dup_bool = nonsig_merged.duplicated(subset=['pair_id'], keep=False)
#nonsig_merged_dup = nonsig_merged.loc[nonsig_merged_dup_bool == True]
print ("nonsig duplicates: ", nonsig_merged[["pair_id"]].shape[0] - nonsig_merged[["pair_id"]].drop_duplicates().shape[0])
print ('removing duplicates...')
sig_merged = sig_merged.drop_duplicates(subset=['pair_id'])
nonsig_merged = nonsig_merged.drop_duplicates(subset=['pair_id'])

npval_compair = pd.DataFrame()
npval_compair["{}_non_signif".format(cellType)] = nonsig_merged["npval"]
npval_compair["{}_signif".format(cellType)] = sig_merged["npval"]
npval_compair.to_csv("npval_compair_{}".format(cellType), sep="\t", index=False, header=True)
npval_compair_melt = pd.melt(npval_compair, value_name="npval", var_name="")

print ("sig pval: ", sig_merged["npval"].describe())
print ("nonsig pval: ", nonsig_merged["npval"].describe())
print(stats.ks_2samp(sig_merged["npval"], nonsig_merged["npval"]))

plt.figure(figsize=(10,7))
sns.violinplot(x="", y="npval", data=npval_compair_melt, palette="Pastel1", cut=0)
plt.savefig(r'C:\Users\libin\UCSF\hfb\eqtl\20200115\{}_{}_pval_compare.pdf'.format(cellType, date), transparent=True)

# post_intersect_eqtl(sig_tss=sys.argv[1], sig_snp=sys.argv[2], nonsig_tss=sys.argv[3], nonsig_snp=sys.argv[4], cellType=sys.argv[5])