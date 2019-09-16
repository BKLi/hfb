# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 04:23:23 2019

@author: libin
"""

import pandas as pd
import numpy as np
import seaborn as sns
from functools import reduce
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


coverage_df_lst = []
#eN, iN, IPC, RG
cellTypes = ["MS075", "MS093", "MS119", "MS121", "MS160", "MS164", 
             "MS079", "MS094", "MS106", "MS122", "MS161", "MS165", 
             "MS074", "MS092", "MS104", "MS118", "MS163", "MS166",
             "MS058", "MS072", "MS091", "MS162"]


for ct in cellTypes:
    coverage_df = pd.read_csv(r'C:\Users\libin\UCSF\hfb\ATAC_cov\all_covs\{}.coverage'.format(ct), sep="\t", \
                              names = ["chr","start","end","count","len_read","len_peak","percentage" ])
    coverage_df = coverage_df.drop(["len_read","len_peak","percentage"], axis=1)
    coverage_df = coverage_df.rename(columns={"count" : "{}".format(ct)})
    coverage_df_lst.append(coverage_df)
    
coverage_merged = reduce(lambda x, y : pd.merge(x ,y, on=["chr", "start", "end"]), coverage_df_lst)
coverage_merged = coverage_merged.rename(columns={"MS075":"eN1", "MS093":"eN2", "MS119":"eN3", "MS121":"eN4", "MS160":"eN5", "MS164":"eN6",\
                                                  "MS079":"iN1", "MS094":"iN2", "MS106":"iN3", "MS122":"iN4", "MS161":"iN5", "MS165":"iN6",\
                                                  "MS074":"IPC1", "MS092":"IPC2", "MS104":"IPC3", "MS118":"IPC4", "MS163":"IPC5", "MS166":"IPC6",\
                                                  "MS058":"RG1", "MS072":"RG2", "MS091":"RG3", "MS162":"RG4"})
coverage_merged.to_csv(r'C:\Users\libin\UCSF\hfb\ATAC_cov\ATAC_cov_byRep', sep="\t", index=False, header=True)


coverge_merged_sum = coverage_merged.copy()
coverge_merged_sum["eN"] = coverge_merged_sum.iloc[:,3:9].sum(axis=1)
coverge_merged_sum["iN"] = coverge_merged_sum.iloc[:,9:15].sum(axis=1)
coverge_merged_sum["IPC"] = coverge_merged_sum.iloc[:,15:21].sum(axis=1)
coverge_merged_sum["RG"] = coverge_merged_sum.iloc[:,21:25].sum(axis=1)
coverge_merged_pooled = coverge_merged_sum[["chr","start","end","eN","iN","IPC","RG"]]
coverge_merged_pooled.to_csv(r'C:\Users\libin\UCSF\hfb\ATAC_cov\ATAC_cov_pooled', sep="\t", index=False, header=True)


coverage_merged_cor = coverage_merged.drop(["chr", "start", "end"], axis = 1)
coverage_merged_cor = coverage_merged_cor.corr()
plt.figure(figsize=(15,15))
sns.clustermap(coverage_merged_cor)
plt.savefig(r'C:\Users\libin\UCSF\hfb\ATAC_cov\coverage_corr.pdf', transparent=True)