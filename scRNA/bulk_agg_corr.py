# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 21:53:42 2019

@author: libin
"""

import pandas as pd
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from functools import reduce
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

large_types = pd.read_csv("C:\\Users\libin\\UCSF\hfb\scRNA\\AllCellSubTypes.csv", sep=",").rename(columns={"Unnamed: 0":"gene_name"})
# all_types = pd.read_table("C:\\Users\libin\\UCSF\hfb\scRNA\\AllCellSubTypes.csv", sep=",").rename(columns={"Unnamed: 0":"gene_name"})

large_type_ct = large_types.columns.drop("gene_name",1).tolist()
for ct1 in large_type_ct:
    large_types[ct1] = large_types[ct1].apply(lambda x: np.log2(x+1))
    
bulk_IN = pd.read_table(r"C:\Users\libin\UCSF\hfb\scRNA\updated_Aug_27\\iNs.average.genes.results", sep="\t")\
[["gene_name","iN"]]\
.set_index("gene_name")\
[["iN"]].apply(lambda x : np.log2(x+1))
# .rename(columns={"TPM":"IN"})\

bulk_IPC = pd.read_table(r"C:\Users\libin\UCSF\hfb\scRNA\updated_Aug_27\\IPCs.average.genes.results", sep="\t")\
[["gene_name","IPC"]]\
.set_index("gene_name")\
[["IPC"]].apply(lambda x : np.log2(x+1))
# .rename(columns={"TPM":"IPC"})\

bulk_neuron = pd.read_table(r"C:\Users\libin\UCSF\hfb\scRNA\updated_Aug_27\\eNs.average.genes.results", sep="\t")\
[["gene_name","eN"]]\
.set_index("gene_name")\
[["eN"]].apply(lambda x : np.log2(x+1))
# .rename(columns={"TPM":"neuron"})\

bulk_RGC = pd.read_table(r"C:\Users\libin\UCSF\hfb\scRNA\updated_Aug_27\\RGs.average.genes.results", sep="\t")\
[["gene_name","RG"]]\
.set_index("gene_name")\
[["RG"]].apply(lambda x : np.log2(x+1))
# .rename(columns={"TPM":"RGC"})\

# merge four dfs on TPM
# https://stackoverflow.com/questions/23668427/pandas-three-way-joining-multiple-dataframes-on-columns
# dfs = [bulk_IN, bulk_IPC, bulk_neuron, bulk_RGC]
# bulk_joined = reduce(lambda x,y: pd.merge(x,y,on="gene_name",how="inner"), dfs)
bulk_joined = pd.concat([bulk_IN,bulk_IPC,bulk_neuron,bulk_RGC], axis=1).reset_index()
bulk_large_types = pd.merge(bulk_joined, large_types, on="gene_name", how="inner").set_index("gene_name")
bulk_large_types["eN_mean"] = \
bulk_large_types[
        ['eN1',
         'eN2',
         'eN3',
         'eN4',
         'eN5',
         'eN6',
         'eN7',
         'eN8',
         'eN9' ]].mean(axis=1)

bulk_large_types["iN_mean"] = \
bulk_large_types[ 
['iN1',
 'iN2',
 'iN3',
 'iN4',
 'iN5']].mean(axis=1)

bulk_large_types["NiN_mean"] = \
bulk_large_types[ 
['NiN1',
 'NiN2',
 'NiN3',
 'NiN4',
 'NiN5']].mean(axis=1)

bulk_large_types["IPC_mean"] = \
bulk_large_types[ 
['IPC1',
 'IPC2',
 'IPC3',
 'IPC4',
 'IPC5']].mean(axis=1)

bulk_large_types["RGC_mean"] = \
bulk_large_types[ 
['RG1',
 'RG2',
 'RG3',
 'RG4',
 'RG5',
 'RG6',
 'RG7']].mean(axis=1)

bulk_large_types["Endothelial"] = \
bulk_large_types[ 
['other1']]

bulk_large_types["Mural"] = \
bulk_large_types[ 
['other2']]

bulk_large_types["Microglia"] = \
bulk_large_types[ 
['other5']]

bulk_large_types["Choroid"] = \
bulk_large_types[ 
['other7']]

bulk_large_types["ERG"] = \
bulk_large_types[ 
['other8']]

bulk_large_types["VP_mean"] = \
bulk_large_types[ 
['VP1',
 'VP2',
 'VP3',
 'VP4',
 'VP5',
 'VP6']].mean(axis=1)

bulk_ctypes = ["iN", "IPC", "RG", "eN"]

# plt.figure(figsize=(20,20))
# sns.clustermap(bulk_large_types.corr(method="pearson"), row_cluster=False, cmap="RdBu_r",xticklabels=True,yticklabels=True)
# plt.savefig('C:\\Users\\libin\\UCSF\\hfb\\scRNA\\bulk_large_type_corr.pdf', transparent=True, bbox_inches = 'tight',
#    pad_inches = 0.1)

corr = bulk_large_types.corr(method="pearson")
corr_part = corr.drop(["iN", "IPC", "RG", "eN"], axis=0)
corr_part = corr_part.loc[:, ["RG", "IPC", "eN", "iN"]]
## !!! March 3rd
corr_part = corr_part.loc[['RGC_mean', 'IPC_mean', 'eN_mean', 'iN_mean', 'NiN_mean', 
                        'VP_mean', 'ERG', 'Endothelial', 'Mural', 'Microglia', 'Choroid'],:].rename({'RGC_mean': "  RG  ", 'IPC_mean':"  IPC  ", 'eN_mean':"  eN  ", 'iN_mean':"  iN  ", 'NiN_mean':"  NiN  ", 
                         'VP_mean':"  VP  "}, axis='index')

corr_part.to_csv(r"C:\Users\libin\UCSF\hfb\scRNA\updated_Aug_27\\mean_corr.csv", sep=",")
# corr_norm_row=corr_part.sub(corr_part.mean(axis=1), axis=0)
# corr_norm_row=corr_norm_row.div(corr_part.std(axis=1), axis=0)
# corr_norm_col=(corr_part-corr_part.mean())/corr_part.std()
'''corr_part_IN = corr_part[["IN"]]
corr_part_IN = corr_part_IN.sort_values(by=["IN"], ascending=False)
plt.figure(figsize=(5,7))
sns.heatmap(corr_part_IN, cmap="RdBu_r", yticklabels=True)
plt.savefig('C:\\Users\\libin\\UCSF\\hfb\\scRNA\\heatmap_IN-mean_allTypes.pdf', transparent=True, bbox_inches = 'tight',
    pad_inches = 0.1)

corr_part_IPC = corr_part[["IPC"]]
corr_part_IPC = corr_part_IPC.sort_values(by=["IPC"], ascending=False)
plt.figure(figsize=(5,7))
sns.heatmap(corr_part_IPC, cmap="RdBu_r", yticklabels=True)
plt.savefig('C:\\Users\\libin\\UCSF\\hfb\\scRNA\\heatmap_IPC_mean_allTypes.pdf', transparent=True, bbox_inches = 'tight',
    pad_inches = 0.1)

corr_part_RGC = corr_part[["RGC"]]
corr_part_RGC = corr_part_RGC.sort_values(by=["RGC"], ascending=False)
plt.figure(figsize=(5,7))
sns.heatmap(corr_part_RGC, cmap="RdBu_r", yticklabels=True)
plt.savefig('C:\\Users\\libin\\UCSF\\hfb\\scRNA\\heatmap_RGC_mean_allTypes.pdf', transparent=True, bbox_inches = 'tight',
    pad_inches = 0.1)

corr_part_neuron = corr_part[["neuron"]]
corr_part_neuron = corr_part_neuron.sort_values(by=["neuron"], ascending=False)
plt.figure(figsize=(5,7))
sns.heatmap(corr_part_neuron, cmap="RdBu_r", yticklabels=True)
plt.yticks(rotation=0) 
plt.savefig('C:\\Users\\libin\\UCSF\\hfb\\scRNA\\heatmap_neuron_mean_allTypes.pdf', transparent=True, bbox_inches = 'tight',
    pad_inches = 0.1)

# find genes co-expressed in IN and eNs
IN_eN4 = bulk_large_types[["IN","eN4"]]
IN_eN4_filtered = IN_eN4[(IN_eN4["IN"] > 0) & (IN_eN4["eN4"] > 0)]\
.sort_values(["IN","eN4"], ascending=[False, False])'''

plt.figure(figsize=(5,5))
sns.heatmap(corr_part, cmap="RdBu_r", yticklabels=True)
plt.yticks(rotation=0) 
plt.savefig(r'C:\Users\libin\UCSF\hfb\scRNA\updated_Aug_27\\heatmap_average_all_fullTypes.pdf', transparent=True, bbox_inches = 'tight',
    pad_inches = 0.1)


