# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 18:32:03 2019

@author: libin
"""

import pandas as pd
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

tpm = pd.read_table("C:\\Users\libin\\UCSF\hfb\scRNA\\IPCs.average.genes.results", sep="\t")
tpm = tpm.rename(index=str, columns={"gene_name": "gene"})
tpm = tpm[["gene",'TPM']]

# raw_read matrix
sc_expression = pd.read_table("C:\\Users\\libin\\UCSF\\hfb\\scRNA\\exprMatrix.tsv",sep="\t")
tpm_sc = pd.merge(tpm, sc_expression, on="gene", how="inner")
tpm_sc = tpm_sc.set_index("gene")

# calculate pairwaise correations
sc_list = tpm_sc.columns.tolist()[1:]
# sc_for_calc = tpm_sc.drop(columns=["gene_id", "cortical"])
sc_for_calc = tpm_sc
correlations = {}
for cell in sc_list: correlations[cell] = pearsonr(sc_for_calc.loc[:, "TPM"].apply(lambda x: np.log2(x+1)), sc_for_calc.loc[:, cell].apply(lambda y: np.log2(y+1)))
# for cell in sc_list: correlations[cell] = pearsonr(sc_for_calc.loc[:, "TPM"], sc_for_calc.loc[:, cell])
correlation_df = pd.DataFrame.from_dict(correlations, orient="index")
correlation_df.columns = ["PCC", "p-value"]

projection = pd.read_table("C:\\Users\\libin\\UCSF\\hfb\\scRNA\\TSNE Projection_C1Data.csv", sep=",", index_col=0)
projection_pcc = pd.merge(projection, correlation_df, left_index=True, right_index=True, how="inner")

# binning PCC into discrete categories
# steps = (projection_pcc["PCC"].max()-abs(projection_pcc["PCC"].min()))/20
# bins = np.arange(projection_pcc["PCC"].min(), projection_pcc["PCC"].max(), steps)
# projection_pcc["PCCrange"] = pd.cut(projection_pcc["PCC"], bins)
projection_pcc["PCC_bin"] = pd.qcut(projection_pcc["PCC"], 30)

plt.figure(figsize=(15,9))
sns.scatterplot(x="tSNE_1", y="tSNE_2", hue="PCC_bin", data=projection_pcc, palette="RdBu_r", s=150, marker=".", linewidth=0)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig('C:\\Users\\libin\\UCSF\\hfb\\scRNA\\tsne_corr_IPCs_pearson_recolor.pdf', transparent=True, bbox_inches = 'tight',
    pad_inches = 0.1)

# projection_pcc.to_csv("C:\\Users\\libin\\UCSF\\hfb\\scRNA\\projection_PCC.csv", sep=",")
# tpm_sc.to_csv("C:\\Users\\libin\\UCSF\\hfb\\scRNA\\bulk_sc_expression.csv", sep=",")