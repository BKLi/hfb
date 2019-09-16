# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 14:42:43 2019

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

fpkm = pd.read_table("C:\\Users\\libin\\R_projects\\DESeq2\\rpkm.data.individual.txt", sep="\t")
translation = pd.read_table("C:\\Users\\libin\\R_projects\\DESeq2\\promoters.strict.txt", sep="\t")
                            # names=["chrom","chromStart","chromEnd","strand","gene_id","gene_type","gene_name"])
                            
# keep only protein cos=ding genes to avoid multiple messy id-name correlation
translation = translation[translation["gene_type"] == "protein_coding"]

# translation = translation
translation_dedup = translation.drop_duplicates(subset=['gene_id','gene_type','gene_name'], keep="first")
# remove chrY par genes to avoid duplication in gene_name
translation_dedup = translation_dedup[~translation_dedup["gene_id"].str.contains("ENSGR")]

fpkm.index.name = "gene_id"
fpkm.reset_index(inplace=True)

# fpkm_outjoin = pd.merge(fpkm, translation_dedup, on="gene_id", how="outer")
# fpkm_debug = fpkm_outjoin[fpkm_outjoin.isnull().any(axis=1)]

fpkm_innerjoin = pd.merge(fpkm, translation_dedup, on ="gene_id", how="inner")
# rename for merging later
fpkm_innerjoin = fpkm_innerjoin.rename(index=str, columns={"gene_name": 'gene'})
# remove unwanted columns
fpkm_innerjoin = fpkm_innerjoin[['gene','gene_id','cortical']]

# raw read matrix
sc_expression = pd.read_table("C:\\Users\\libin\\UCSF\\hfb\\scRNA\\exprMatrix.tsv",sep="\t")
fpkm_sc = pd.merge(fpkm_innerjoin, sc_expression, on="gene", how="inner")
fpkm_sc = fpkm_sc.set_index("gene")
# calculate pairwaise correations
sc_list = fpkm_sc.columns.tolist()[2:]
sc_for_calc = fpkm_sc.drop(columns=["gene_id","cortical"])
correlations = {}
for cell in sc_list: correlations[cell] = pearsonr(fpkm_sc.loc[:, "cortical"].apply(lambda x: np.log2(x+1)), sc_for_calc.loc[:, cell].apply(lambda y: np.log2(y+1)))
correlation_df = pd.DataFrame.from_dict(correlations, orient="index")
correlation_df.columns = ["PCC", "p-value"]

projection = pd.read_table("C:\\Users\\libin\\UCSF\\hfb\\scRNA\\TSNE Projection_C1Data.csv", sep=",", index_col=0)
projection_pcc = pd.merge(projection, correlation_df, left_index=True, right_index=True, how="inner")

# binning PCC into discrete categories
# steps = (projection_pcc["PCC"].max()-abs(projection_pcc["PCC"].min()))/20
# bins = np.arange(projection_pcc["PCC"].min(), projection_pcc["PCC"].max(), steps)
# projection_pcc["PCCrange"] = pd.cut(projection_pcc["PCC"], bins)
projection_pcc["PCC_bin"] = pd.qcut(projection_pcc["PCC"], 10)

plt.figure(figsize=(15,9))
sns.scatterplot(x="tSNE_1", y="tSNE_2", hue="PCC_bin", data=projection_pcc, palette="RdBu_r", s=150, marker=".")
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig('C:\\Users\\libin\\UCSF\\hfb\\scRNA\\tsne_corr_pchic_log_pearson.pdf', transparent=True, bbox_inches = 'tight',
    pad_inches = 0.1)

# projection_pcc.to_csv("C:\\Users\\libin\\UCSF\\hfb\\scRNA\\projection_PCC.csv", sep=",")
# fpkm_sc.to_csv("C:\\Users\\libin\\UCSF\\hfb\\scRNA\\bulk_sc_expression.csv", sep=",")
