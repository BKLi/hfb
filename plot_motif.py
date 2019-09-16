# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 17:07:08 2019

@author: libin
"""

import pandas as pd
import numpy as np
import seaborn as sns
from functools import reduce

cellType = ['RG', 'IPC', 'eN', 'iN']
homer_extra = ["Eomes", "Tbr1", "PAX6", "Brn2", "Brn1", "Sox2"]
jaspar_extra = ["EOMES", "TBR1", "Pax6", "DLX6", "POU3F2", "POU3F3", "Sox2"]

# EOMES : JASPAR,HOMER
# TBR1: JASPAR, HOMER
# PAX6: JASPAR, HOMER
# DLX5: NA
# DLX6: JASPAR
# POU3F2: JASPAR; BRN2: HOMER
# POU3F3: JASPAR; BRN1: HOMER
# SOX2: HOMER, JASPAR


jaspar_df_lst = []
jaspar_cleaned_df_lst = []
jaspar_merged_df_lst = []
jaspar_motif_list = []
for ct in cellType:
# ct = 'eN'
    jaspar = pd.read_csv(r'C:\Users\libin\UCSF\motif_analysis\Aug 20\jaspar_known_results\knownResults_{}.txt'.format(ct), sep="\t")
    jaspar["CellType"] = ct
    jaspar["Motif"] = jaspar['Motif Name'].str.extract(r'(.+?)[\/]', expand=True)
    # jaspar["Motif"] = jaspar['Motif Name'].str.upper()
    jaspar = jaspar[["Motif", "P-value"]]
    jaspar = jaspar.rename(columns={"P-value" : "{}".format(ct)})
    jaspar_df_lst.append(jaspar)
    jaspar_cleaned = jaspar.iloc[:15,:]
    jaspar_cleaned_df_lst.append(jaspar_cleaned)
    jaspar_motif_list = jaspar_motif_list +  jaspar_cleaned['Motif'].tolist()
jaspar_motif_list = jaspar_motif_list + jaspar_extra
jaspar_motif_list_dedup = list(dict.fromkeys(jaspar_motif_list))
jaspar_motif_list_dedup = pd.DataFrame(jaspar_motif_list_dedup, columns=["Motif"])
for jaspar_df in jaspar_df_lst:
    jaspar_merged = pd.merge(jaspar_motif_list_dedup, jaspar_df, on=["Motif"], how="inner")
    jaspar_merged_df_lst.append(jaspar_merged)
jaspar_merged_all = reduce(lambda x, y : pd.merge(x ,y, on=["Motif"]), jaspar_merged_df_lst)
jaspar_merged_all.to_csv(r'C:\Users\libin\R_projects\motif_analysis\jaspar_merged', index=False, header=True, sep="\t")
    

homer_df_lst = []
homer_cleaned_df_lst = []
homer_merged_df_lst = []
homer_motif_list = []
for ct in cellType:
    homer = pd.read_csv(r'C:\Users\libin\UCSF\motif_analysis\Aug_21\homer_known_results\knownResults_{}.txt'.format(ct), sep="\t")
    homer["CellType"] = ct
    # remove 'SeqBias' motifs
    homer = homer[~homer["Motif Name"].str.contains("Bias")]
    homer["Motif"] = homer['Motif Name'].str.extract(r'(.+?)(?:\/|$)', expand=True)
    # homer["Motif"] = homer['Motif Name'].str.upper()
    homer["Motif_cleaned"] = homer['Motif'].str.extract(r'(.+?)(?:\(|$)', expand=True)
    homer = homer[["Motif_cleaned", "P-value"]]
    homer = homer.rename(columns={"P-value" : "{}".format(ct), "Motif_cleaned" : "Motif"})
    homer_df_lst.append(homer)
    homer_cleaned = homer.iloc[:15,:]
    homer_cleaned_df_lst.append(homer_cleaned)
    homer_motif_list = homer_motif_list +  homer_cleaned['Motif'].tolist()
homer_motif_list = homer_motif_list + homer_extra
homer_motif_list_dedup = list(dict.fromkeys(homer_motif_list))
homer_motif_list_dedup = pd.DataFrame(homer_motif_list_dedup, columns=["Motif"])
for homer_df in homer_df_lst:
    homer_merged = pd.merge(homer_motif_list_dedup, homer_df, on=["Motif"], how="inner")
    homer_merged_df_lst.append(homer_merged)
homer_merged_all = reduce(lambda x, y : pd.merge(x ,y, on=["Motif"]), homer_merged_df_lst)
homer_merged_all.to_csv(r'C:\Users\libin\R_projects\motif_analysis\homer_merged', index=False, header=True, sep="\t")
