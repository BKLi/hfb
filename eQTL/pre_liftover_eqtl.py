# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 18:20:28 2020

@author: bingkun
@project: hfb -- eQTL

format eQTL dataset from Geschwind Lab / prepare for lifting over to hg38
* S1 of alternative pipeline
* next step: eqtl_hg19tohg38.R
"""

import sys
import pandas as pd
import numpy as np


def reformat_eqtl(infile):
        
        eqtl_hg19 = pd.read_csv(infile, delim_whitespace=True)
        print ("number of total eqtl-gene pair: {}".format(eqtl_hg19.shape[0]))
        
        # add ID to eqtl-gene pairs
        eqtl_hg19["pairid"] = ["pair_{}".format(i) for i in range(eqtl_hg19.shape[0])]
        
        # filter out snp-gene pair < 10kb
        filtered_eqtl_hg19 = eqtl_hg19[abs(eqtl_hg19["dist"]) > 10000]
        print ("number of pairs with dist > 10kb: {}".format(filtered_eqtl_hg19.shape[0]))
        
        # calc tss & take abs of dist
        filtered_eqtl_hg19["tss"] = filtered_eqtl_hg19["pos"] - filtered_eqtl_hg19["dist"]
        filtered_eqtl_hg19["dist"] = filtered_eqtl_hg19["dist"].abs()
        
        # reformatting for lifting over
        filtered_eqtl_hg19["pos_end"] = filtered_eqtl_hg19["pos"] + 1
        filtered_eqtl_hg19["tss_end"] = filtered_eqtl_hg19["tss"] + 1
        
        filtered_eqtl_hg19 = filtered_eqtl_hg19[['chr', 'pos', 'pos_end', 'tss', 'tss_end', 'pairid', 'ref', 'alt', 'snpid', 'geneid', 'dist', 'npval', 'slope']]
        filtered_eqtl_hg19_eqtl = filtered_eqtl_hg19[['chr', 'pos', 'pos_end', 'pairid', 'geneid', 'npval', 'dist']]
        filtered_eqtl_hg19_eqtl.to_csv("HFB_eqtl_snp.hg19", sep="\t", index=False, header=True)
        
        filtered_eqtl_hg19_tss = filtered_eqtl_hg19[['chr', 'tss', 'tss_end', 'pairid', 'geneid', 'npval', 'dist']]
        filtered_eqtl_hg19_tss.to_csv("HFB_eqtl_tss.hg19", sep="\t", index=False, header=True)

reformat_eqtl(infile=sys.argv[1])