# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 18:16:41 2020

@author: bingkun
@project: hfb -- H-MAGMA

* process of MAGMA output
* s3 of pipeline
"""

import pandas as pd

cellType="neuron"
disease= "AD"

geneID2geneName = pd.read_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.geneID2geneNameType', sep="\t")

#%% magma output using FB annotation from won lab
magma_out_fb = pd.read_csv(r'C:\Users\libin\UCSF\MAGMA\output\AD_FB.genes.out', delim_whitespace=True)
magma_out_fb = magma_out_fb.rename(columns={"GENE":"gene_id", "P":"FB_P"})
magma_out_fb_withName = pd.merge(geneID2geneName, magma_out_fb, on=["gene_id"], how="inner")
magma_out_fb_withName = magma_out_fb_withName[magma_out_fb_withName["gene_type"] == "protein_coding"]

#%% magma output using CT-specific annotation
magma_out_ct = pd.read_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\output\{}_{}.genes.out'.format(disease, cellType),
                        delim_whitespace=True)
magma_out_ct = magma_out_ct.rename(columns={"GENE":"gene_id", "P":"{}_P".format(cellType)})
magma_out_ct_withName = pd.merge(geneID2geneName, magma_out_ct, on=["gene_id"], how="inner")
magma_out_ct_withName = magma_out_ct_withName[magma_out_ct_withName["gene_type"] == "protein_coding"]

#%% compare
magma_out_fb_withName_cut = magma_out_fb_withName[['gene_id','gene_name','gene_type','FB_P']]
magma_out_ct_withName_cut = magma_out_ct_withName[['gene_id','gene_name','gene_type',"{}_P".format(cellType)]]

magma_compare = pd.merge(magma_out_ct_withName_cut, magma_out_fb_withName_cut, on=['gene_id','gene_name','gene_type'], how="inner")