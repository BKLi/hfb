# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 18:16:41 2020

@author: bingkun
@project: hfb -- H-MAGMA

* process of MAGMA output
* s3 of pipeline
"""

import pandas as pd
import itertools

cellType_list=["IN", "IPC", "RGC", "neuron"]
disease_list= ["ASD", "AD", "ADHD", "BIP", "IQ", "MDD", "SCZ"]
cutoff = 0.05

geneID2geneName = pd.read_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.geneID2geneNameType', sep="\t")

for x, y in list(itertools.product(cellType_list, disease_list)):
        cellType = x
        disease = y
        print (cellType, disease)

        #%% magma output using FB annotation from won lab
        #magma_out_fb = pd.read_csv(r'C:\Users\libin\UCSF\MAGMA\output\AD_FB.genes.out', delim_whitespace=True)
        #magma_out_fb = magma_out_fb.rename(columns={"GENE":"gene_id", "P":"FB_P"})
        #magma_out_fb_withName = pd.merge(geneID2geneName, magma_out_fb, on=["gene_id"], how="inner")
        #magma_out_fb_withName = magma_out_fb_withName[magma_out_fb_withName["gene_type"] == "protein_coding"]
        #magma_out_fb_withName_sig = magma_out_fb_withName[magma_out_fb_withName["FB_P"] < 0.05]
        #magma_out_fb_withName_sig_nameList = magma_out_fb_withName_sig[["gene_name"]]
        #magma_out_fb_withName_sig_nameList.to_csv(r'C:\Users\libin\UCSF\MAGMA\output\magma_out_fb_{}_geneList'.format(disease), header=False, index=False)
        
        #%% magma output using CT-specific annotation
        magma_out_ct = pd.read_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\output\{}_{}_all.genes.out'.format(disease, cellType),
                                delim_whitespace=True)
        magma_out_ct = magma_out_ct.rename(columns={"GENE":"gene_id", "P":"{}_P".format(cellType)})
        magma_out_ct_withName = pd.merge(geneID2geneName, magma_out_ct, on=["gene_id"], how="inner")
        magma_out_ct_withName = magma_out_ct_withName[magma_out_ct_withName["gene_type"] == "protein_coding"]
        magma_out_ct_withName_sig = magma_out_ct_withName[magma_out_ct_withName["{}_P".format(cellType)] < 0.05]
        magma_out_ct_withName_sig_nameList = magma_out_ct_withName_sig[["gene_name"]]
        magma_out_ct_withName_sig_nameList.to_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\output\magma_out_{}_{}_geneList'.format(cellType, disease), header=False, index=False)
        
        
        #%% magma output using target-only CT-specific annotation
        magma_out_ct = pd.read_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\output\{}_{}_target.genes.out'.format(disease, cellType),
                                delim_whitespace=True)
        magma_out_ct = magma_out_ct.rename(columns={"GENE":"gene_id", "P":"{}_P".format(cellType)})
        magma_out_ct_withName = pd.merge(geneID2geneName, magma_out_ct, on=["gene_id"], how="inner")
        magma_out_ct_withName = magma_out_ct_withName[magma_out_ct_withName["gene_type"] == "protein_coding"]
        magma_out_ct_withName_sig = magma_out_ct_withName[magma_out_ct_withName["{}_P".format(cellType)] < 0.05]
        magma_out_ct_withName_sig_nameList = magma_out_ct_withName_sig[["gene_name"]]
        magma_out_ct_withName_sig_nameList.to_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\output\magma_out_{}_{}_target_geneList'.format(cellType, disease), header=False, index=False)
        
        #%% compare
        #magma_out_fb_withName_cut = magma_out_fb_withName[['gene_id','gene_name','gene_type','FB_P']]
        #magma_out_ct_withName_cut = magma_out_ct_withName[['gene_id','gene_name','gene_type',"{}_P".format(cellType)]]
        
        #magma_compare = pd.merge(magma_out_ct_withName_cut, magma_out_fb_withName_cut, on=['gene_id','gene_name','gene_type'], how="inner")