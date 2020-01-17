# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:01:10 2020

@author: bingkun

@project hfb -- eQTL

Script for processing complete (pospoisson) post-liftover interations

* S2 of pipeline
* collect boths ends of interactions after liftover
* data cleaning
* interactions sampling
* output file for intersection

!!!!! deprecated -- lifting over eQTL to hg38 instead for simplicity 
"""

import pandas as pd

cellType='neuron'

#%%
################################## PART ONE #####################################

# read in interactions before liftover, with ID
raw_all = pd.read_csv(r'C:\Users\libin\UCSF\hfb\eqtl\hfb\MAPS_allPeaks\{}.xor.MAPS2_pospoisson.reform.withID'.format(cellType), sep="\t",
                      names=['ID','bin1_mid','bin2_mid','count','X1D_peak_bin1','X1D_peak_bin2','effective_length1',
                             'gc1','mappability1','short_count1','effective_length2','gc2','mappability2',
                             'short_count2','dist','logl','loggc','logm','logdist','logShortCount','chr','expected',
                             'p_val','expected2','ratio2','p_val_reg2','p_bonferroni','fdr'])
print ("number of raw interactions: {}".format(raw_all.shape[0]))
#%%
# read in both ends of interactions after lifting over
end1 = pd.read_csv(r'C:\Users\libin\UCSF\hfb\eqtl\hfb\MAPS_allPeaks\{}.xor.MAPS2_pospoisson.reform.withID.bin1.hg19'.format(cellType), sep="\t")
end2 = pd.read_csv(r'C:\Users\libin\UCSF\hfb\eqtl\hfb\MAPS_allPeaks\{}.xor.MAPS2_pospoisson.reform.withID.bin2.hg19'.format(cellType), sep="\t")
#%%
end_both = pd.merge(end1, end2, how="inner", on=['ID', 'bin1', 'bin2', 'fdr'])
end_both = end_both[['chr1', 'start1', 'chr2', 'start2', 'ID', 'fdr', 'bin1', 'bin2']]
#%%
# remove interactions whose ends got mapped to diff chr after liftover
misplaced_index = end_both[end_both["chr1"] != end_both["chr2"]].index
end_both = end_both.drop(misplaced_index)
#%%
end_both["end1"] = end_both["start1"] + 5000
end_both["end2"] = end_both["start2"] + 5000
print("number of interactons after liftover: {}".format(end_both.shape[0]))

#%% 
############################ PART TWO -- LABEL PROMOTERS ############################
end_both_bin1_anchor = end_both[end_both["bin1"] == 1]
end_both_bin2_anchor = end_both[end_both["bin2"] == 1]
#%%
# rename so start1 is always anchor bin
end_both_bin2_anchor = end_both_bin2_anchor.rename(columns={"chr1":"chr2", "chr2":"chr1", "start1":"start2", "start2":"start1", "end1":"end2", "end2":"end1"})
#%%
end_both_bin2_anchor = end_both_bin2_anchor[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'ID', 'fdr', 'bin1', 'bin2']]
end_both_bin1_anchor = end_both_bin1_anchor[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'ID', 'fdr', 'bin1', 'bin2']]
end_both_concat = pd.concat([end_both_bin1_anchor, end_both_bin2_anchor])
#%%
# sampling based on fdr
print(end_both_concat["fdr"].describe())




