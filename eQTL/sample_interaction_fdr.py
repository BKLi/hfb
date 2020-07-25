# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 16:02:37 2020

@author: libin
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 13:43:10 2020

@author: bingkun
@project: hfb -- eQTL

sample nonsig interactions from MAPS pospoisson file
take first and last n interactions (sorted by fdr)

* S2 of alternative pipeline
"""

#%%
import pandas as pd
import numpy as np

#%%
####### read in complete MAPS interaction sets
cellType='IPC'
#sample_power
power = 1

inter_all = pd.read_csv(r'C:\Users\libin\UCSF\hfb\eqtl\hfb\MAPS_allPeaks\{}.xor.MAPS2_pospoisson.reform.withID'.format(cellType), sep="\t",
                      names=['ID','bin1_mid','bin2_mid','count','X1D_peak_bin1','X1D_peak_bin2','effective_length1',
                             'gc1','mappability1','short_count1','effective_length2','gc2','mappability2',
                             'short_count2','dist','logl','loggc','logm','logdist','logShortCount','chr','expected',
                             'p_val','expected2','ratio2','p_val_reg2','p_bonferroni','fdr'])
inter_all = inter_all[['ID', 'chr', 'bin1_mid', 'bin2_mid', 'X1D_peak_bin1', 'X1D_peak_bin2', 'fdr', "dist"]]
inter_all["dist"] = inter_all["dist"]*5000
print ("number of all XOR interactions: {}".format(inter_all.shape[0]))

inter_all_bin1_anchor = inter_all[inter_all["X1D_peak_bin1"] == 1]
inter_all_bin2_anchor = inter_all[inter_all["X1D_peak_bin2"] == 1]
# rename so bin1 is always anchor bin
inter_all_bin2_anchor = inter_all_bin2_anchor.rename(columns={"bin1_mid":"bin2_mid", "bin2_mid":"bin1_mid"})
inter_all_bin2_anchor = inter_all_bin2_anchor[['ID', 'chr', 'bin1_mid', 'bin2_mid', 'X1D_peak_bin1', 'X1D_peak_bin2', 'fdr', "dist"]]

inter_all_reformed = pd.concat([inter_all_bin1_anchor, inter_all_bin2_anchor])
inter_all_reformed['bin1_mid'] = inter_all_reformed['bin1_mid'].astype('int')
inter_all_reformed['bin2_mid'] = inter_all_reformed['bin2_mid'].astype('int')

inter_all_reformed['bin1_end'] = inter_all_reformed['bin1_mid'] + 5000
inter_all_reformed['bin2_end'] = inter_all_reformed['bin2_mid'] + 5000

inter_all_reformed["chr1"] = inter_all_reformed["chr"]
inter_all_reformed["chr2"] = inter_all_reformed["chr"]

inter_all_reformed = inter_all_reformed.rename(columns={"bin1_mid":"start1", "bin2_mid":"start2", "bin1_end":"end1", "bin2_end":"end2"})
inter_all_reformed = inter_all_reformed[['ID','chr','start1','start2','X1D_peak_bin1','X1D_peak_bin2','fdr','chr1','chr2','end1','end2',"dist"]]
inter_all_reformed = inter_all_reformed[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'fdr', 'ID', 'X1D_peak_bin1', 'X1D_peak_bin2', 'dist']]

# sort by fdr
inter_all_reformed = inter_all_reformed.sort_values(by='fdr')

#%%
### ---------------------------------------------------------------------------------------------------------------
distance_bins = list(np.concatenate(
        (np.arange(0, 20000, 10000), np.arange(20000, 80000, 30000), np.arange(80000, 200000, 20000), 
         np.arange(200000, 300000, 25000), np.arange(300000, 500000, 50000), 
         np.arange(500000, 700000, 100000),np.arange(700000, 1000001, 150000)
         )))

####### sample most significant and non-sig MAPS interaction sets
inter_sig = inter_all_reformed.head(n=10000)

inter_sig_binned = inter_sig.groupby(pd.cut(inter_sig["dist"].abs(), distance_bins))
print ("distance distribution (sig interaction): ")
print(inter_sig_binned.count()["dist"])

inter_nonsig = inter_all_reformed.tail(n=10000)

inter_nonsig_binned = inter_nonsig.groupby(pd.cut(inter_nonsig["dist"].abs(), distance_bins))
print ("distance distribution (non-sig interaction): ")
print(inter_nonsig_binned.count()["dist"])


#plt.figure(figsize=(10,7))
#sns.distplot(inter_sig_XOR_reformed["dist_sig"], hist=False, color="blue", label="sig")
#sns.distplot(inter_all_reformed_nonsig_sampled_concat["dist"], hist=False, color="orange", label="nonsig")
#plt.xlabel('distance distribution ({})'.format(cellType), fontsize=15)
#plt.savefig(r'C:\Users\libin\UCSF\hfb\eqtl\20200126\{}_dist_compare.pdf'.format(cellType), transparent=True)

#%%

inter_sig_anchor = inter_sig[['chr1','start1','end1','fdr','ID', "dist"]]
inter_sig_anchor.to_csv(r'C:\Users\libin\UCSF\hfb\eqtl\20200126\interaction_sampled_sig_10k_anchor_{}'.format(cellType), sep="\t", index=False, header=False)
inter_sig_target = inter_sig[['chr2','start2','end2','fdr','ID', "dist"]]
inter_sig_target.to_csv(r'C:\Users\libin\UCSF\hfb\eqtl\20200126\interaction_sampled_sig_10k_target_{}'.format(cellType), sep="\t", index=False, header=False)


inter_nonsig_anchor = inter_nonsig[['chr1','start1','end1','fdr','ID', "dist"]]
inter_nonsig_anchor.to_csv(r'C:\Users\libin\UCSF\hfb\eqtl\20200126\interaction_sampled_nonsig_10k_anchor_{}'.format(cellType), sep="\t", index=False, header=False)
inter_nonsig_target = inter_nonsig[['chr2','start2','end2','fdr','ID', "dist"]]
inter_nonsig_target.to_csv(r'C:\Users\libin\UCSF\hfb\eqtl\20200126\interaction_sampled_nonsig_10k_target_{}'.format(cellType), sep="\t", index=False, header=False)











