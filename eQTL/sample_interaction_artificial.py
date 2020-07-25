# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 11:55:40 2020

@author: libin
@project: hfb -- eQTL

create artificial interactions with sampled distances from MAPS pospoisson file

only sample distances from interactions with the same anchors as sig-interactions

including both AND and XOR interactions

* not really part of the pipeline
"""

#%%
import random
from random import seed
from random import randint
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['figure.dpi']= 200


cellType='IPC'
#sample_power
# power = 1
iteration = 5
date = 20200218
brk_num = 4

#%%
### read in chrom sizes
chrom_size = pd.read_csv(r'C:\Users\libin\UCSF\hg38_general\hg38.chrom.sizes.noContig', sep="\t", names=["chr1", "size"])
chrom_size.loc[:, "chr_start"] = 1
chrom_size.loc[:, "chr_end"] = chrom_size["size"]
chrom_size = chrom_size[["chr1", "chr_start", "chr_end"]]

####### read in & process significant MAPS interaction sets
inter_sig = pd.read_csv(r'C:\Users\libin\UCSF\hfb\interactions_final\{}.5k.downsample.bedpe'.format(cellType), sep="\t")
inter_sig["sig_ID"] = ["sig_{}_{}".format(cellType, i) for i in range(inter_sig.shape[0])]
print ("number of signif interactions: {}".format(inter_sig.shape[0]))

anchor_bins = pd.read_csv(r'C:\Users\libin\UCSF\hfb\interactions_final\MACS2_peak_5KbBin_Union4.txt', sep="\t")
anchor_bins["anchor_ID"] = ["anchor_{}".format(i) for i in range(anchor_bins.shape[0])]

inter_sig_left_anchor = pd.merge(inter_sig, anchor_bins, left_on=['chr1', 'start1', 'end1'],
                                    right_on = ['chr', 'start', 'end'], how="inner")
left_anchors_list = inter_sig_left_anchor[['sig_ID', 'anchor_ID']]

inter_sig_right_anchor = pd.merge(inter_sig, anchor_bins, left_on=['chr2', 'start2', 'end2'],
                                     right_on = ['chr', 'start', 'end'], how="inner")
right_anchors_list = inter_sig_right_anchor[['sig_ID', 'anchor_ID']]

anchor_anchor_list = pd.merge(left_anchors_list, right_anchors_list, on="sig_ID", how="inner")

inter_sig_XOR = inter_sig[~inter_sig["sig_ID"].isin(anchor_anchor_list["sig_ID"])]
print ("number of signif XOR interactions: {}".format(inter_sig_XOR.shape[0]))

inter_sig_AND = inter_sig[inter_sig["sig_ID"].isin(anchor_anchor_list["sig_ID"])]
inter_sig_AND.loc[:, "anchor"] = "both"

print ("number of signif AND interactions: {}".format(inter_sig_AND.shape[0]))


inter_sig_left_anchor_XOR = pd.merge(inter_sig_XOR, anchor_bins, left_on=['chr1', 'start1', 'end1'],
                                    right_on = ['chr', 'start', 'end'], how="inner")
inter_sig_left_anchor_XOR = inter_sig_left_anchor_XOR[['chr1','start1','end1','chr2','start2','end2','count','expected','fdr','ClusterLabel','ClusterSize',
                                                         'ClusterType','ClusterNegLog10P','ClusterSummit','sig_ID']]
inter_sig_left_anchor_XOR.loc[:, "anchor"] = "left"
print ("number of signif interactions with bin1 as anchor bin: {}".format(inter_sig_left_anchor_XOR.shape[0]))

inter_sig_right_anchor_XOR = pd.merge(inter_sig_XOR, anchor_bins, left_on=['chr2', 'start2', 'end2'],
                                     right_on = ['chr', 'start', 'end'], how="inner")
print ("number of signif interactions with bin2 as anchor bin: {}".format(inter_sig_right_anchor_XOR.shape[0]))
# rename so bin1(start1) is always anchor bin
inter_sig_right_anchor_XOR = inter_sig_right_anchor_XOR.rename(columns={"chr1":"chr2", "start1":"start2", "end1":"end2", "chr2":"chr1", "start2":"start1", "end2":"end1"})
inter_sig_right_anchor_XOR = inter_sig_right_anchor_XOR[['chr1','start1','end1','chr2','start2','end2','count','expected','fdr','ClusterLabel','ClusterSize',
                                                         'ClusterType','ClusterNegLog10P','ClusterSummit','sig_ID']]
inter_sig_right_anchor_XOR.loc[:, "anchor"] = "right"

inter_sig_reformed = pd.concat([inter_sig_AND, inter_sig_left_anchor_XOR, inter_sig_right_anchor_XOR])
inter_sig_reformed = inter_sig_reformed[['chr1','start1','end1','chr2','start2','end2','count','expected','fdr','ClusterLabel','ClusterSize',
                                                         'ClusterType','ClusterNegLog10P','ClusterSummit','sig_ID', "anchor"]]

inter_sig_reformed.loc[:, "dist_sig"] = inter_sig_reformed["start2"] - inter_sig_reformed["start1"]
#sort from low to high
inter_sig_reformed = inter_sig_reformed.sort_values(by=['chr1', 'start1'], ascending=True)

# extract list of anchors (keep dups)
# inter_sig_anchor_list = inter_sig_reformed[['chr1','start1','end1']]

### group by chromosome
inter_sig_reformed_chr = inter_sig_reformed.groupby(['chr1', 'chr2'])
inter_sig_reformed_chr_df_list = [chrs for _, chrs in inter_sig_reformed_chr]

# function for shuffling distances 
def shuf_dist(inter_sig, brk_num, random_num): 
        
        inter_sig = inter_sig.reset_index().drop(['index'], axis=1)
        
        break1 = int(inter_sig["start1"].quantile([.25, .5, .75]).tolist()[0])
        break2 = int(inter_sig["start1"].quantile([.25, .5, .75]).tolist()[1])
        break3 = int(inter_sig["start1"].quantile([.25, .5, .75]).tolist()[2])
        # print (break1, break2, break3)
        
        inter_sig_1 = inter_sig[inter_sig["start1"] <= break1]
        #print (inter_sig_1["dist_sig"].min())
        inter_sig_2 = inter_sig[(inter_sig["start1"] > break1) & (inter_sig["start1"] <= break2)]
        #print (inter_sig_2["dist_sig"].min())

        inter_sig_3 = inter_sig[(inter_sig["start1"] > break2) & (inter_sig["start1"] <= break3)]
        #print (inter_sig_3["dist_sig"].min())

        inter_sig_4 = inter_sig[inter_sig["start1"] > break3]
        #print (inter_sig_4["dist_sig"].min())


        print ("block size: ", inter_sig_1.shape[0], inter_sig_2.shape[0], inter_sig_3.shape[0], inter_sig_4.shape[0])
        
        dist_list = inter_sig["dist_sig"]
        #print ("minimal all: ", dist_list.min())
        # sort from high to low
        dist_list_sorted = sorted(dist_list, reverse=True)
        #print ('max and min of sorted: ', dist_list_sorted[0], dist_list_sorted[-1])
        #print (dist_list_sorted[0:10])
        #print (dist_list_sorted[-10:])

        dist_list_sorted_1 = dist_list_sorted[0:inter_sig_1.shape[0]]
        #print (min(dist_list_sorted_1))
        #print (dist_list_sorted_1[0:10])
        dist_list_sorted_2 = dist_list_sorted[inter_sig_1.shape[0]: inter_sig_1.shape[0] + inter_sig_2.shape[0]]
        #print (min(dist_list_sorted_2))

        dist_list_sorted_3 = dist_list_sorted[inter_sig_1.shape[0] + inter_sig_2.shape[0]: inter_sig_1.shape[0] + inter_sig_2.shape[0] + inter_sig_3.shape[0]]
        #print (min(dist_list_sorted_3))

        dist_list_sorted_4 = dist_list_sorted[inter_sig_1.shape[0] + inter_sig_2.shape[0] + inter_sig_3.shape[0]: inter_sig_1.shape[0] + inter_sig_2.shape[0] + inter_sig_3.shape[0] + inter_sig_4.shape[0]]
        #print (min(dist_list_sorted_4))
        #print (dist_list_sorted_4[-10:])

        # print(len(dist_list_sorted_1), len(dist_list_sorted_2), len(dist_list_sorted_3), len(dist_list_sorted_4))
        
        dist_list_sorted_1_shuf = random.sample(dist_list_sorted_1, len(dist_list_sorted_1)) 
        #print (min(dist_list_sorted_1_shuf))
        dist_list_sorted_2_shuf = random.sample(dist_list_sorted_2, len(dist_list_sorted_2)) 
        #print (min(dist_list_sorted_2_shuf))
        dist_list_sorted_3_shuf = random.sample(dist_list_sorted_3, len(dist_list_sorted_3)) 
        #print (min(dist_list_sorted_2_shuf))
        dist_list_sorted_4_shuf = random.sample(dist_list_sorted_4, len(dist_list_sorted_4)) 
        #print (min(dist_list_sorted_4_shuf))

        #print (dist_list_sorted_1_shuf[0:10])
        inter_shuf_1 = inter_sig_1.loc[:, ['chr1','start1','end1','chr2','start2','end2', 'dist_sig']].reset_index().drop(["index"], axis=1)
        inter_shuf_1.loc[:, "dist_shuf"] = dist_list_sorted_1_shuf
        inter_shuf_1.loc[:, "start2_shuf"] = inter_shuf_1["start1"] + inter_shuf_1["dist_shuf"]
        inter_shuf_1.loc[:, "end2_shuf"] = inter_shuf_1["start2_shuf"] + 5000
        
        inter_shuf_2 = inter_sig_2.loc[:, ['chr1','start1','end1','chr2','start2','end2', 'dist_sig']].reset_index().drop(["index"], axis=1)
        inter_shuf_2.loc[:, "dist_shuf"] = dist_list_sorted_2_shuf
        inter_shuf_2.loc[:, "start2_shuf"] = inter_shuf_2["start1"] + inter_shuf_2["dist_shuf"]
        inter_shuf_2.loc[:, "end2_shuf"] = inter_shuf_2["start2_shuf"] + 5000
        
        inter_shuf_3 = inter_sig_3.loc[:, ['chr1','start1','end1','chr2','start2','end2', 'dist_sig']].reset_index().drop(["index"], axis=1)
        inter_shuf_3.loc[:, "dist_shuf"] = dist_list_sorted_3_shuf
        inter_shuf_3.loc[:, "start2_shuf"] = inter_shuf_3["start1"] + inter_shuf_3["dist_shuf"]
        inter_shuf_3.loc[:, "end2_shuf"] = inter_shuf_3["start2_shuf"] + 5000
        
        inter_shuf_4 = inter_sig_4.loc[:, ['chr1','start1','end1','chr2','start2','end2', 'dist_sig']].reset_index().drop(["index"], axis=1)
        inter_shuf_4.loc[:, "dist_shuf"] = dist_list_sorted_4_shuf
        inter_shuf_4.loc[:, "start2_shuf"] = inter_shuf_4["start1"] + inter_shuf_4["dist_shuf"]
        inter_shuf_4.loc[:, "end2_shuf"] = inter_shuf_4["start2_shuf"] + 5000
        
        inter_shuf_all = pd.concat([inter_shuf_1, inter_shuf_2, inter_shuf_3, inter_shuf_4])
        # print (inter_shuf_all.head())
        # print (list(inter_shuf_all))
        inter_shuf_wSize = pd.merge(inter_shuf_all, chrom_size, on="chr1", how="inner")
        inter_shuf_wSize.loc[:, "chr_end_diff"] = inter_shuf_wSize["chr_end"] - inter_shuf_wSize["end2_shuf"]
        # check if start position > 0 
        inter_shuf_wSize_minStart = inter_shuf_wSize['start2_shuf'].min()
        # check if end position < chrom size
        inter_shuf_wSize_minDiff = inter_shuf_wSize['chr_end_diff'].min()
        return(inter_shuf_all, inter_shuf_wSize, inter_shuf_wSize_minStart, inter_shuf_wSize_minDiff)
        

chr_df_list = []
for chr_df in inter_sig_reformed_chr_df_list:
        
        #print (chr_df["])
        chr_num = chr_df["chr1"].tolist()[0]
        print (chr_num)
        i = iteration
        seed(i)
        ran_num = randint(1, 1000000)
        print (ran_num)
        inter_shuf_out, inter_shuf_wSize_out, inter_shuf_wSize_minStart_out, inter_shuf_wSize_minDiff_out = shuf_dist(chr_df, brk_num, ran_num)
        #shuf_dist(chr_df, brk_num, ran_num)
        print (inter_shuf_wSize_minStart_out, inter_shuf_wSize_minDiff_out)
        
        if inter_shuf_wSize_minStart_out > 0 and inter_shuf_wSize_minDiff_out > 0:
                print ("successful at run1")
                plt.figure(figsize=(7,4))
                sns.distplot(inter_shuf_wSize_out["dist_sig"], hist=False, color="blue", label="sig")
                sns.distplot(inter_shuf_wSize_out["dist_shuf"], hist=False, color="orange", label="nonsig")
                plt.xlabel('distance distribution ({})'.format(chr_num), fontsize=15)
                #plt.savefig(r'C:\Users\libin\UCSF\hfb\eqtl\20200218\{}_dist_compare_run{}_{}.pdf'.format(cellType, iteration, chr_num), transparent=True)

                chr_df_list.append(inter_shuf_wSize_out)       
                    
        else:
                print ("starting iterations...")
                while inter_shuf_wSize_minStart_out < 0 or inter_shuf_wSize_minDiff_out < 0:
                        i += 1
                        seed(i)
                        ran_num = randint(1, 1000000)
                        inter_shuf_out, inter_shuf_wSize_out, inter_shuf_wSize_minStart_out, inter_shuf_wSize_minDiff_out =  shuf_dist(chr_df, brk_num, ran_num)                          
                        print (ran_num)
                        print ("iteration: ", i)
                        print (inter_shuf_wSize_minStart_out, inter_shuf_wSize_minDiff_out)
                plt.figure(figsize=(7,4))
                sns.distplot(inter_shuf_wSize_out["dist_sig"], hist=False, color="blue", label="sig")
                sns.distplot(inter_shuf_wSize_out["dist_shuf"], hist=False, color="orange", label="nonsig")
                plt.xlabel('distance distribution ({})'.format(chr_df["chr1"].tolist()[0]), fontsize=15)
                #plt.savefig(r'C:\Users\libin\UCSF\hfb\eqtl\20200218\{}_dist_compare_run{}_{}.pdf'.format(cellType, iteration, chr_num), transparent=True)

                chr_df_list.append(inter_shuf_wSize_out)

inter_shuf_final = pd.concat(chr_df_list)
inter_shuf_final = inter_shuf_final[['chr1','start1','end1','chr2','start2_shuf','end2_shuf', 'dist_shuf', 'chr_end_diff', 'start2', 'end2', 'dist_sig']]
inter_shuf_final.to_csv(r'C:\Users\libin\UCSF\hfb\eqtl\20200218\interaction_artificial_{}_run{}'.format(cellType, iteration), sep="\t", index=False, header=True)

intersect_shuf_raw = pd.merge(inter_shuf_final, inter_sig, left_on=['chr1','start1','end1','chr2','start2_shuf','end2_shuf'], 
                              right_on=['chr1','start1','end1','chr2','start2','end2'], how="inner")
print ("number of overlap: ", intersect_shuf_raw.shape[0])


#### -------------- incomplete test code below ------------
'''
        break_list = list(np.arange(0,1,1/brk_num))
        print (break_list)
        
        break_point_list = []
        for b in break_list:
                break_i = inter_sig["start1"].quantile(b)
                break_point_list.append(break_i)
        print (break_point_list)
        break_df_list = []
        for i in range(len(break_point_list)):
                intersig_break = 
                
                #print (break_i)
                #inter_sig_brk = inter_sig[inter_sig["start1"] <= break_i]
        # determine breaking points of each stratification       
        '''












