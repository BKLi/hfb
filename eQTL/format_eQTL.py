# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 18:29:21 2019

@author: bingkun

@project hfb -- eQTL

* Pre-process of HFB eQTL set from Geschwind Lab; output files for intersection.
* S3 of pipeline
* only need to run once 
"""


import sys
import pandas as pd
import numpy as np

#import seaborn as sns
#import matplotlib
#import matplotlib.pyplot as plt
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.rcParams['figure.dpi']= 100

def reformat_eqtl(sig_threshold, infile, random_gene, power):
        
        # sig_threshold : string; npval cutoff for sig pairs
        # infile: HFB eQTL file
        # random_gene : boolean; whether to control for exact same list of genes
        #               or only control for same number of genes
        # power: string; sampling power (1X sampling means sample same number of pairs as sig pair)
        
        ## example input:
        # random_gene = 0
        # sig_threshold = 1e-5
        # infile = r'C:\Users\libin\UCSF\hfb\eqtl\all_snp_gene_eqtl.head.50k'
        # power = 1
        
        random_gene = int(random_gene)
        sig_threshold = float(sig_threshold)
        
        
        raw_eqtl = pd.read_csv(infile, delim_whitespace=True)
        print ("number of total eqtl-gene pair: {}".format(raw_eqtl.shape[0]))
        raw_eqtl["pairid"] = ["pair_{}".format(i) for i in range(raw_eqtl.shape[0])]
        
        # filter out snp-gene pair < 10kb
        filtered_eqtl = raw_eqtl[abs(raw_eqtl["dist"]) > 10000]
        print ("number of pairs with dist > 10kb: {}".format(filtered_eqtl.shape[0]))
        
        # calc tss & take abs of dist
        filtered_eqtl["tss"] = filtered_eqtl["pos"] - filtered_eqtl["dist"]
        filtered_eqtl["dist"] = filtered_eqtl["dist"].abs()
        
        # list of all genes involved
        total_gene_list = filtered_eqtl[["geneid"]].drop_duplicates()
        print ("number of total genes involved in pairs: {}".format(total_gene_list.shape[0]))
        
        # expand eqtl region to 5kb
        filtered_eqtl["eqtl_start"] = filtered_eqtl["pos"] - 2500
        filtered_eqtl["eqtl_end"] = filtered_eqtl["pos"] + 2500
        
        # expand tss region to 5kb
        filtered_eqtl["gene_start"] = filtered_eqtl["tss"] - 2500
        filtered_eqtl["gene_end"] = filtered_eqtl["tss"] + 2500
        filtered_eqtl = filtered_eqtl[(filtered_eqtl["eqtl_start"] > 0) & (filtered_eqtl["gene_start"] > 0)]
        print ("number of pairs with start pos > 2500: {}".format(filtered_eqtl.shape[0]))
        
        #### filter for significant eqtl-gene pairs
        print ("sig_threshold: {}".format(sig_threshold))
        sig_eqtl = filtered_eqtl[filtered_eqtl["npval"] <= sig_threshold].reset_index()
        
        sig_eqtl.to_csv("HFB_sig_eqtl_tss_pair", index=False, header=True, sep="\t")
        sig_eqtl_snp_region = sig_eqtl[['chr', 'eqtl_start', 'eqtl_end', 'snpid', 'npval', 'pairid']]
        sig_eqtl_snp_region.to_csv("sig_eqtl_snp_region", index=False, header=False, sep="\t")
        sig_eqtl_gene_region = sig_eqtl[['chr', 'gene_start', 'gene_end', 'geneid', 'npval', 'pairid']]
        sig_eqtl_gene_region.to_csv("sig_eqtl_gene_region", index=False, header=False, sep="\t")

        print ("number of sig eqtl-gene pairs: {}".format(sig_eqtl.shape[0]))
        
        #### list of genes involved in sig snp-gene pairs
        sig_gene_list = sig_eqtl[["geneid"]].drop_duplicates()
        sig_gene_list.to_csv("HFB_sig_gene_list", index=False, header=True)
        print ("number of sig-genes: {}".format(sig_gene_list.shape[0]))
        
        #%% stratify distance distribution
        #print (abs(sig_eqtl["dist"].describe()))
        
        distance_bins = list(np.concatenate(
        (np.arange(10000, 80000, 10000), np.arange(80000, 200000, 20000), np.arange(200000, 300000, 50000), np.arange(300000, 600000, 150000), np.arange(600000, 1000001, 400000)
         )))
        sig_eqtl_binned = sig_eqtl.groupby(pd.cut(sig_eqtl["dist"].abs(), distance_bins))
        print ("distance distribution (sig pairs): ")
        print(sig_eqtl_binned.count()["dist"])
        
        # determine size of each bin in distance-matched sampling
        sample_size = [i*int(power) for i in sig_eqtl_binned.count()["dist"].tolist()]

        #%% filter for non-sig pairs

        nonsig_eqtl = filtered_eqtl[filtered_eqtl["npval"] > sig_threshold].reset_index()
        nonsig_eqtl_snp_region = nonsig_eqtl[['chr', 'eqtl_start', 'eqtl_end', 'snpid', 'npval', 'pairid']]
        nonsig_eqtl_snp_region.to_csv("nonsig_eqtl_snp_region", index=False, header=False, sep="\t")
        nonsig_eqtl_gene_region = nonsig_eqtl[['chr', 'gene_start', 'gene_end', 'geneid', 'npval', 'pairid']]
        nonsig_eqtl_gene_region.to_csv("nonsig_eqtl_gene_region", index=False, header=False, sep="\t")

        print ("number of non-sig eqtl-gene pairs: {}".format(nonsig_eqtl.shape[0]))
        
        
        #%% random gene sampling
        
        if random_gene:
                
                print ("controlling on rand genes...")
                #### control for same number but random list of genes
                random_gene_list = total_gene_list.sample(n=sig_gene_list.shape[0])
                random_gene_list.to_csv("HFB_random_sampled_gene_list", index=False, header=True)
                
                #### all pairs controlled on number of genes
                rand_genes_all_pair = pd.merge(filtered_eqtl, random_gene_list, on="geneid", how="inner")
                print ("number of pairs controlled on rand genes: {}".format(rand_genes_all_pair.shape[0]))
                rand_genes_all_pair.to_csv("HFB_rand_genes_all_pair", index=False, header=True, sep="\t")

                #### non-distance matched sampling
                # rand_genes_all_pair = rand_genes_all_pair.sample(n=sig_eqtl.shape[0]).reset_index()
                rand_genes_all_pair_snp_region = rand_genes_all_pair[['chr', 'eqtl_start', 'eqtl_end', 'snpid', 'npval', 'pairid']]
                rand_genes_all_pair_snp_region.to_csv("rand_genes_all_pair_snp_region", index=False, header=False, sep="\t")
                rand_genes_all_pair_gene_region = rand_genes_all_pair[['chr', 'gene_start', 'gene_end', 'geneid', 'npval', 'pairid']]
                rand_genes_all_pair_gene_region.to_csv("rand_genes_all_pair_gene_region", index=False, header=False, sep="\t")
                
                #### distance-matched sampling
                ## check distance distribution first               
                rand_genes_all_pair_binned = rand_genes_all_pair.groupby(pd.cut(rand_genes_all_pair["dist"], distance_bins))
                print ("distance distribution before sampling (all pairs): ")
                print (rand_genes_all_pair_binned.count()["dist"])
                
                rand_genes_all_pair_binned_df_list = [group.reset_index().drop(["index"], axis=1) for _, group in rand_genes_all_pair_binned]    
                rand_genes_all_pair_sampled_list = [rand_genes_all_pair_binned_df_list[i].sample(n=sample_size[i]) for i in range(len(sample_size))]                    
                rand_genes_all_pair_sampled_concat = pd.concat(rand_genes_all_pair_sampled_list).reset_index().drop(["index"], axis=1)
                rand_genes_all_pair_sampled_binned = rand_genes_all_pair_sampled_concat.groupby(pd.cut(rand_genes_all_pair_sampled_concat["dist"], distance_bins))
                print ("distance distribution after sampling (all pairs): ")
                print(rand_genes_all_pair_sampled_binned.count()["dist"])
                
                rand_genes_all_pair_sampled_concat.to_csv('HFB_rand_genes_all_pair_sampled_dist_matched', index=False, header=True, sep="\t")
                rand_genes_all_pair_sampled_concat_snp_region = rand_genes_all_pair_sampled_concat[['chr', 'eqtl_start', 'eqtl_end', 'snpid', 'npval', 'pairid']]
                rand_genes_all_pair_sampled_concat_snp_region.to_csv("rand_genes_all_pair_snp_region_dist_matched", index=False, header=False, sep="\t")
                rand_genes_all_pair_sampled_concat_gene_region = rand_genes_all_pair_sampled_concat[['chr', 'gene_start', 'gene_end', 'geneid', 'npval', 'pairid']]
                rand_genes_all_pair_sampled_concat_gene_region.to_csv("rand_genes_all_pair_gene_region_dist_matched", index=False, header=False, sep="\t")
                
                              
                #### non-sig pairs controlled on number of genes
                rand_genes_nonsig_pair = pd.merge(nonsig_eqtl, random_gene_list, on="geneid", how="inner")
                print ("number of non-sig pairs controlled on rand genes: {}".format(rand_genes_nonsig_pair.shape[0]))
                rand_genes_nonsig_pair.to_csv("HFB_rand_genes_nonsig_pair", index=False, header=True, sep="\t")
                #### non-distance matched sampling
                #rand_genes_nonsig_pair = rand_genes_nonsig_pair.sample(n=sig_eqtl.shape[0]).reset_index()
                rand_genes_nonsig_pair_snp_region = rand_genes_nonsig_pair[['chr', 'eqtl_start', 'eqtl_end', 'snpid', 'npval', 'pairid']]
                rand_genes_nonsig_pair_snp_region.to_csv("rand_genes_nonsig_pair_snp_region", index=False, header=False, sep="\t")
                rand_genes_nonsig_pair_gene_region = rand_genes_nonsig_pair[['chr', 'gene_start', 'gene_end', 'geneid', 'npval', 'pairid']]
                rand_genes_nonsig_pair_gene_region.to_csv("rand_genes_nonsig_pair_gene_region", index=False, header=False, sep="\t")

                #### distance-matched sampling
                ## check distance distribution first
                rand_genes_nonsig_pair_binned = rand_genes_nonsig_pair.groupby(pd.cut(rand_genes_nonsig_pair["dist"], distance_bins))
                print ("distance distribution before sampling (nonsig pairs): ")
                print (rand_genes_nonsig_pair_binned.count()["dist"])
                
                rand_genes_nonsig_pair_binned_df_list = [group.reset_index().drop(["index"], axis=1) for _, group in rand_genes_nonsig_pair_binned]    
                rand_genes_nonsig_pair_sampled_list = [rand_genes_nonsig_pair_binned_df_list[i].sample(n=sample_size[i]) for i in range(len(sample_size))]                    
                rand_genes_nonsig_pair_sampled_concat = pd.concat(rand_genes_nonsig_pair_sampled_list).reset_index().drop(["index"], axis=1)
                rand_genes_nonsig_pair_sampled_binned = rand_genes_nonsig_pair_sampled_concat.groupby(pd.cut(rand_genes_nonsig_pair_sampled_concat["dist"], distance_bins))
                print ("distance distribution after sampling (nonsig pairs): ")
                print(rand_genes_nonsig_pair_sampled_binned.count()["dist"])
                
                rand_genes_nonsig_pair_sampled_concat.to_csv('HFB_rand_genes_nonsig_pair_sampled_dist_matched', index=False, header=True, sep="\t")
                rand_genes_nonsig_pair_sampled_concat_snp_region = rand_genes_nonsig_pair_sampled_concat[['chr', 'eqtl_start', 'eqtl_end', 'snpid', 'npval', 'pairid']]
                rand_genes_nonsig_pair_sampled_concat_snp_region.to_csv("rand_genes_nonsig_pair_snp_region_dist_matched", index=False, header=False, sep="\t")
                rand_genes_nonsig_pair_sampled_concat_gene_region = rand_genes_nonsig_pair_sampled_concat[['chr', 'gene_start', 'gene_end', 'geneid', 'npval', 'pairid']]
                rand_genes_nonsig_pair_sampled_concat_gene_region.to_csv("rand_genes_nonsig_pair_gene_region_dist_matched", index=False, header=False, sep="\t")
                
                                   
                ### plot eqtl-tss distance distribution         
                df_dist_for_plot = pd.DataFrame()
                df_dist_for_plot["rand_gene_all"] = abs(rand_genes_all_pair_sampled_concat["dist"]).apply(float)
                df_dist_for_plot["rand_gene_nonsig"] = abs(rand_genes_nonsig_pair_sampled_concat["dist"]).apply(float)
                df_dist_for_plot["sig"] = abs(sig_eqtl["dist"]).apply(float)
                df_dist_for_plot.to_csv("HFB_rand_gene_dist_matched_for_plot", index=False, sep="\t")
                
                #df_dist_for_plot_melt = pd.melt(df_dist_for_plot, value_name="dist", var_name="")
                #sns.violinplot(x="", y="dist", data=df_dist_for_plot_melt, palette="Pastel1")
        
        #%% sig gene sampling
                       
        else:
                print ("controlling on sig genes ...")
                #### control for same genes
                #### all pairs controlled on sig-genes
                sig_genes_all_pair = pd.merge(filtered_eqtl, sig_gene_list, on="geneid", how="inner")
                print ("number of pairs controlled on sig-genes: {}".format(sig_genes_all_pair.shape[0]))
                sig_genes_all_pair.to_csv("HFB_sig_genes_all_pair", index=False, header=True, sep="\t")

                # non-distance matched sampling
                #sig_genes_all_pair = sig_genes_all_pair.sample(n=sig_eqtl.shape[0]).reset_index()
                sig_genes_all_pair_snp_region = sig_genes_all_pair[['chr', 'eqtl_start', 'eqtl_end', 'snpid', 'npval', 'pairid']]
                sig_genes_all_pair_snp_region.to_csv("sig_genes_all_pair_snp_region", index=False, header=False, sep="\t")
                sig_genes_all_pair_gene_region = sig_genes_all_pair[['chr', 'gene_start', 'gene_end', 'geneid', 'npval', 'pairid']]
                sig_genes_all_pair_gene_region.to_csv("sig_genes_all_pair_gene_region", index=False, header=False, sep="\t")
                
                #### distance-matched sampling
                ## check distance distribution first
                sig_genes_all_pair_binned = sig_genes_all_pair.groupby(pd.cut(sig_genes_all_pair["dist"], distance_bins))
                print ("distance distribution before sampling (all pairs): ")
                print (sig_genes_all_pair_binned.count()["dist"])
                
                sig_genes_all_pair_binned_df_list = [group.reset_index().drop(["index"], axis=1) for _, group in sig_genes_all_pair_binned]    
                sig_genes_all_pair_sampled_list = [sig_genes_all_pair_binned_df_list[i].sample(n=sample_size[i]) for i in range(len(sample_size))]                    
                sig_genes_all_pair_sampled_concat = pd.concat(sig_genes_all_pair_sampled_list).reset_index().drop(["index"], axis=1)
                sig_genes_all_pair_sampled_binned = sig_genes_all_pair_sampled_concat.groupby(pd.cut(sig_genes_all_pair_sampled_concat["dist"], distance_bins))
                print ("distance distribution after sampling (all pairs): ")
                print(sig_genes_all_pair_sampled_binned.count()["dist"])
                
                sig_genes_all_pair_sampled_concat.to_csv('HFB_sig_genes_all_pair_sampled_dist_matched', index=False, header=True, sep="\t")
                sig_genes_all_pair_sampled_concat_snp_region = sig_genes_all_pair_sampled_concat[['chr', 'eqtl_start', 'eqtl_end', 'snpid', 'npval', 'pairid']]
                sig_genes_all_pair_sampled_concat_snp_region.to_csv("sig_genes_all_pair_snp_region_dist_matched", index=False, header=False, sep="\t")
                sig_genes_all_pair_sampled_concat_gene_region = sig_genes_all_pair_sampled_concat[['chr', 'gene_start', 'gene_end', 'geneid', 'npval', 'pairid']]
                sig_genes_all_pair_sampled_concat_gene_region.to_csv("sig_genes_all_pair_gene_region_dist_matched", index=False, header=False, sep="\t")
                
                
                #### non-sig pairs controlled on sig-genes
                sig_genes_nonsig_pair = pd.merge(nonsig_eqtl, sig_gene_list, on="geneid", how="inner")
                print ("number of non-sig pairs controlled on sig-genes: {}".format(sig_genes_nonsig_pair.shape[0]))
                sig_genes_nonsig_pair.to_csv('HFB_sig_genes_nonsig_pair', index=False, header=True, sep="\t")                               
                # non-distance matched sampling
                #sig_genes_nonsig_pair = sig_genes_nonsig_pair.sample(n=sig_eqtl.shape[0]).reset_index()
                sig_genes_nonsig_pair_snp_region = sig_genes_nonsig_pair[['chr', 'eqtl_start', 'eqtl_end', 'snpid', 'npval', 'pairid']]
                sig_genes_nonsig_pair_snp_region.to_csv("sig_genes_nonsig_pair_snp_region", index=False, header=False, sep="\t")
                sig_genes_nonsig_pair_gene_region = sig_genes_nonsig_pair[['chr', 'gene_start', 'gene_end', 'geneid', 'npval', 'pairid']]
                sig_genes_nonsig_pair_gene_region.to_csv("sig_genes_nonsig_pair_gene_region", index=False, header=False, sep="\t")
                
                #### distance-matched sampling
                ## check distance distribution first
                sig_genes_nonsig_pair_binned = sig_genes_nonsig_pair.groupby(pd.cut(sig_genes_nonsig_pair["dist"], distance_bins))
                print ("distance distribution before sampling (nonsig pairs): ")
                print (sig_genes_nonsig_pair_binned.count()["dist"])
                
                sig_genes_nonsig_pair_binned_df_list = [group.reset_index().drop(["index"], axis=1) for _, group in sig_genes_nonsig_pair_binned]    
                sig_genes_nonsig_pair_sampled_list = [sig_genes_nonsig_pair_binned_df_list[i].sample(n=sample_size[i]) for i in range(len(sample_size))]                    
                sig_genes_nonsig_pair_sampled_concat = pd.concat(sig_genes_nonsig_pair_sampled_list).reset_index().drop(["index"], axis=1)
                sig_genes_nonsig_pair_sampled_binned = sig_genes_nonsig_pair_sampled_concat.groupby(pd.cut(sig_genes_nonsig_pair_sampled_concat["dist"], distance_bins))
                print ("distance distribution after sampling (nonsig pairs): ")
                print(sig_genes_nonsig_pair_sampled_binned.count()["dist"])
                
                sig_genes_nonsig_pair_sampled_concat.to_csv('HFB_sig_genes_nonsig_pair_sampled_dist_matched', index=False, header=True, sep="\t")
                sig_genes_nonsig_pair_sampled_concat_snp_region = sig_genes_nonsig_pair_sampled_concat[['chr', 'eqtl_start', 'eqtl_end', 'snpid', 'npval', 'pairid']]
                sig_genes_nonsig_pair_sampled_concat_snp_region.to_csv("sig_genes_nonsig_pair_snp_region_dist_matched", index=False, header=False, sep="\t")
                sig_genes_nonsig_pair_sampled_concat_gene_region = sig_genes_nonsig_pair_sampled_concat[['chr', 'gene_start', 'gene_end', 'geneid', 'npval', 'pairid']]
                sig_genes_nonsig_pair_sampled_concat_gene_region.to_csv("sig_genes_nonsig_pair_gene_region_dist_matched", index=False, header=False, sep="\t")
                
                ### plot eqtl-tss distance distributione
                df_dist_for_plot = pd.DataFrame()
                df_dist_for_plot["sig_gene_all"] = abs(sig_genes_all_pair_sampled_concat["dist"]).apply(float)
                df_dist_for_plot["sig_gene_nonsig"] = abs(sig_genes_nonsig_pair_sampled_concat["dist"]).apply(float)
                df_dist_for_plot["sig"] = abs(sig_eqtl["dist"]).apply(float)
                df_dist_for_plot.to_csv("HFB_sig_gene_dist_matched_for_plot", index=False, sep="\t")

        
reformat_eqtl(sig_threshold=sys.argv[1], infile=sys.argv[2], random_gene=sys.argv[3], power=sys.argv[4])