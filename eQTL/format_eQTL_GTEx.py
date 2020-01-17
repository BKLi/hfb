# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 13:39:49 2020

@author: bingkun
@project hfb -- eQTL

* Pre-processing of GTEx eQTL data
* Using signif sets provided by GTEx
* sig-gene defined as genes appeared in signif pairs sets
* Generate control sets
* Output file ready for intersecting
"""


import sys
import pandas as pd
import numpy as np

def reform_eQTL(data_dir, tissue, process_raw_pair, HFB_sig_gene_file, power):
        # data_dir: path for eQTL data
        # tissue: tissue name
        # process_raw_pair: boolean. Whether or not to process complete dataset 
        # HFB_sig_gene_file: list of HFB eQTL sig genes
        # power: string; sampling power (1X sampling means sample same number of pairs as sig pair)

        ## example input
        #data_dir = r'C:\Users\libin\UCSF\hfb\eqtl\\'
        #tissue = "Brain_Cortex"
        #HFB_sig_gene_file = r'C:\Users\libin\UCSF\hfb\eqtl\HFB_sig_gene_list'
        #power = 1
        
        sig_file = r'{}{}.v7.signif_variant_gene_pairs.txt'.format(data_dir, tissue)
        all_file = r'{}{}.allpairs.txt'.format(data_dir, tissue)
               
        ##
        HFB_sig_gene = pd.read_csv(HFB_sig_gene_file)
        ## ------------------------------- signif data --------------------------------
        eQTL_raw_sig = pd.read_table(sig_file, delim_whitespace=True)
            
        # keep only eQTL-TSS pairs with a tss_distance of at least 10kb
        eQTL_raw_sig = eQTL_raw_sig[eQTL_raw_sig["tss_distance"].abs() > 10000]
        print ("number of signif eQTL-gene pairs (dist > 10kb) : {}".format(eQTL_raw_sig.shape[0]))
        
        eQTL_raw_sig["chr"] = "chr" + eQTL_raw_sig["variant_id"].str.split("_", expand=True)[0]
        eQTL_raw_sig["pos_eqtl"] = eQTL_raw_sig["variant_id"].str.split("_", expand=True)[1]
        eQTL_raw_sig["pos_eqtl"] = eQTL_raw_sig["pos_eqtl"].apply(int)
        eQTL_raw_sig["pos_gene"] = eQTL_raw_sig["pos_eqtl"] - eQTL_raw_sig["tss_distance"]
        
        # take abs of distance
        eQTL_raw_sig["tss_distance"] = eQTL_raw_sig["tss_distance"].abs()
        
        # expand region to 5kb
        eQTL_raw_sig["eqtl_start"] = eQTL_raw_sig["pos_eqtl"] - 2500
        eQTL_raw_sig["eqtl_end"] = eQTL_raw_sig["pos_eqtl"] + 2500
        
        eQTL_raw_sig["gene_start"] = eQTL_raw_sig["pos_gene"] - 2500
        eQTL_raw_sig["gene_end"] = eQTL_raw_sig["pos_gene"] + 2500
        
        eQTL_raw_sig = eQTL_raw_sig[(eQTL_raw_sig["eqtl_start"] > 0) & (eQTL_raw_sig["gene_start"] > 0)]
        print ("number of signif eQTL-gene pairs with start pos > 2500 : {}".format(eQTL_raw_sig.shape[0]))
        eQTL_raw_sig["pair_id"] = ["pair_{}".format(i) for i in range(eQTL_raw_sig.shape[0])]
        
        # pval stats
        # print(eQTL_raw_sig["pval_nominal"].describe())
        # pval_cutoff = eQTL_raw_sig["pval_nominal"].max()
        
        # reformat gene ID
        eQTL_raw_sig["geneID_new"] = eQTL_raw_sig["gene_id"].str.extract(r'(.+)\.+')
        sig_gene_list = eQTL_raw_sig[["geneID_new"]].drop_duplicates()
        print ("number of sig genes: {}".format(sig_gene_list.shape[0]))
        #!!! sample numbers of sig genes so it's the same size as HFB sig genes
        sig_gene_list_shuf = sig_gene_list.sample(n=HFB_sig_gene.shape[0])
        # sig_gene_overlap = pd.merge(sig_gene_list, HFB_sig_gene, left_on=["geneID_new"], right_on=["geneid"], how="inner")
        #### sig eqtl-gene pairs controlled on number of sig genes
        eqtl_sig_shuf = pd.merge(sig_gene_list_shuf, eQTL_raw_sig, on="geneID_new", how="inner")
        print("eqtl shuf pval stats\n", eqtl_sig_shuf["pval_nominal"].describe())
        pval_cutoff_shuf = eqtl_sig_shuf["pval_nominal"].max()
        
        # eQTL_sig_eqtl_region = eQTL_raw_sig[["chr", "eqtl_start", "eqtl_end", "variant_id", "pair_id", "pval_nominal"]]
        # eQTL_sig_gene_region = eQTL_raw_sig[["chr", "gene_start", "gene_end", "pair_id", "geneID_new", "pval_nominal"]]

        eQTL_sig_shuf_eqtl_region = eqtl_sig_shuf[["chr", "eqtl_start", "eqtl_end", "variant_id", "pair_id", "pval_nominal"]]
        eQTL_sig_shuf_gene_region = eqtl_sig_shuf[["chr", "gene_start", "gene_end", "pair_id", "geneID_new", "pval_nominal"]]
        eQTL_sig_shuf_eqtl_region.to_csv("{}_sig_eqtl_shuf_snp_region".format(tissue), index=False, header=False, sep="\t")
        eQTL_sig_shuf_gene_region.to_csv("{}_sig_eqtl_shuf_gene_region".format(tissue), index=False, header=False, sep="\t")
        
        # !!! distance distribution
        # sns.violinplot(abs(eqtl_sig_shuf["tss_distance"]), orient="v")
        distance_bins = list(np.concatenate(
        (np.arange(10000, 80000, 10000), np.arange(80000, 200000, 20000), np.arange(200000, 300000, 50000), np.arange(300000, 600000, 150000), np.arange(600000, 1000001, 400000)
         )))
        eqtl_sig_shuf_binned = eqtl_sig_shuf.groupby(pd.cut(eqtl_sig_shuf["tss_distance"], distance_bins))
        print ("distance distribution of signif pairs: ")
        print (eqtl_sig_shuf_binned.count()["tss_distance"])
        # determine size of each bin in distance-matched sampling
        sample_size = [i*int(power) for i in eqtl_sig_shuf_binned.count()["tss_distance"].tolist()]
        
## ----------------------------- complete data --------------------------------
        if process_raw_pair:
                eQTL_raw_all = pd.read_table(all_file, delim_whitespace=True)
                    
                # keep only eQTL-TSS pairs with a tss_distance of at least 10kb
                eQTL_raw_all = eQTL_raw_all[eQTL_raw_all["tss_distance"].abs() > 10000]
                print ("number of all eQTL-gene pairs (dist > 1kb) : {}".format(eQTL_raw_all.shape[0]))
                
                eQTL_raw_all["chr"] = "chr" + eQTL_raw_all["variant_id"].str.split("_", expand=True)[0]
                eQTL_raw_all["pos_eqtl"] = eQTL_raw_all["variant_id"].str.split("_", expand=True)[1]
                eQTL_raw_all["pos_eqtl"] = eQTL_raw_all["pos_eqtl"].apply(int)
                eQTL_raw_all["pos_gene"] = eQTL_raw_all["pos_eqtl"] - eQTL_raw_all["tss_distance"]
                
                # takes abs of distance
                eQTL_raw_all["tss_distance"] = eQTL_raw_all["tss_distance"].abs()
                
                # expand region to 5kb
                eQTL_raw_all["eqtl_start"] = eQTL_raw_all["pos_eqtl"] - 2500
                eQTL_raw_all["eqtl_end"] = eQTL_raw_all["pos_eqtl"] + 2500
                
                eQTL_raw_all["gene_start"] = eQTL_raw_all["pos_gene"] - 2500
                eQTL_raw_all["gene_end"] = eQTL_raw_all["pos_gene"] + 2500
                
                eQTL_raw_all = eQTL_raw_all[(eQTL_raw_all["eqtl_start"] > 0) & (eQTL_raw_all["gene_start"] > 0)]
                print ("number of all eQTL-gene pairs with start pos > 2500 : {}".format(eQTL_raw_all.shape[0]))
                eQTL_raw_all["pair_id"] = ["pair_{}".format(i) for i in range(eQTL_raw_all.shape[0])]
                
                # reformat gene ID
                eQTL_raw_all["geneID_new"] = eQTL_raw_all["gene_id"].str.extract(r'(.+)\.+')
                # control on sig genes
                sig_genes_all_pair = pd.merge(sig_gene_list_shuf, eQTL_raw_all, on=["geneID_new"], how="inner")
                sig_genes_all_pair_snp_region = sig_genes_all_pair[["chr", "eqtl_start", "eqtl_end", "variant_id", "pair_id", "pval_nominal"]]
                sig_genes_all_pair_gene_region = sig_genes_all_pair[["chr", "gene_start", "gene_end", "pair_id", "geneID_new", "pval_nominal"]]
                print ("number of pairs controlled on sampled sig-genes: {}".format(sig_genes_all_pair.shape[0]))
                sig_genes_all_pair_snp_region.to_csv('{}{}.sig_genes_all_pair_snp_region'.format(data_dir, tissue),index=False, header=False, sep="\t")
                sig_genes_all_pair_gene_region.to_csv('{}{}.sig_genes_all_pair_gene_region'.format(data_dir, tissue),index=False, header=False, sep="\t")
                
                #!!! distance-matched sampling
                sig_genes_all_pair_binned = sig_genes_all_pair.groupby(pd.cut(sig_genes_all_pair["tss_distance"], distance_bins))
                print ("distance distribution before sampling (all pairs): ")
                print (sig_genes_all_pair_binned.count()["tss_distance"])
                
                sig_genes_all_pair_binned_df_list = [group.reset_index().drop(["index"], axis=1) for _, group in sig_genes_all_pair_binned]    
                sig_genes_all_pair_sampled_list = [sig_genes_all_pair_binned_df_list[i].sample(n=sample_size[i]) for i in range(len(sample_size))]                    
                sig_genes_all_pair_sampled_concat = pd.concat(sig_genes_all_pair_sampled_list).reset_index().drop(["index"], axis=1)
                sig_genes_all_pair_sampled_binned = sig_genes_all_pair_sampled_concat.groupby(pd.cut(sig_genes_all_pair_sampled_concat["tss_distance"], distance_bins))
                print ("distance distribution after sampling (all pairs): ")
                print(sig_genes_all_pair_sampled_binned.count()["tss_distance"])
    
                sig_genes_all_pair_sampled_concat_snp_region = sig_genes_all_pair_sampled_concat[["chr", "eqtl_start", "eqtl_end", "variant_id", "pair_id", "pval_nominal"]]
                sig_genes_all_pair_sampled_concat_snp_region.to_csv("{}{}.sig_genes_all_pair_snp_region_dist_matched".format(data_dir, tissue), index=False, header=False, sep="\t")
                sig_genes_all_pair_sampled_concat_gene_region = sig_genes_all_pair_sampled_concat[["chr", "gene_start", "gene_end", "pair_id", "geneID_new", "pval_nominal"]]
                sig_genes_all_pair_sampled_concat_gene_region.to_csv("{}{}.sig_genes_all_pair_gene_region_dist_matched".format(data_dir, tissue), index=False, header=False, sep="\t")
                
                
                ### filter for non-sig pairs
                eQTL_nonsig = eQTL_raw_all[eQTL_raw_all["pval_nominal"] > pval_cutoff_shuf]
                # control on sig genes
                sig_genes_nonsig_pair = pd.merge(sig_gene_list_shuf, eQTL_nonsig, on=["geneID_new"], how="inner")
                sig_genes_nonsig_pair_snp_region = sig_genes_all_pair[["chr", "eqtl_start", "eqtl_end", "variant_id", "pair_id", "pval_nominal"]]
                sig_genes_nonsig_pair_gene_region = sig_genes_all_pair[["chr", "gene_start", "gene_end", "pair_id", "geneID_new", "pval_nominal"]]
                print ("number of nonsig pairs controlled on sampled sig-genes: {}".format(sig_genes_nonsig_pair.shape[0]))
                sig_genes_nonsig_pair_snp_region.to_csv('{}{}.sig_genes_nonsig_pair_snp_region'.format(data_dir, tissue),index=False, header=False, sep="\t")
                sig_genes_nonsig_pair_gene_region.to_csv('{}{}.sig_genes_nonsig_pair_gene_region'.format(data_dir, tissue),index=False, header=False, sep="\t")


                #!!! distance-matched sampling
                sig_genes_nonsig_pair_binned = sig_genes_nonsig_pair.groupby(pd.cut(sig_genes_nonsig_pair["tss_distance"], distance_bins))
                print ("distance distribution before sampling (nonsig pairs): ")
                print (sig_genes_nonsig_pair_binned.count()["tss_distance"])
                
                sig_genes_nonsig_pair_binned_df_list = [group.reset_index().drop(["index"], axis=1) for _, group in sig_genes_nonsig_pair_binned]    
                sig_genes_nonsig_pair_sampled_list = [sig_genes_nonsig_pair_binned_df_list[i].sample(n=sample_size[i]) for i in range(len(sample_size))]                    
                sig_genes_nonsig_pair_sampled_concat = pd.concat(sig_genes_nonsig_pair_sampled_list).reset_index().drop(["index"], axis=1)
                sig_genes_nonsig_pair_sampled_binned = sig_genes_nonsig_pair_sampled_concat.groupby(pd.cut(sig_genes_nonsig_pair_sampled_concat["tss_distance"], distance_bins))
                print ("distance distribution after sampling (nonsig pairs): ")
                print(sig_genes_nonsig_pair_sampled_binned.count()["tss_distance"])
    
                sig_genes_nonsig_pair_sampled_concat_snp_region = sig_genes_nonsig_pair_sampled_concat[["chr", "eqtl_start", "eqtl_end", "variant_id", "pair_id", "pval_nominal"]]
                sig_genes_nonsig_pair_sampled_concat_snp_region.to_csv("{}{}.sig_genes_nonsig_pair_snp_region_dist_matched".format(data_dir, tissue), index=False, header=False, sep="\t")
                sig_genes_nonsig_pair_sampled_concat_gene_region = sig_genes_nonsig_pair_sampled_concat[["chr", "gene_start", "gene_end", "pair_id", "geneID_new", "pval_nominal"]]
                sig_genes_nonsig_pair_sampled_concat_gene_region.to_csv("{}{}.sig_genes_nonsig_pair_gene_region_dist_matched".format(data_dir, tissue), index=False, header=False, sep="\t")
                    
                df_dist_for_plot = pd.DataFrame()
                df_dist_for_plot["sig_gene_all"] = abs(sig_genes_all_pair_sampled_concat["tss_distance"]).apply(float)
                df_dist_for_plot["sig_gene_nonsig"] = abs(sig_genes_nonsig_pair_sampled_concat["tss_distance"]).apply(float)
                df_dist_for_plot["sig"] = abs(eqtl_sig_shuf["tss_distance"]).apply(float)
                df_dist_for_plot.to_csv("{}_sig_gene_shuf_dist_matched_for_plot".format(tissue), index=False, sep="\t")
                
reform_eQTL(data_dir=sys.argv[1], tissue=sys.argv[2], process_raw_pair=sys.argv[3], HFB_sig_gene_file=sys.argv[4], power=sys.argv[5])