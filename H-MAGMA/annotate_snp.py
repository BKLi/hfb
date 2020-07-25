# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 19:10:17 2020

@author: bingkun

@project  H-MAGMA (hfb / BEST1)

* annotate SNP to related genes
* output annotation file required by H-MAGMA
* s2 of pipeline

next step : sed -i 's/,/\t/g' *.genes.annot
"""

import pandas as pd
pd.__version__
from collections import OrderedDict


cellType_1 = "KJ184_KJ211_KJ189_KJ212"
cellType_2 = "RPE"
reformat_snp = False
target_only = False
data_dir = r"C:\Users\libin\UCSF\MAGMA\for_kirsty"

#promoter_reformed = pd.read_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\gencode_v19_promoters.reform.bed', sep="\t",
#                                names=["chr", "start", "end", "gene_id", "gene_type", "gene_name"])
#promoter_raw = pd.read_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\gencode_v19_promoters.bed', sep="\t",
#                                names=["chr", "start", "end", "starnd", "gene_id", "gene_type", "gene_name"])

# fetal_annotation_snp_list = pd.read_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\Fetal_brain.genes.annot.reform.reform', sep="\t", names=["snp_id"])
# fetal_annotation_snp_list_dedup = fetal_annotation_snp_list.drop_duplicates()


transcript2gene = pd.read_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.transcript2gene', sep="\t")
ct_interactions = pd.read_csv(r'{}\{}_final_filtered_interactions.bed.withID'.format(data_dir, cellType_1))


if reformat_snp:
        snp = pd.read_csv(r'C:\Users\libin\UCSF\MAGMA\for_kirsty\efotraits_EFO_0001365-associations-2020-07-20.csv', sep=",")
        #rs71507014-<b>?</b>
        #6:43860845
        snp = snp[snp["Location"] != "Mapping not available"]
        snp = snp[~snp["Location"].str.contains("MT")]
        snp["snpID"] = snp['Variant and risk allele'].str.extract(r'(.+)-')
        snp["P-value"] = snp["P-value"].str.replace(" ","")
        snp["P-value"] = snp["P-value"].str.replace("x10","E")
        snp["P-value"] = snp["P-value"].astype("float64")

        snp["chr"] = "chr" + snp["Location"].str.extract(r'(.*):')
        snp["start"] = snp["Location"].str.extract(r':(.*)').astype(int)
        snp["end"] = snp["start"] + 1
        snp = snp[["chr", "start", "end", "snpID", 'P-value']]
        snp = snp.sort_values(by=["snpID", "P-value"])
        snp = snp.drop_duplicates(subset=["snpID"], keep='first')
        snp.to_csv(r'C:\Users\libin\UCSF\MAGMA\for_kirsty\AMD_snp.hg38.bed', sep="\t", index=False, header=True)
        


if target_only:
        #### target bin (interactions) annotations
        # cell-type specific
        anchor_promoter = pd.read_csv(r'{}\{}_anchor.promoter.intersect'.format(data_dir, cellType_2), sep="\t",
                                      names=["ac_chr", "ac_start", "ac_end", "inter_id", "ac_id", "pr_chr", "pr_start", "pr_end", "gene_id", "gene_type", "gene_name"])
        # tx-based to gene-based
        anchor_promoter_dedup = anchor_promoter.drop_duplicates(subset=["inter_id", "ac_id", "gene_id"])
        anchor_promoter_cut = anchor_promoter_dedup[["inter_id", "ac_id", "gene_id"]]
        
        snp_target = pd.read_csv(r'{}\AMD_snp.{}_target.intersect'.format(data_dir,cellType_2), sep="\t",
                                 names=["snp_chr", "snp_pos", "snp_end", "snp_id", "snp_pval", "tg_chr", "tg_start", "tg_end", "inter_id"])
        snp_target_cut = snp_target[["inter_id", "snp_id"]]
        snp_target_cut = snp_target_cut.groupby(['inter_id'], as_index=False).agg(",".join)
        # each interaction may intersect more than one promoter
        targetSNP2gene_tmp = pd.merge(anchor_promoter_cut, snp_target_cut, on=["inter_id"], how="inner")
        targetSNP2gene_tmp = targetSNP2gene_tmp[["snp_id", "gene_id"]]
        # collapse on gene ID
        targetSNP2gene = targetSNP2gene_tmp.groupby(['gene_id'], as_index=False).agg(",".join)
        print ("gene duplicates: ", targetSNP2gene[["gene_id"]].drop_duplicates().shape[0] - targetSNP2gene.shape[0])
        print ("number of genes involved (target only): ", targetSNP2gene.shape[0])
        
        targetSNP2gene["snp_id_dedup"] = targetSNP2gene["snp_id"].str.split(",").apply(lambda x: OrderedDict.fromkeys(x).keys()).str.join('\t')
        targetSNP2gene["snp_id_dedup_list"] = targetSNP2gene["snp_id_dedup"].str.split("\t")
        
        snp_list_tmp = targetSNP2gene.set_index(['gene_id', 'snp_id', 'snp_id_dedup']).apply(lambda x: x.explode()).reset_index()
        snp_list = snp_list_tmp['snp_id_dedup_list'].drop_duplicates().values.tolist()
        print ("number of SNPs annotated (target only): ", len(snp_list))
        
        geneID2genePOS = pd.read_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.geneID2genePOS', sep="\t")
        targetSNP2gene_with_pos = pd.merge(geneID2genePOS, targetSNP2gene, on=["gene_id"], how="inner")
        targetSNP2gene_with_pos = targetSNP2gene_with_pos[['gene_id', 'gene_pos_merge', 'snp_id_dedup']]
        targetSNP2gene_with_pos.to_csv(r'{}\{}.genes.target.annot'.format(data_dir, cellType_1), sep=",", index=False, header=False)
        
        geneID2geneName = pd.read_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.geneID2geneNameType', sep="\t")
        targetSNP2gene_with_pos_name = pd.merge(targetSNP2gene_with_pos, geneID2geneName, on=["gene_id"], how="inner")


#### first priority: exon annotations
# same for each cell type
        
else:
        exon_intersection = pd.read_csv(r'{}\AMD_snp.exon.intersect'.format(data_dir), sep="\t",
                                        names=["snp_chr", "snp_pos", "snp_end", "snp_id", "snp_pval", "exon_chr", "exon_start", "exon_end", "exon_id", "notImportant", "strand"])
        exon_intersection["transcript_id"] = exon_intersection["exon_id"].str.extract('(.+)\.\d+_exon_\d+_\d+_.*')
        exon_intersection_cut = exon_intersection[['snp_id', 'transcript_id']]
        # collapse on transcript ID
        exonSNP2transcript = exon_intersection_cut.groupby(['transcript_id'], as_index=False).agg(",".join)
        
        ##### check for dups -- shouldn't be any at this point
        # each snp can only be overlapped only once by each transcript
        exonSNP2transcript["num_snp"] = exonSNP2transcript["snp_id"].str.split(",").apply(lambda x: len(x))
        # https://stackoverflow.com/questions/47316783/python-dataframe-remove-duplicate-words-in-the-same-cell-within-a-column-in-pyt
        exonSNP2transcript["snp_id_dedup"] = exonSNP2transcript["snp_id"].str.split(",").apply(lambda x: OrderedDict.fromkeys(x).keys()).str.join(',')
        exonSNP2transcript["num_snp_dedup"] = exonSNP2transcript["snp_id_dedup"].str.split(",").apply(lambda x: len(x))
        exonSNP2transcript["num_snp_diff"] = exonSNP2transcript["num_snp_dedup"] - exonSNP2transcript["num_snp"]
        print ("duplicates: ", exonSNP2transcript[exonSNP2transcript["num_snp_diff"] !=0].shape[0]) #should be 0
        
        exonSNP2gene_tmp = pd.merge(exonSNP2transcript, transcript2gene, on=["transcript_id"], how="inner")
        exonSNP2gene_tmp = exonSNP2gene_tmp[['snp_id', 'gene_id']]
        # collapse on gene ID
        exonSNP2gene = exonSNP2gene_tmp.groupby(['gene_id'], as_index=False).agg(",".join)
        # check for duplicates
        print ("duplicates: ", exonSNP2gene[["gene_id"]].drop_duplicates().shape[0] - exonSNP2gene.shape[0])
        
        
        #### second priority: promoter annotations
        # same for each cell type
        promoter_intersections = pd.read_csv(r'{}\AMD_snp.promoter.intersect'.format(data_dir), sep="\t",
                                            names=["snp_chr", "snp_pos", "snp_end", "snp_id", "snp_pval", "gene_chr", "gene_start", "gene_end", "gene_id", "gene_type", "gene_name"])
        # there are dups in gencode promoter file
        # bc annotation is transcript-based instead of gene-based
        # so each gene may have more than one promoters
        # and those promoters may have same or diff positions
        # I'm removing dups subsetting on ["snp_id", 'gene_id']
        # so as to transform annotation into gene-based from transcript-based
        promoter_intersections_dedup = promoter_intersections.drop_duplicates(subset=['snp_id', 'gene_id'])
        promoter_intersections_cut = promoter_intersections_dedup[['snp_id', 'gene_id']]
        promoterSNP2gene = promoter_intersections_cut.groupby(['gene_id'], as_index=False).agg(",".join)
        print ("duplicates: ", promoterSNP2gene[["gene_id"]].drop_duplicates().shape[0] - promoterSNP2gene.shape[0])
        
        
        #### target bin (interactions) annotations
        # cell-type specific
        anchor_promoter = pd.read_csv(r'{}\{}_anchor.promoter.intersect'.format(data_dir, cellType_2), sep="\t",
                                      names=["ac_chr", "ac_start", "ac_end", "inter_id", "ac_id", "pr_chr", "pr_start", "pr_end", "gene_id", "gene_type", "gene_name"])
        # tx-based to gene-based
        anchor_promoter_dedup = anchor_promoter.drop_duplicates(subset=["inter_id", "ac_id", "gene_id"])
        anchor_promoter_cut = anchor_promoter_dedup[["inter_id", "ac_id", "gene_id"]]
        
        snp_target = pd.read_csv(r'{}\AMD_snp.{}_target.intersect'.format(data_dir,cellType_2), sep="\t",
                                 names=["snp_chr", "snp_pos", "snp_end", "snp_id", "snp_pval", "tg_chr", "tg_start", "tg_end", "inter_id"])
        snp_target_cut = snp_target[["inter_id", "snp_id"]]
        snp_target_cut = snp_target_cut.groupby(['inter_id'], as_index=False).agg(",".join)
        # each interaction may intersect more than one promoter
        targetSNP2gene_tmp = pd.merge(anchor_promoter_cut, snp_target_cut, on=["inter_id"], how="inner")
        targetSNP2gene_tmp = targetSNP2gene_tmp[["snp_id", "gene_id"]]
        # collapse on gene ID
        targetSNP2gene = targetSNP2gene_tmp.groupby(['gene_id'], as_index=False).agg(",".join)
        print ("gene duplicates: ", targetSNP2gene[["gene_id"]].drop_duplicates().shape[0] - targetSNP2gene.shape[0])
 
        
        # concat three types of annotations
        SNP2gene_all_tmp = pd.concat([exonSNP2gene, promoterSNP2gene, targetSNP2gene], axis=0)
        SNP2gene_all = SNP2gene_all_tmp.groupby(['gene_id'], as_index=False).agg(",".join)
        print ("number of genes involved: ", SNP2gene_all.shape[0])
        
        SNP2gene_all["snp_id_dedup"] = SNP2gene_all["snp_id"].str.split(",").apply(lambda x: OrderedDict.fromkeys(x).keys()).str.join('\t')
        SNP2gene_all["snp_id_dedup_list"] = SNP2gene_all["snp_id_dedup"].str.split("\t")
        
        snp_list_tmp = SNP2gene_all.set_index(['gene_id', 'snp_id', 'snp_id_dedup']).apply(lambda x: x.explode()).reset_index()
        snp_list = snp_list_tmp['snp_id_dedup_list'].drop_duplicates().values.tolist()
        print ("number of SNPs annotated: ", len(snp_list))
        
        
        geneID2genePOS = pd.read_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.geneID2genePOS', sep="\t")
        SNP2gene_all_with_pos = pd.merge(geneID2genePOS, SNP2gene_all, on=["gene_id"], how="inner")
        SNP2gene_all_with_pos = SNP2gene_all_with_pos[['gene_id', 'gene_pos_merge', 'snp_id_dedup']]
        SNP2gene_all_with_pos.to_csv(r'{}\{}.genes.annot'.format(data_dir, cellType_1), sep=",", index=False, header=False)
        
        geneID2geneName = pd.read_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.geneID2geneNameType', sep="\t")
        SNP2gene_all_with_pos_name = pd.merge(SNP2gene_all_with_pos, geneID2geneName, on=["gene_id"], how="inner")
