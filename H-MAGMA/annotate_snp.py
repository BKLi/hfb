# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 19:10:17 2020

@author: bingkun

@project hfb -- H-MAGMA

* annotate SNP with related genes
* output annotation file required by H-MAGMA
* s2 of pipeline
"""

import pandas as pd
from collections import OrderedDict

cellType = "RGC"

#promoter_reformed = pd.read_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\gencode_v19_promoters.reform.bed', sep="\t",
#                                names=["chr", "start", "end", "gene_id", "gene_type", "gene_name"])
#promoter_raw = pd.read_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\gencode_v19_promoters.bed', sep="\t",
#                                names=["chr", "start", "end", "starnd", "gene_id", "gene_type", "gene_name"])

# fetal_annotation_snp_list = pd.read_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\Fetal_brain.genes.annot.reform.reform', sep="\t", names=["snp_id"])
# fetal_annotation_snp_list_dedup = fetal_annotation_snp_list.drop_duplicates()

transcript2gene = pd.read_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.transcript2gene', sep="\t")
ct_interactions = pd.read_csv(r'C:\Users\libin\UCSF\MAGMA\{}.5k.downsample.bedpe.withID'.format(cellType))

#%%
#### first priority: exon annotations
# same for each cell type
exon_intersection = pd.read_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\snp.exon.intersect', sep="\t",
                                names=["snp_chr", "snp_pos", "snp_end", "snp_id", "exon_chr", "exon_start", "exon_end", "exon_id", "notImportant", "strand"])
exon_intersection["transcript_id"] = exon_intersection["exon_id"].str.extract('(.+)\.\d+_exon_\d+_\d+_.*')
exon_intersection_cut = exon_intersection[['snp_id', 'transcript_id']]
# collapse on transcript ID
exonSNP2transcript = exon_intersection_cut.groupby(['transcript_id'], as_index=False).agg(",".join)

##### check for dups -- shouldn't be any at this point
# each snp get overlapped once by each transcript
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

#%%
#### second priority: promoter annotations
# same for each cell type
promoter_intersections = pd.read_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\snp.promoter.intersect', sep="\t",
                                    names=["snp_chr", "snp_pos", "snp_end", "snp_id", "gene_chr", "gene_start", "gene_end", "gene_id", "gene_type", "gene_name"])
# there are dups in gencode promoter file
# bc annotation is transcript-based instead of gene-based
# so each gene may have more than one promoters
# and those promoters may have same or diff position
# I'm removing dups subsetting on ["snp_id", 'gene_id']
# so as to transform annotation into gene-based from transcript-based
promoter_intersections_dedup = promoter_intersections.drop_duplicates(subset=['snp_id', 'gene_id'])
promoter_intersections_cut = promoter_intersections_dedup[['snp_id', 'gene_id']]
promoterSNP2gene = promoter_intersections_cut.groupby(['gene_id'], as_index=False).agg(",".join)
print ("duplicates: ", promoterSNP2gene[["gene_id"]].drop_duplicates().shape[0] - promoterSNP2gene.shape[0])

#%%
#### target bin (interactions) annotations
# cell-type specific
anchor_promoter = pd.read_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\{}.anchor_promoter.intersect'.format(cellType), sep="\t",
                              names=["ac_chr", "ac_start", "ac_end", "inter_id", "ac_id", "pr_chr", "pr_start", "pr_end", "gene_id", "gene_type", "gene_name"])
anchor_promoter_dedup = anchor_promoter.drop_duplicates(subset=["inter_id", "gene_id"])
anchor_promoter_cut = anchor_promoter_dedup[["inter_id", "ac_id", "gene_id"]]

snp_target = pd.read_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\snp.{}_target.intersect'.format(cellType), sep="\t",
                         names=["snp_chr", "snp_pos", "snp_end", "snp_id", "tg_chr", "tg_start", "tg_end", "inter_id", "ac_id"])
snp_target_cut = snp_target[["inter_id", "ac_id", "snp_id"]]
snp_target_cut = snp_target_cut.groupby(['inter_id', 'ac_id'], as_index=False).agg(",".join)
# each interaction may intersect more than one promoter
targetSNP2gene_tmp = pd.merge(anchor_promoter_cut, snp_target_cut, on=["inter_id", "ac_id"], how="inner")
targetSNP2gene_tmp = targetSNP2gene_tmp[["snp_id", "gene_id"]]
# collapse on gene ID
targetSNP2gene = targetSNP2gene_tmp.groupby(['gene_id'], as_index=False).agg(",".join)
print ("duplicates: ", targetSNP2gene[["gene_id"]].drop_duplicates().shape[0] - targetSNP2gene.shape[0])

#%% concat three types of annotations
SNP2gene_all_tmp = pd.concat([exonSNP2gene, promoterSNP2gene, targetSNP2gene], axis=0)
SNP2gene_all = SNP2gene_all_tmp.groupby(['gene_id'], as_index=False).agg(",".join)
SNP2gene_all["snp_id_dedup"] = SNP2gene_all["snp_id"].str.split(",").apply(lambda x: OrderedDict.fromkeys(x).keys()).str.join('\t')

geneID2genePOS = pd.read_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.geneID2genePOS', sep="\t")
SNP2gene_all_with_pos = pd.merge(geneID2genePOS, SNP2gene_all, on=["gene_id"], how="inner")
SNP2gene_all_with_pos = SNP2gene_all_with_pos[['gene_id', 'gene_pos_merge', 'snp_id_dedup']]
SNP2gene_all_with_pos.to_csv(r'C:\Users\libin\UCSF\hfb\HMAGMA\{}.genes.annot'.format(cellType), sep=",", index=False, header=False)
