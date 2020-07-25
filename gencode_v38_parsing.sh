# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 00:46:56 2020

@author: bingkun
@ some basic formatting for GENCODE.V32 GTF
"""


import pandas as pd

#%% transcript_id to gene_id
gtf_v38 = pd.read_csv(r'C:\Users\libin\UCSF\hg38_general\gencode.v32.annotation.gtf', sep="\t", comment="#",
                      names=["chr", "source", "type", "start", "end", "score", "strand", "genomic_phase", "annotation"])
gtf_v38 = gtf_v38[gtf_v38["type"] == "transcript"]

gtf_v38["gene_id"] = gtf_v38["annotation"].str.extract(r'gene_id "(.+?)\.\d+";')
gtf_v38["transcript_id"] = gtf_v38["annotation"].str.extract(r'transcript_id "(.+?)\.\d+";')
gtf_v38["gene_name"] = gtf_v38["annotation"].str.extract(r'gene_name "(.+?)";')

gtf_v38 = gtf_v38[["gene_id", "transcript_id", "gene_name"]]
gtf_v38.to_csv(r'C:\Users\libin\UCSF\hg38_general\gencode.v32.transcript2gene', sep="\t", index=False, header=True)


#%% gene_id to gene_pos
gtf_v38 = pd.read_csv(r'C:\Users\libin\UCSF\hg38_general\gencode.v32.annotation.gtf', sep="\t", comment="#",
                      names=["chr", "source", "type", "start", "end", "score", "strand", "genomic_phase", "annotation"])
gtf_v38 = gtf_v38[gtf_v38["type"] == "gene"]
gtf_v38["gene_id"] = gtf_v38["annotation"].str.extract(r'gene_id "(.+?)\.\d+";')
gtf_v38 = gtf_v38[["chr", "start", "end", "gene_id"]]
gtf_v38["chr_num"] = gtf_v38["chr"].str.extract('chr(.+)')
gtf_v38["gene_pos_merge"] = gtf_v38[["chr_num", "start", "end"]].apply(lambda x: ":".join(x.map(str)), axis=1)
gtf_v38 = gtf_v38[["gene_id", "gene_pos_merge"]]
gtf_v38.to_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v32.geneID2genePOS', sep="\t", index=False, header=True)


#%% gene_id to gene_nameType
gtf_v38 = pd.read_csv(r'C:\Users\libin\UCSF\hg38_general\gencode.v32.annotation.gtf', sep="\t", comment="#",
                      names=["chr", "source", "type", "start", "end", "score", "strand", "genomic_phase", "annotation"])
gtf_v38 = gtf_v38[gtf_v38["type"] == "gene"]
gtf_v38["gene_id"] = gtf_v38["annotation"].str.extract(r'gene_id "(.+?)\.\d+";')
gtf_v38["gene_name"] = gtf_v38["annotation"].str.extract(r'gene_name "(.+?)";')
gtf_v38["gene_type"] = gtf_v38["annotation"].str.extract(r'gene_type "(.+?)";')
gtf_v38 = gtf_v38[["gene_id", "gene_name", "gene_type"]]
gtf_v38.to_csv(r'C:\Users\libin\UCSF\hg38_general\gencode.v32.geneID2geneNameType', sep="\t", index=False, header=True)
