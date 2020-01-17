# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 00:46:56 2020

@author: bingkun
@ some basic formatting for GENCODE.V19 GTF
"""


import pandas as pd

#%% transcript_id to gene_id
gtf_v19 = pd.read_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.annotation.gtf', sep="\t", comment="#",
                      names=["chr", "source", "type", "start", "end", "score", "strand", "genomic_phase", "annotation"])
gtf_v19 = gtf_v19[gtf_v19["type"] == "transcript"]

gtf_v19["gene_id"] = gtf_v19["annotation"].str.extract(r'gene_id "(.+?)\.\d+";')
gtf_v19["transcript_id"] = gtf_v19["annotation"].str.extract(r'transcript_id "(.+?)\.\d+";')
gtf_v19["gene_name"] = gtf_v19["annotation"].str.extract(r'gene_name "(.+?)";')

gtf_v19 = gtf_v19[["gene_id", "transcript_id", "gene_name"]]
gtf_v19.to_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.transcript2gene', sep="\t", index=False, header=True)


#%% gene_id to gene_pos
gtf_v19 = pd.read_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.annotation.gtf', sep="\t", comment="#",
                      names=["chr", "source", "type", "start", "end", "score", "strand", "genomic_phase", "annotation"])
gtf_v19 = gtf_v19[gtf_v19["type"] == "gene"]
gtf_v19["gene_id"] = gtf_v19["annotation"].str.extract(r'gene_id "(.+?)\.\d+";')
gtf_v19 = gtf_v19[["chr", "start", "end", "gene_id"]]
gtf_v19["chr_num"] = gtf_v19["chr"].str.extract('chr(.+)')
gtf_v19["gene_pos_merge"] = gtf_v19[["chr_num", "start", "end"]].apply(lambda x: ":".join(x.map(str)), axis=1)
gtf_v19 = gtf_v19[["gene_id", "gene_pos_merge"]]
gtf_v19.to_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.geneID2genePOS', sep="\t", index=False, header=True)


#%% gene_id to gene_nameType
gtf_v19 = pd.read_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.annotation.gtf', sep="\t", comment="#",
                      names=["chr", "source", "type", "start", "end", "score", "strand", "genomic_phase", "annotation"])
gtf_v19 = gtf_v19[gtf_v19["type"] == "gene"]
gtf_v19["gene_id"] = gtf_v19["annotation"].str.extract(r'gene_id "(.+?)\.\d+";')
gtf_v19["gene_name"] = gtf_v19["annotation"].str.extract(r'gene_name "(.+?)";')
gtf_v19["gene_type"] = gtf_v19["annotation"].str.extract(r'gene_type "(.+?)";')
gtf_v19 = gtf_v19[["gene_id", "gene_name", "gene_type"]]
gtf_v19.to_csv(r'C:\Users\libin\UCSF\hg19_general\gencode.v19.geneID2geneNameType', sep="\t", index=False, header=True)
