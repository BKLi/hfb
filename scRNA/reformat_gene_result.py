# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 22:13:57 2019

@author: libin
"""

import pandas as pd


hg38_gtf = pd.read_csv("C:\\Users\libin\\UCSF\gene_annotation\gencode.v29.annotation.gtf", sep="\t", comment="#",
                         names=["chr", "source", "type", "start", "end", "score", "strand", "genomic_phase", "annotation"])
# non-greedy regex matching -- smallest number possible
# https://stackoverflow.com/questions/9519734/python-regex-to-find-a-string-in-double-quotes-within-a-string
hg38_gtf["gene_id"] = hg38_gtf["annotation"].str.extract(r'gene_id "(.+?)\.\S+";')
hg38_gtf["gene_name"] = hg38_gtf["annotation"].str.extract(r'gene_name "(.+?)";')
hg38_gtf_cut = hg38_gtf[["gene_id","gene_name"]]
hg38_gtf_cut_gene = hg38_gtf_cut[hg38_gtf["type"] == "gene"]
# duplication due to same gene on X/Y chromosome
hg38_gtf_cut_gene = hg38_gtf_cut_gene.drop_duplicates()

gene_results = pd.read_csv("C:\\Users\libin\\UCSF\hfb\\rpkm.data.individual.txt", sep="\t")
gene_results = gene_results.reset_index().rename(columns={"index":"gene_id"})
gene_results_merged = pd.merge(hg38_gtf_cut_gene, gene_results, on=["gene_id"], how="inner")

RG_out = gene_results_merged[["gene_name", "RG"]]
IPC_out = gene_results_merged[["gene_name", "IPC"]]
eN_out = gene_results_merged[["gene_name", "eN"]]
iN_out = gene_results_merged[["gene_name", "iN"]]

RG_out.to_csv("C:\\Users\\libin\\UCSF\\hfb\\scRNA\\updated_Aug_27\\RGs.average.genes.results", sep="\t", index=False, header=True)
IPC_out.to_csv("C:\\Users\\libin\\UCSF\\hfb\\scRNA\\updated_Aug_27\\IPCs.average.genes.results", sep="\t", index=False, header=True)
eN_out.to_csv("C:\\Users\\libin\\UCSF\\hfb\\scRNA\\updated_Aug_27\\IPCs.average.genes.results", sep="\t", index=False, header=True)