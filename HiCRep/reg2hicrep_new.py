# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 17:39:48 2019

@author: libin
"""
#start = "105000"
#end = "58585000"
#res = "5000"
#chr="chr19"
#infile = r'C:\Users\libin\Scripts\hfb\MS084.chr19.merged.cut'
#outfile = r'C:\Users\libin\Scripts\hfb\reg_raw.chr19.MS084.matrix'

import pandas as pd
import numpy as np
import sys

def reg2hicrep(start, end, res, chr, infile, outfile):
    
    start = int(start)
    end = int(end)
    res = int(res)
    map_out = pd.read_csv(infile, sep="\t", names=["bin1_mid", "bin2_mid", "count"],
    dtype = {"bin1_mid": np.float64, "bin2_mid": np.float64, "count": np.float64})
    #cols = ["bin1_mid", "bin2_mid", "count"]
    #map_out = map_out[~map_out["bin1_mid"].str.contains("bin2_mid")]
    #for col in cols:
        #map_out[col] = pd.to_numeric(map_out[col])
    # map_out = map_out[['bin1_mid','bin2_mid','count']]
    # map_out = map_out.astype(float)
    map_out = map_out.astype(int)
    print(map_out.head())

    #map_out["bin1_mid"] = map_out["bin1_mid"].astype('float64').astype('int64')
    #map_out["bin2_mid"] = map_out["bin2_mid"].astype('float64').astype('int64')
    
    # start = map_out['bin1_mid'].min().astype('int64')
    # print(start, type(start))
    # end = map_out['bin2_mid'].max().astype('int64')
    # print(end, type(end))
    
    dim = int((end - start)/res) + 1
    
    bins = []
    i = 0
    while i < dim:
        bins.append(np.arange(i * int(res), i * int(res) + int(res) + 1, int(res)).tolist())
        i += 1
        
    bins_df = pd.DataFrame({'bins': bins})
    bins_df[["start","end"]] = pd.DataFrame(bins_df["bins"].values.tolist(), index=bins_df.index)
    bins_df = bins_df.drop("bins", 1)
    bins_df["chr"] = chr
    bins_df = bins_df[["chr", "start", "end"]]
    
    map_out_convert = pd.DataFrame()
    map_out_convert["bin1_mid"] = map_out['bin1_mid'].apply(lambda x: (x-start)/res).astype('int64')
    map_out_convert["bin2_mid"] = map_out['bin2_mid'].apply(lambda x: (x-start)/res).astype('int64')
    map_out_convert["count"] = map_out['count']
    shape=map_out_convert.shape[0]
    
    initial_matrix = pd.DataFrame(np.zeros(shape=(dim,dim), dtype=int))
    
    for i in range(0,shape):
        bin1 = map_out_convert["bin1_mid"].iloc[i]
        bin2 = map_out_convert["bin2_mid"].iloc[i]
        count = map_out_convert["count"].iloc[i]
        initial_matrix.iat[bin1, bin2] = count
        
    map_out_full = pd.merge(bins_df, initial_matrix, left_index=True, right_index=True)
    map_out_full.to_csv(outfile, sep="\t", header=False, index=False)

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("usage : reg2hicrep")
        sys.exit()

reg2hicrep(start=sys.argv[1], end=sys.argv[2], res=sys.argv[3], chr=sys.argv[4], infile=sys.argv[5], outfile=sys.argv[6])