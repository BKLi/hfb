# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 20:32:58 2019

@author: libin
"""

import pandas as pd
import numpy as np
import os
import glob


path = r"C:\Users\libin\UCSF\hfb\reg\matrix"


min_list = []
max_list = []
for filename in glob.glob(os.path.join(path, '*.cut')):
    print(filename)
    infile = pd.read_csv(filename, sep="\t", names=["bin1_mid", "bin2_mid", "count"],
    dtype = {"bin1_mid": np.float64, "bin2_mid": np.float64, "count": np.float64})
    bin1_mid_min = infile["bin1_mid"].min()
    bin2_mid_min = infile["bin2_mid"].min()
    bin1_mid_max = infile["bin1_mid"].max()
    bin2_mid_max = infile["bin2_mid"].max()
    min_list.append(bin1_mid_min)
    min_list.append(bin2_mid_min)
    max_list.append(bin1_mid_max)
    max_list.append(bin2_mid_max)

print(int(min(min_list)))
print(int(max(max_list)))