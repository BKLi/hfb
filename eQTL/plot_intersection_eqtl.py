# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 00:33:56 2020

@author: libin
"""

import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['figure.dpi']= 200

cellType="neuron"
npval_compair = pd.read_csv(r"C:\Users\libin\UCSF\hfb\eqtl\20200115\npval_compair_{}".format(cellType), sep="\t")
npval_compair_melt = pd.melt(npval_compair, value_name="npval", var_name="")

#plt.figure(figsize=(10,7))
sns.violinplot(x="", y="npval", data=npval_compair_melt, palette="Pastel1", cut=0)
