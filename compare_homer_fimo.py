# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 09:56:48 2019

@author: libin
"""

import pandas as pd
import numpy as np

cell = "iN"
homer = pd.read_csv("C:\\Users\\libin\\UCSF\\motif_analysis\\FIMO&HOMER\\{}_homer.txt".format(cell), sep="\t")
jaspar = pd.read_csv("C:\\Users\\libin\\UCSF\\motif_analysis\\FIMO&HOMER\\{}_jaspar.txt".format(cell), sep="\t")
