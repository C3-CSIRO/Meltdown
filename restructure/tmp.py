# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 14:19:09 2015

@author: ris01e
"""

import matplotlib.pyplot as plt
import pandas as pd

dataFilePath = "C:/Users/RIS01E/Documents/GitHub/Meltdown/restructure/dataSample1.txt"
controlnoprotein = "C:/Users/RIS01E/Documents/GitHub/Meltdown/restructure/data/noProteinControl.csv"
controlnodye = "C:/Users/RIS01E/Documents/GitHub/Meltdown/restructure/data/noDyeControl.csv"

data = pd.DataFrame.from_csv(dataFilePath, sep='\t', index_col='Temperature')
noprotein = pd.Series.from_csv(controlnoprotein)
tmps = noprotein.index

temperatures = list(data.index)
wells = ["H10", "H11", "H12"]

for well in wells:
    area=0
    curve = list(data[well])
    for point in curve:
        area += point*0.5
    curve = [x/area for x in curve]
    plt.plot(temperatures, curve)
plt.show()

#plt.plot(tmps, list(noprotein))
#plt.show()
