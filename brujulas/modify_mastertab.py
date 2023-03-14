# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 18:46:06 2023

@author: Pablo
"""



import pandas as pd
from tqdm import tqdm

FILENAME = "master_tab.tsv"
ID_COLUMN_NAME = "MP"

df = pd.read_csv(FILENAME, sep=",")

df=df.drop( df.columns[0] , axis="columns")


contexts = []

bathy_list  = ["MP0313","MP0315","MP0530","MP0532","MP0534","MP0780","MP0782",
               "MP784", # <<---- missing 0?
               "MP0784","MP0878","MP0880","MP0882","MP1154","MP1162","MP1409",
               "MP1411","MP1519","MP1521","MP1676","MP1674","MP1847","MP1845",
               "MP2231","MP2233bis","MP2809","MP2811" ]


meso_list = ["MP0317","MP0319","MP0536","MP0538","MP0786","MP0788","MP0884",
             "MP0886","MP1164","MP1166","MP1178","MP1413","MP1415","MP1417",
             "MP1523","MP1164","MP1525","MP1677","MP1678","MP1680","MP1681",
             "MP1682","MP1849","MP1851","MP1853","MP2235","MP2237","MP2817",
             "MP2813","MP2815" ]

epi_list= ["MP2239","MP2241","MP0323", "MP2819", "MP1857","MP0311","MP1419",
           "MP1421","MP1517","MP1527","MP1529","MP0790","MP0888","MP0778",
           "MP0528","MP0540","MP1176", "MP1174", "MP2821","MP1855" , "MP1684",
           "MP0321","MP2243","MP1672"]

for jj, row in tqdm(df.iterrows(), total=len(df)):
# row = df.iloc[4]
    if row[ID_COLUMN_NAME] in bathy_list:
        contexts.append("Bathypelagic")
    elif row[ID_COLUMN_NAME] in meso_list:
        contexts.append("Mesopelagic")
    elif row[ID_COLUMN_NAME] in epi_list:
        contexts.append("Epipelagic")
    else:
        print( f"Couldnt find a context for sample id: {row['Sample_ID']}")
    

df["Context"] = contexts
df.to_csv("mastertable_w_ctxt.tsv", sep="\t", index=False)





