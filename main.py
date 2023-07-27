"""
Código para identificação de ilhas observáveis, implementação em python
"""

#%%
from datetime import timedelta

from classes import *
from readfiles import *
from networkstruc import *
from networkcalc import *
import numpy as np
import time

sys = "IEEE14"

dfDBAR, dfDBRAN, dfDMED = read_files(sys)

[bars, nbars, pv, pq, ind_i] = creat_bar(dfDBAR)

[ram, nbran] = create_bran(dfDBRAN, ind_i)

graph = create_graph(bars, ram)

## filtra DMED


#%%

HT = montaH(graph, dfDMED,ind_i)




Hfat,L,dfPseudo,Permutacao=fatoraH(HT,graph,dfDMED)



# analiza caminhos
lstvar=list(range(graph))
lstpercorridas=[]
caminhos={}
for var in lstvar:

    if var not in lstpercorridas:
        caminhos[var]=[var]
        lstpercorridas.append(var)
        i=0
        for element in L[var+1:,var]:
            if  np.abs(element)>0:
                teta=var+1+i
                caminhos[var].append(teta)
                lstpercorridas.append(teta)










# %%
