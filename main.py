"""
Código para identificação de ilhas observáveis, implementação em python
"""

#%%
from datetime import timedelta

from classes import *
from readfiles import *
from networkstruc import *
from networkcalc import *
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import time

def rgb_to_hex(r, g, b):
    return '#{:02x}{:02x}{:02x}'.format(r, g, b)



def encontracaminhos(graph,L):
    
    
    flag=0
    # analiza caminhos
    lstvar=list(range(len(graph)))
    lstpercorridas=[]
    lstemcaminho=[]
    caminhos={}
    for var in lstvar:

        if var not in lstpercorridas:
            caminhos[var]=[var]
            lstpercorridas.append(var)
            
            var2=var
            zerocol=0


            flagcruzou=False
            
            while var2<(len(graph)+flag):
   
                if zerocol==1:
                    break
                i=0
                for element in L[var2+1:,var2]:
                    zerocol=1
                    if  np.abs(element)>0:
                        teta=var2+1+i
                        caminhos[var].append(teta)
                        lstpercorridas.append(teta)
                        if teta in lstemcaminho:
                            flagcruzou=True
                        lstemcaminho.append(teta)
                        var2=teta
                        zerocol=0
                        break
                    i=i+1
                if flagcruzou==True:
                    break
    
    return caminhos


def plota_grafo(graph,caminhos,ind_i,i_ind):
    g=nx.Graph()

    for no in graph:
        g.add_node(no.id)



    for key, item in caminhos.items():
        cor=list(np.random.choice(range(256), size=3))
        hex_col=rgb_to_hex(cor[0],cor[1],cor[2])
        if len(item)>1:
            for i in range(1,len(item)):
                g.add_edge(item[i-1],item[i],color=hex_col,weight=1)


    edges = g.edges()
    colors = [g[u][v]['color'] for u,v in edges]
    weights = [g[u][v]['weight'] for u,v in edges]

    nx.draw(g,labels=i_ind,with_labels=True,edge_color=colors,width=weights)


    plt.savefig("Graph.pdf")
    plt.savefig("Graph.png")
    
    # make an undirected copy of the digraph
    subgraphs=nx.connected_components(g)

    return list(subgraphs)


def verifica_medidas_caminhos(graph,HT,dfDMED,subgraphs):
    HT[:,0]
    for i in range(len(dfDMED)):
        
        c=HT[:,i]
        nodesinvolved=[]
        for j in range(len(c)):
            if np.abs(c[j])>0:
                nodesinvolved.append(j)
        
        emquantosgrafos=0
        for graph in subgraphs:
            a=set(graph)
            b=set(nodesinvolved)
            a.intersection(b)
            if len(a.intersection(b))>0:
                emquantosgrafos=emquantosgrafos+1
        
        if emquantosgrafos>1:
            print(dfDMED.iloc[i])
        

sys = "IEEE14"

dfDBAR, dfDBRAN, dfDMED = read_files(sys)

[bars, nbars, pv, pq, ind_i] = creat_bar(dfDBAR)

[ram, nbran] = create_bran(dfDBRAN, ind_i)


i_ind= {v: k for k, v in ind_i.items()}

graph = create_graph(bars, ram)

## filtra DMED


#%%

HT = montaH(graph, dfDMED,ind_i)




Hfat,L,dfPseudo,Permutacao=fatoraH(HT,graph,dfDMED)

caminhos=encontracaminhos(graph,L)

#%%




subgraphs=plota_grafo(graph,caminhos,ind_i,i_ind)
verifica_medidas_caminhos(graph,HT,dfDMED,subgraphs)




#%%
         
                










# %%
