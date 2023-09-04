"""
Código para identificação de ilhas observáveis, implementação em python
"""

#%%
from datetime import timedelta

from classes import *
from readfiles import *
from networkstruc import *
from networkcalc import *
from mpl_toolkits.mplot3d import axes3d
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import time

from mpl_toolkits.mplot3d import Axes3D 

# New axis settings


def rgb_to_hex(r, g, b):
    return '#{:02x}{:02x}{:02x}'.format(r, g, b)

def _format_axes(ax):
    """Visualization options for the 3D axes."""
    # Turn gridlines off
    ax.grid(True,color="k",linestyle="-",linewith=20)
    # Suppress tick labels
    for dim in (ax.xaxis, ax.yaxis, ax.zaxis):
        dim.set_ticks([])
    # Set axes labels
    # ax.set_xlabel("x")
    # ax.set_ylabel("y")
    # ax.set_zlabel("z")


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


def plota_grafo(graph,caminhos,ind_i,i_ind,flagPMU):
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

    
    fixedpos={}
    for no in graph:
        if no.xy != None:
            fixedpos[no.id]=no.xy
    nodePos = nx.spring_layout(g)
    labelList=list(i_ind.values())
    # nx.draw(g,labels=i_ind,with_labels=True,edge_color=colors,width=weights)
    nx.draw_networkx(g,
                node_color =['#74CCF4' for i in nodePos], 
                edge_color=colors,
                pos=nodePos, 
                labels=i_ind,
                node_size=[(len(str(labelList[i]))+1)**2 * 60 for i in nodePos]
                )

    plt.box(False)
    plt.savefig("Graph.pdf")
    plt.savefig("Graph.png")
    
    # make an undirected copy of the digraph
    subgraphs=nx.connected_components(g)

    return list(subgraphs)


def plota_grafo_posSys(graph,caminhos,ind_i,i_ind,flagPMU):
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

    
    fixedpos={}
    for no in graph:
        if no.xy != None:
            fixedpos[no.id]=no.xy
    nodePos = nx.spring_layout(g,pos=fixedpos,fixed=fixedpos.keys(),seed=10)
    
    x=0
    y=0
    for key, pos in fixedpos.items():
        x=x+pos[0]
        y=y+pos[1]
    
    x=x/len(fixedpos)
    y=y/len(fixedpos)
    nodePos[ind_i["gps"]][0]=x
    nodePos[ind_i["gps"]][1]=y


    labelList=list(i_ind.values())
    # nx.draw(g,labels=i_ind,with_labels=True,edge_color=colors,width=weights)
    nx.draw_networkx(g,
                node_color =['#74CCF4' for i in nodePos], 
                edge_color=colors,
                pos=nodePos, 
                labels=i_ind,
                node_size=[(len(str(labelList[i]))+1)**2 * 60 for i in nodePos]
                )

    plt.box(False)
    plt.savefig("Graph_possys.pdf")
    plt.savefig("Graphpossys.png")

    
    # make an undirected copy of the digraph
    subgraphs=nx.connected_components(g)

    return list(subgraphs)



def plota_grafo_3d(graph,caminhos,ind_i,i_ind,flagPMU):
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

    
    fixedpos={}
    for no in graph:
        if no.xy != None:
            fixedpos[no.id]=no.xy+tuple([0])

    x=0
    y=0
    for key, pos in fixedpos.items():
        x=x+pos[0]
        y=y+pos[1]
    
    x=x/len(fixedpos)
    y=y/len(fixedpos)
    
    nodePos = nx.spring_layout(g,dim=3,pos=fixedpos,fixed=fixedpos.keys(),seed=10)

    # nx.draw(g,labels=i_ind,with_labels=True,edge_color=colors,width=weights,pos=nodePos)
    nodePos[ind_i["gps"]][2]=-0.1
    nodePos[ind_i["gps"]][0]=x
    nodePos[ind_i["gps"]][1]=y

    node_xyz = np.array([nodePos[v] for v in sorted(g)])
    edge_xyz = np.array([(nodePos[u], nodePos[v]) for u, v in g.edges()])

    # Create the 3D figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Plot the nodes - alpha is scaled by "depth" automatically
    ax.scatter(*node_xyz.T, s=100, ec="w")

    # Plot the edges
    i=0
    for vizedge in edge_xyz:
        ax.plot(*vizedge.T, color=colors[i])
        i=i+1

    _format_axes(ax)
    fig.tight_layout()
    plt.show()

    # plt.savefig("Graph.pdf")
    # plt.savefig("Graph.png")
    
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
        
def plota_grafo_2d_e_3d(graph,caminhos,ind_i,i_ind,flagPMU):
    g=nx.Graph()

    for no in graph:
        g.add_node(no.id)


    for key, item in caminhos.items():
        cor=list(np.random.choice(range(256), size=3))
        hex_col=rgb_to_hex(cor[0],cor[1],cor[2])
        if len(item)>1:
            for i in range(1,len(item)):
                g.add_edge(item[i-1],item[i],color=hex_col,weight=1)
                
    subgraphs=list(nx.connected_components(g))

    nilhas=len(subgraphs)

    cores=[]
    for _ in range(nilhas):
        cor=list(np.random.choice(range(256), size=3))
        cores.append(rgb_to_hex(cor[0],cor[1],cor[2]))

    edges = g.edges()
    for u,v in edges:
        i=0
        for gr in subgraphs:
            if u in gr:
                g[u][v]['color']=cores[i]
                break
            i=i+1
    for u,v in edges:
        i=0
        if (u == ind_i["gps"] ) |(v == ind_i["gps"] ):
            g[u][v]['weight']=1
        else:
            g[u][v]['weight']=3


    edges = g.edges()
    colors = [g[u][v]['color'] for u,v in edges]
    weights = [g[u][v]['weight'] for u,v in edges]

    
    fixedpos={}
    for no in graph:
        if no.xy != None:
            fixedpos[no.id]=no.xy
    if len(fixedpos)>1:
        nodePos = nx.spring_layout(g,pos=fixedpos,fixed=fixedpos.keys(),seed=10)
    else:
        nodePos = nx.spring_layout(g,k=1/np.sqrt(g.order()),iterations=20,seed=5)
    
    if len(fixedpos)>1:
        x=0
        y=0
        for key, pos in fixedpos.items():
            x=x+pos[0]
            y=y+pos[1]
        
        x=x/len(fixedpos)
        y=y/len(fixedpos)
        nodePos[ind_i["gps"]][0]=x
        nodePos[ind_i["gps"]][1]=y

    fig, ax = plt.subplots(figsize=(16,16))
    labelList=list(i_ind.values())
    # nx.draw(g,labels=i_ind,with_labels=True,edge_color=colors,width=weights)
    
    nx.draw_networkx(g,
                node_color =['#74CCF4' for i in nodePos], 
                edge_color=colors,
                pos=nodePos, 
                labels=i_ind,
                width=weights,
                node_size=[(len(str(labelList[i]))+1)**2*50  for i in nodePos],
                ax=ax
                )

    
    
    plt.savefig("Graph_possys.pdf")
    plt.savefig("Graphpossys.png")
    plt.close()

    
    if len(fixedpos)>1:
        fixedpos={}
        for no in graph:
            if no.xy != None:
                fixedpos[no.id]=no.xy+tuple([0])

        x=0
        y=0
        for key, pos in fixedpos.items():
            x=x+pos[0]
            y=y+pos[1]
        
        x=x/len(fixedpos)
        y=y/len(fixedpos)
        
        nodePos = nx.spring_layout(g,dim=3,pos=fixedpos,fixed=fixedpos.keys(),seed=5)
    else:
        nodePos = nx.spring_layout(g,dim=3,seed=5)
    # nx.draw(g,labels=i_ind,with_labels=True,edge_color=colors,width=weights,pos=nodePos)
    
    if len(fixedpos)>1:
        nodePos[ind_i["gps"]][2]=-0.1
        nodePos[ind_i["gps"]][0]=x
        nodePos[ind_i["gps"]][1]=y

    node_xyz = np.array([nodePos[v] for v in sorted(g)])
    edge_xyz = np.array([(nodePos[u], nodePos[v]) for u, v in g.edges()])

    # Create the 3D figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection="3d")
    ax.view_init(elev=30., azim=-60)

    plt.rcParams['grid.linewidth'] = 4   # change linwidth
    plt.rcParams['grid.color'] = "black" # change color
    ax.set_proj_type('ortho')

    # Plot the nodes - alpha is scaled by "depth" automatically
    
    ax.scatter(*node_xyz.T, s=100, ec="w")
    
    # Plot the edges
    i=0
    for vizedge in edge_xyz:
        ax.plot(*vizedge.T, color=colors[i])
        i=i+1

    _format_axes(ax)
    fig.tight_layout()
    

    plt.savefig("Graph3dsys.pdf")
    plt.savefig("Graph3dsys.png")
    
    # make an undirected copy of the digraph
    subgraphs=nx.connected_components(g)

    return list(subgraphs)


sys = "IEEE118"

dfDBAR, dfDBRAN, dfDMED = read_files(sys)

[bars, nbars, pv, pq, ind_i] = creat_bar(dfDBAR)

[ram, nbran] = create_bran(dfDBRAN, ind_i)



graph = create_graph(bars, ram)

try:
    dfDPOS=pd.read_csv(sys+"/DPOS.csv",header=None)
    dfDPOS.columns=["barra","x","y"]
except:
    dfDPOS=pd.DataFrame()


for idx,row in dfDPOS.iterrows():
    id=row["barra"]
    x=row["x"]
    y=row["y"]
    graph[ind_i[id]].xy=(x,y)



flagPMU=1
if flagPMU==1:
    ind_i["gps"]=len(ind_i)
#%%

i_ind= {v: k for k, v in ind_i.items()}
## filtra DMED


#%%

HT = montaH(graph, dfDMED,ind_i)

Hpseudo,dfPseudo= montaHpseudo(graph,dfDMED,ind_i)




#%%

# Hfat,L,dfPseudo,Permutacao=fatoraH(HT,graph,dfDMED)

Hfat,L,dfPseudoangulo,Permutacao, fluxos_descartaveis=fatoraH_2(HT,graph,dfDMED,Hpseudo,dfPseudo)

print("pseudo fluxo descartavel")
print(dfPseudo.iloc[fluxos_descartaveis])

#%%
print("injeção descartavel")
mask1=dfDMED["de"].isin(dfPseudo.iloc[fluxos_descartaveis]["de"].to_list()+dfPseudo.iloc[fluxos_descartaveis]["para"].to_list())
mask2=dfDMED["type"]==0
print(dfDMED[(mask1) & (mask2)])
#%%
caminhos=encontracaminhos(graph,L)

#%%





subgraphs=plota_grafo_2d_e_3d(graph,caminhos,ind_i,i_ind,flagPMU=1)
verifica_medidas_caminhos(graph,HT,dfDMED,subgraphs)


#%%
ilhas=[]
for subgraph in subgraphs:
    il=list(subgraph)
    ilha=[]
    for barra in il:
        if barra<len(graph):
            ilha.append(graph[barra].bar.id)
        else:
            ilha.append("GPS") 
    ilhas.append(ilha)





#%%
