from classes import *
from readfiles import *
from networkstruc import *
import numpy as np




tolpiv=1e-9


def reorder_columns(matrix, permutation_vector):
    # Convert matrix and permutation vector to numpy arrays
    matrix = np.array(matrix)
    permutation_vector = np.array(permutation_vector)
    
    # Get the number of columns in the matrix
    num_columns = matrix.shape[1]
    
    # Check if permutation vector is valid
    if len(permutation_vector) != num_columns:
        raise ValueError("Permutation vector length does not match the number of columns in the matrix.")
    
    # Reorder the columns of the matrix
    reordered_matrix = matrix[:, permutation_vector]
    
    return reordered_matrix

def montaH(graph,dfDMED,ind_i):

    """
    Função para montar a matriz H
    @param: graph lista de objetos do tipo nó com todas as indormações sobre a rede
    @param: DMEDS dataframe com as informações de todas as medidas possiveis
    @return: HT Matriz Jacobiana transposta
    @return: Oorder vetor com a ordenação das medidas   
    """
    
    i=0
    flagPMUV=sum(dfDMED["type"]==5)>0
            
    HT=np.zeros((len(graph)+flagPMUV,len(dfDMED)+len(graph)+flagPMUV))

    dfDMED["i_de"]=int(-9)
    dfDMED["i_para"]=int(-9)

    nmeds=len(dfDMED)

    for idx,med in dfDMED.iterrows():
        dfDMED.at[idx,"i_de"]=int(ind_i[med["de"]])
        if ~((med["type"]==0) | (med["type"]==5)):
            dfDMED.at[idx,"i_para"]=int(ind_i[med["para"]])


    

    dfDMED=dfDMED.astype({"i_de":'int32',"i_para":'int32'},)
    dfDMED["inst"]=1
    for idx,med in dfDMED.iterrows():
        # if the measurement is a Flow or a current fwo
        if (med.type==2)| (med.type==6):
            de=int(med.i_de)
            para=int(med.i_para)
            if para in  graph[de].ladjk+graph[de].ladjm:
                HT[de][i]=1
                HT[para][i]=-1
        # if the measurement is a injection 
        elif (med.type==0):
            k=int(med.i_de)
            m=0
            for item in graph[k].ladjk:
                HT[item][i]=-1
                m=m+1
            for item in graph[k].ladjm:
                HT[item][i]=-1
                m=m+1
            HT[k][i]=m
        elif (med.type==5):
            k=int(med.i_de)
            HT[k][i]=1
            HT[-1][i]=-1
        i=i+1

    
    return HT

def montaHpseudo(graph,dfDMED,ind_i):

    """
    Função para montar a matriz H
    @param: graph lista de objetos do tipo nó com todas as indormações sobre a rede
    @param: DMEDS dataframe com as informações de todas as medidas possiveis
    @return: HT Matriz Jacobiana transposta
    @return: Oorder vetor com a ordenação das medidas   
    """
    


    d={}
    d["type"]=[]
    d["de"]=[]
    d["para"]=[]
    d["zmed"]=[]
    d["prec"]=[]
    d["i_de"]=[]
    d["i_para"]=[]

    for index, row in dfDMED[dfDMED["type"]==0].iterrows():
        k=int(row["i_de"])
        for key, item in graph[k].adjk.items():
            d["type"].append(2)
            d["de"].append(graph[item.de].bar.id)
            d["para"].append(graph[item.para].bar.id)
            d["zmed"].append(0)
            d["prec"].append(0)
            d["i_de"].append(item.de)
            d["i_para"].append(item.para)
        for key, item in graph[k].adjm.items():
            d["type"].append(2)
            d["de"].append(graph[item.para].bar.id)
            d["para"].append(graph[item.de].bar.id)
            d["zmed"].append(0)
            d["prec"].append(0)
            d["i_de"].append(item.para)
            d["i_para"].append(item.de)
    
    dfPseudoMedidas=pd.DataFrame(data=d)



    for index, row in dfDMED[(dfDMED["type"]==6)|(dfDMED["type"]==2)].iterrows():
        listind=dfPseudoMedidas[(dfPseudoMedidas["de"]==row["de"]) & (dfPseudoMedidas["para"]==row["para"])].index
        dfPseudoMedidas=dfPseudoMedidas.drop(listind)
    dfPseudoMedidas.reset_index(inplace=True)

    i=0
    flagPMUV=sum(dfDMED["type"]==5)>0
            
    Hpseu=np.zeros((len(graph)+flagPMUV,len(dfPseudoMedidas)))
    







    dfPseudoMedidas=dfPseudoMedidas.astype({"i_de":'int32',"i_para":'int32'})
    dfPseudoMedidas["inst"]=1
    for idx,med in dfPseudoMedidas.iterrows():
        # if the measurement is a Flow or a current fwo
        if (med.type==2)| (med.type==6):
            de=int(med.i_de)
            para=int(med.i_para)
            Hpseu[de][i]=1
            Hpseu[para][i]=-1
        # if the measurement is a injection 
        elif (med.type==0):
            k=int(med.i_de)
            m=0
            for item in graph[k].ladjk:
                Hpseu[item][i]=-1
                m=m+1
            for item in graph[k].ladjm:
                Hpseu[item][i]=-1
                m=m+1
            Hpseu[k][i]=m
        elif (med.type==5):
            k=int(med.i_de)
            Hpseu[k][i]=1
            Hpseu[-1][i]=-1
        i=i+1


    
    return Hpseu,dfPseudoMedidas

def fatoraH(HT,graph,dfDMED):
    """
    Função para fatorar a matriz H
    @param: matriz jacobiana HT array by dimensional de float do numpy
    @param: graph lista de objetos do tipo nó com todas as indormações sobre a rede
    @param: Oorder vetor com a ordenação das medidas
    @param: ind objeto do tipo inividuo com o plano de medições e os indicadores de fit
    @return: H_fat matriz fatorada
    @return: Numero de medidas criticas
    @return: N medidas adicionadas 
    melhorar, na permutacao pode mandar uma medida instalada lá para trás 
    """

    #%%-------------------------------------------inicializa variáveis-----------------------%%#


    tollfill=1e-10

    Hfat=HT.copy()


    nmeds=len(dfDMED)
    flag=sum(dfDMED["type"]==5)>0
    Permutacao=np.array(dfDMED.index,dtype=int) # vetor de permutacoes tem a mesma ordem da Htriang final e em suas poiscoes tem a posicao da medida no plano
    
    L=np.identity((len(graph)+1))*-9999

    dfPseudo=pd.DataFrame()
    
    
    # O vetor Permutacao permite acessar de acordo com a coluna i, a medida referente no plano original
    # ele traduz a ordem das colunas com a ordem do plano data frame
    for  i in range(len(graph)-1+flag):
     # Verifica se o pivo é nulo

        if np.abs(Hfat[i][i])<tolpiv:
            # se for ele procura a próxima medida com pivo não nulo dentro das medidas
            obs=0
            for j in range(i+1,nmeds):
                if np.abs(Hfat[i][j])>tolpiv: # se o pivo é não nulo
                    #encontrou medida dentro do conjunto que existe
                    permutaMedida(Hfat,Permutacao,i,j) # permuta a medida i com j, j dá a informcao
                    obs=1
                    break
            if obs==0:
                # não encontrou, precisa instalar pseudo medida de angulo
                barra=i
                Permutacao=np.append(Permutacao,nmeds)
                Hfat[i][nmeds]=1
                d=pd.DataFrame({"index":[nmeds],"type":[7],"i_de":[barra]})
                dfPseudo=pd.concat([dfPseudo,d])
                permutaMedida(Hfat,Permutacao,i,nmeds)
                nmeds=nmeds+1                
                    
        
        for j in range(i+1,len(graph)+flag):
                if np.abs(Hfat[j][i])>tollfill:
                    L[j,i]=-(Hfat[j][i]/Hfat[i,i])
                    Hfat[j,:]=Hfat[j,:]-(Hfat[j][i]/Hfat[i,i])*(Hfat[i,:])    
                else:  
                    Hfat[j][i]=0  

    # #termina de fazer a Hdelta
    # for i in range(1,(len(graph)-1+flag)):
    #     for j in range(0,i):
    #         if np.abs(Hfat[j][i])>tollfill:
    #             Hfat[j,:] = Hfat[j,:]-Hfat[j,i]*Hfat[i,:]
    #         else:
    #             Hfat[j][i]=0




    
    return Hfat,L,dfPseudo,Permutacao

def fatoraH_2(HT,graph,dfDMED,Hpseudo,dfPseudos):
    """
    Função para fatorar a matriz H
    @param: matriz jacobiana HT array by dimensional de float do numpy
    @param: graph lista de objetos do tipo nó com todas as indormações sobre a rede
    @param: Oorder vetor com a ordenação das medidas
    @param: ind objeto do tipo inividuo com o plano de medições e os indicadores de fit
    @return: H_fat matriz fatorada
    @return: Numero de medidas criticas
    @return: N medidas adicionadas 
    melhorar, na permutacao pode mandar uma medida instalada lá para trás 
    """

    #%%-------------------------------------------inicializa variáveis-----------------------%%#


    tollfill=1e-10

    Hfat=HT.copy()


    nmeds=len(dfDMED)
    flag=sum(dfDMED["type"]==5)>0
    Permutacao=np.array(dfDMED.index,dtype=int) # vetor de permutacoes tem a mesma ordem da Htriang final e em suas poiscoes tem a posicao da medida no plano
    
    L=np.identity((len(graph)+1))*-9999

    dfPseudo=pd.DataFrame()
    
    
    # O vetor Permutacao permite acessar de acordo com a coluna i, a medida referente no plano original
    # ele traduz a ordem das colunas com a ordem do plano data frame
    for  i in range(len(graph)-1+flag):
     # Verifica se o pivo é nulo

        if np.abs(Hfat[i][i])<tolpiv:
            # se for ele procura a próxima medida com pivo não nulo dentro das medidas
            obs=0
            for j in range(i+1,nmeds):
                if np.abs(Hfat[i][j])>tolpiv: # se o pivo é não nulo
                    #encontrou medida dentro do conjunto que existe
                    permutaMedida(Hfat,Permutacao,i,j) # permuta a medida i com j, j dá a informcao
                    obs=1
                    break
            if obs==0:
                # não encontrou, precisa instalar pseudo medida de angulo
                fluxosdescartaveis=[]
                for j in range(Hpseudo.shape[1]):
                    if np.abs(Hpseudo[i][j])>tollfill:
                        fluxosdescartaveis.append(j)
                if len(fluxosdescartaveis)>0:
                    print("existem medidas de fluxo descartaveis")
                    return Hfat,L,dfPseudo,Permutacao,fluxosdescartaveis      

                barra=i
                Permutacao=np.append(Permutacao,nmeds)
                Hfat[i][nmeds]=1
                d=pd.DataFrame({"index":[nmeds],"type":[7],"i_de":[barra]})
                dfPseudo=pd.concat([dfPseudo,d])
                permutaMedida(Hfat,Permutacao,i,nmeds)
                nmeds=nmeds+1                
                    
        
        for j in range(i+1,len(graph)+flag):
                if np.abs(Hfat[j][i])>tollfill:
                    L[j,i]=-(Hfat[j][i]/Hfat[i,i])
                    Hfat[j,:]=Hfat[j,:]-(Hfat[j][i]/Hfat[i,i])*(Hfat[i,:])
                    Hpseudo[j,:]=Hpseudo[j,:]-(Hfat[j][i]/Hfat[i,i])*Hpseudo[j,:]   
                else:  
                    Hfat[j][i]=0  

    # #termina de fazer a Hdelta
    # for i in range(1,(len(graph)-1+flag)):
    #     for j in range(0,i):
    #         if np.abs(Hfat[j][i])>tollfill:
    #             Hfat[j,:] = Hfat[j,:]-Hfat[j,i]*Hfat[i,:]
    #         else:
    #             Hfat[j][i]=0




    
    return Hfat,L,dfPseudo,Permutacao,fluxosdescartaveis





def permutaMedida(Hfat,Permutacao,x,y):
    """
    Troca posicao de x com y na matriz e no vetor 
    """
    
    Hfat[:,[x,y]]=Hfat[:,[y,x]]# altera a ordem da Jacobiana
    aux=Permutacao[y]  
    Permutacao[y]=Permutacao[x] # coloca ela no fim do arquivo de candidatas
    Permutacao[x]=aux # coloca a primeira PMU candidata pro lugar dela
