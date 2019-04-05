# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 17:53:00 2019

@author: User
"""
from sympy import *   # "printing" é um pacote para imprimir bonitinho, 
init_printing(use_latex='mathjax')
x,y, z, a, b, c, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10=symbols("x, y, z,a ,b, c, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10")


def OpElem(tipo,subs,*args):
    """Matrizes de tipo elementar. 
    1) As matrizes de tipo 1: São úteis para permutar duas colunas ou linhas.
    Por exemplo
    MatA = Matrix([[-1, -1, 4, 1], [6, -4, 7, -4], [5, -5, 0, -5], [7, 5, -2, 5]])
    OpElem(1,4,2,4)* MatA          troca as linhas 2 com 4.
    MatA*OpElem(1,4,2,4)           troca as colunas 2 com 4.                   
    
    2)As matrizes de tipo 2: São úteis para multiplicar uma  coluna ou linha por um escalar
    Por exemplo
    OpElem(2,4,-5,4)* MatA          Multiplica a linhaa 4 pelo escalar -5.
    MatA*OpElem(2,4,-5,4)            Multiplica a coluna 4 pelo escalar -5.
    
    3)As matrizes de tipo 2: São úteis para multiplicar uma  coluna ou linha por um escalar
    Por exemplo
    OpElem(2,4,-5,4)* MatA          Multiplica a linhaa 4 pelo escalar -5.
    MatA*OpElem(2,4,-5,4) """
    OpE=eye(args[-1])
    if tipo==3:                                  # Linha"subs" = Linha"subs" + args[1] Linha"args[0]"
        OpE[subs-1,args[0]-1] = args[1]
    elif tipo ==2:                          # Linha"subs" = args[0] Linha"subs" 
        OpE[subs-1,subs-1] = args[0]
    else:
        OpE[args[0]-1,args[0]-1], OpE[subs-1,subs-1], OpE[args[0]-1,subs-1], OpE[subs-1,args[0]-1]= 0, 0, 1, 1            
    return OpE

def RedGaussManual(lista,matriz):
    """Faz manualmente a eliminação de Gauss a variável "lista" tem como entrada uma lista de listas
    lista =[lista0,lista1,...], O resultado é a multiplicacao 
    OpElem(Lk)*...*OpElem(L1)*OpElem(L0) * Matrix *OpElem(C0)*OpElem(C1)*...*OpElem(Cp)
    >>Exemplo:
    MatA = Matrix([[-1, -1, 4, 1], [6, -4, 7, -4], [5, -5, 0, -5], [7, 5, -2, 5]])
    RedGaussManual([[3,1,0,6,4],[3,2,0,5,4],[3,3,0,7,4],[-1,1,3,4],[-3,1,2,-4,4],[3,3,1,-6,4],[3,3,2,8,4],[2,3,-RAT(1)/22,4],[2,2,-RAT(1)/10,4],[2,1,RAT(1)/2,4]],MatA)    """
    foldlist=[matriz]
    for ii in lista:
         if ii[0]>0:
                foldlist.append(OpElem(*ii,matriz.rows)*foldlist[-1])
         else:
                ii[0] = -ii[0]
                foldlist.append(foldlist[-1]*OpElem(*ii,matriz.cols))          
                
    return foldlist

####
####
####
####
####
####    
def oquefazer(matriz):
    
    elmnt=next(((ii,jj) for  ii in range(matriz.rows) for jj in range(ii)  if matriz[ii,jj]!=0),-1)
    if elmnt ==-1:
        print("A matriz está escalonada, pode interpretar as equações ou continuar reduzindo até chegar numa matriz diagonal.")
        elmntSUP=next(((ii,jj) for  ii in range(matriz.rows) for jj in range(ii+1,matriz.rows)  if matriz[ii,jj]!=0),-1)
        if elmntSUP==-1:
            print("A matriz já está em forma Gauss-Jordan, pode interpretar as equações")
        else:
            elmntSUPAux=(elmntSUP[0]+1,elmntSUP[1]+1)
            return print("\n Se quiser zerar a diagonal superior, você pode zerar o elemento ",elmntSUPAux, ". Como?\n Substitua a  linha ",elmntSUPAux[0], 
                 "multiplicando a  linha  ", elmntSUPAux[1], "por ",- matriz[elmntSUP]/matriz[elmntSUP[1],elmntSUP[1]],  "vezes e somando com a linha", elmntSUPAux[0],
                 "isto é;\n linha",elmntSUPAux[0],"-> linha",elmntSUPAux[0],"+ ",- matriz[elmntSUP]/matriz[elmntSUP[1],elmntSUP[1]] ,"x linha",elmntSUPAux[1], 
                     "\n Com qual comando? adicione na lista ",[3,elmntSUPAux[0],elmntSUPAux[1],- matriz[elmntSUP]/matriz[elmntSUP[1],elmntSUP[1]]])
        
    elif  matriz[elmnt[1],elmnt[1]]==0:
        print("para zerar esta coluna precisa permutar duas linhas pois o elemento da diagonal é zero")
    else:
        elemAux=(elmnt[0]+1,elmnt[1]+1)
        return print("você deveria zerar o elemento ",elemAux, ". Como?\n Substitua a  linha ",elemAux[0], 
                 "multiplicando a  linha  ", elemAux[1], "por ", - matriz[elmnt]/matriz[elmnt[1],elmnt[1]],  "vezes e somando com a linha", elemAux[0],
                 "isto é;\n linha",elemAux[0],"-> linha",elemAux[0],"+ ",- matriz[elmnt]/matriz[elmnt[1],elmnt[1]] ,"x linha",elemAux[1], 
                     "\n Com qual comando? adicione na lista ",[3,elemAux[0],elemAux[1],- matriz[elmnt]/matriz[elmnt[1],elmnt[1]]])

####
####
####
####
####
#### 

        

def InterpretarEquacoes(matrizaumentada):
    MATAUG= Matrix(matrizaumentada)
    listainput = input("Insira o nome das variáveis (na  ordem correta) separadas por um espaço:  ")
    varliista=list(map(sympify, listainput.split()))
    SISTEMA=MATAUG[:,:len(varliista)]*Matrix(varliista)
    diferenca=len(varliista)-MATAUG.cols
    param=MATAUG.rank()-len(varliista)
    partenaohomogenea=MATAUG[:,len(varliista):]
    
    if diferenca!=0:
        SISTEMANH=list(map(lambda y: SISTEMA-partenaohomogenea.col(y), range(partenaohomogenea.cols) ))
    else:
        print("Sistema homogêneo sempre tem pelo menos uma solução (0,0,0) ", MATAUG.rank(),"igual ao número de variáveis ", len(varliista), "variáveis, portanto tem "
              , diferenca , "parámetro(s) livre(s)")

    
    if param==0:
        print("A matrix tem posto ", MATAUG.rank(),"igual ao número de variáveis ", len(varliista), "variáveis, portanto tem "
          , param , "parámetro(s) livre(s)")
    else:
            print("A matrix tem posto ", MATAUG.rank(),"e ", len(varliista), "variáveis, portanto tem "
            , param , "parámetro(s) livre(s)")
    print("Sistema: ")
    
    return  SISTEMANH

def SaoLI(Vetores):
    Rango=Matrix(Vetores).rank()
    dimn=Matrix(Vetores).rows
    print("Os vetores são LI se a matriz formada de vetores linha tem posto igual ao número de vetores.\n")
    if dimn!=Rango:
         print("Neste caso escalonando a matriz  [ver comando   MATRIZ.rref() ], pode-se ver que pelo menos uma das linhas é zero; pelo qual:\n \n posto(A)=",Rango,"<",dimn,
              "= número de vetores .\n \n Portanto os vetores NÃO são LI; isto é, um dos vetores pode-se escrever como soma dos outros.\n",
              "Além disso, podemos concluir que a dimensão do espaço gerado é: ", Rango )
    else:
        print("Neste caso escalonando a matriz  [ver comando   MATRIZ.rref() ], pode-se ver que todas as linhas são distintas de zero; pelo qual:\n \n posto(A)=",Rango,
              "=",dimn, "= número de vetores.\n \n Em consequência os vetores são LI; e portanto a dimensão do espaço gerado pelos vetores é: ",  Rango )
        

def EspacoGerado(Vetores):
    Rango=Matrix(Vetores).rank()
    dimn=Matrix(Vetores).rows
    print("O espaço gerado é o conjunto de combinações lineares.\n")
    if dimn!=Rango:
         print(" ", Rango )
    else:
        print("Neste caso escalonando a matriz  [ver comando   MATRIZ.rref() ], pode-se ver que todas as linhas são distintas de zero; pelo qual:\n \n posto(A)=",Rango,
              "=",dimn, "= número de vetores.\n \n Em consequência os vetores são LI; e portanto a dimensão do espaço gerado pelos vetores é: ",  Rango )