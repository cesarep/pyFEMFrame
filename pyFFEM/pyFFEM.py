# -*- coding: utf-8 -*-
"""
Módulo para cálculo de Pórticos em Elementos Finitos

@author: César Eduardo Petersen

@date: Tue Jul 20 13:25:50 2021

Descrição

"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg

from .elementos import *

def deslocCurva(ptos, desloc, R):
    """
    Retorna os pontos de uma curva rotacionados e deslocados para plotagem

    Args:
        ptos (float[2][]): Conjunto de pontos x, y que descrevem a curva.
        desloc (float[2]): vetor de deslocamento da curva.
        R (float[2][2]): Matriz de rotação.

    Returns:
        ptos (float[2][]): Pontos x, y deslocados e rotacionados.
    """
    # rotaciona os pontos  
    ptos = R[0:2, 0:2].transpose() @ ptos
    
    # desloca os pontos
    ptos += desloc.reshape((2, 1))
    
    return ptos

class Estrutura:
    """
    Classe para descrever uma estrutura, agrupando nós e elementos
    """
    def __init__(self, nos, elems):
        self.nos = nos
        self.elementos = elems
        
    def __iadd__(self, obj):
        """
        Adiciona um nó ou Elemento a estrutura

        Args:
            obj (Node|Barra): Nó ou Elemento a adicionar.

        Returns:
            None.

        """
        if isinstance(obj, Node):
            self.nos.append(obj)
            
        elif isinstance(obj, Barra):
            self.elementos.append(obj)
            
        else:
            raise TypeError("Precisa adicionar objeto do tipo Node ou Barra")
        
        return self
            
            
    def calcular(self):
        """
        Calcula a estrutura.

        Returns:
            None.

        """
        # 3 graus de liberdade por nó
        self.ngl = len(self.nos)*3 
        # prepara a matriz de rigidez global
        self.Kg = np.zeros((self.ngl, self.ngl))
        self.F = np.zeros((self.ngl))
        self.U = np.zeros((self.ngl))
        
        # itera cada elemento para montar a matriz de rigidez global   
        for elem in self.elementos:
            # transforma a matriz local para sistema global
            K = elem.T.transpose() @ elem.K @ elem.T

            # adiciona a matriz do elemento na matriz global
            self.Kg[np.ix_(elem.gls, elem.gls)] += K
            
        # itera cada nó para montar o vetor F e separar os graus de liberdade livres e fixos
        glF = []
        glL = []
        for no in self.nos:
            glF += [3*no.id + i for i in range(0,3) if no.apoio[i]]
            glL += [3*no.id + i for i in range(0,3) if not no.apoio[i]]
            self.F[3*no.id:3*no.id+2] = no.forca + no.forcaEqv[0:2]
            self.F[3*no.id+2] = no.momento + no.forcaEqv[2]
            
        self.glF = glF
        self.glL = glL
            
        # Separando as submatrizes
        KFF = self.Kg[np.ix_(glF, glF)]
        KFL = self.Kg[np.ix_(glF, glL)]
        KLL = self.Kg[np.ix_(glL, glL)]
        
        FL = self.F[np.ix_(glL)]
        
        # calcula os deslocamentos livres
        self.U[np.ix_(glL)] = linalg.solve(KLL, FL)
        
        # calcula as reações de apoio
        self.F[np.ix_(glF)] = KFL @ self.U[np.ix_(glL)] - self.F[np.ix_(glF)]
        
        # atualiza a infomação dos desloc
        for no in self.nos:
            no.coordDesloc = self.U[3*no.id:3*no.id+3]
        
        with np.printoptions(formatter={'float': lambda x: "% 8.2f"%x}):
            print("Reações de apoio:")
            print(self.F[np.ix_(glF)].reshape((len(glF),1)))
            
        print("\nDeslocamentos:")
        print(self.U[np.ix_(glL)].reshape((len(glL),1)))
            
      
    ##### PARTE GRÁFICA #####
    
    def exibir(self):
        """
        Plota a imagem da estrutura.

        Returns:
            None.

        """
        plt.title("Esquema da Estrutura")
        
        # plota os nós
        plt.scatter([x.coord[0] for x in self.nos], [x.coord[1] for x in self.nos], c='k', zorder=0)
        
        # plota as forças nos nós
        for no in self.nos:
            nox, noy = no.coord
           
            # se há forças
            if (no.forca != 0).any():
                
                # coordenada deslocada da força para plotagem
                fcoord = no.coord - 1.5*norm(no.forca)
                
                # se há força horizontal
                if no.forca[0] != 0: 
                    plt.arrow(fcoord[0], noy, norm(no.forca)[0], 0, width=.05, color="r")
                    plt.annotate("%.1f"%no.forca[0], (fcoord[0], noy), 
                                 textcoords='offset points', xytext=(-15*sgn(no.forca[0]),0), ha='center', va='center')
                
                # se há força vertical
                if no.forca[1] != 0: 
                    plt.arrow(nox, fcoord[1], 0, norm(no.forca)[1], width=.05, color="r")
                    plt.annotate("%.1f"%no.forca[1], (nox, fcoord[1]), 
                                 textcoords='offset points', xytext=(0, -10*sgn(no.forca[1])), ha='center', va='center')
            
            # se há momento
            if no.momento != 0:
                plt.plot(nox, noy, marker=r'$\circlearrow%s$'%('left' if no.momento > 0 else 'right'),ms=25, color="r")
                plt.annotate("%.1f"%no.momento, no.coord, textcoords='offset points', xytext=(15, 15), ha='center', va='center')

            
            # plota apoios
            no.plotaApoio(plt)
                
                
        # plota as barras
        for barra in self.elementos:
            plt.plot(barra.x, barra.y, 'k', zorder=0)
            no1 = barra.no1.coord
            no2 = barra.no2.coord
            # plot as forças distribuidas na barra, caso existam
            if isinstance(barra.forca, ForcaDistribuida):
                qxy = deslocCurva(barra.forca.plotForca(barra), barra.no1.coord, barra.T)
                q1 = qxy[:, 0]
                q2 = qxy[:,-1]
                d1 = no1 - q1
                d2 = no2 - q2
                plt.plot(qxy[0], qxy[1], 'b')
                plt.arrow(*q1, *.5*d1, width=.05, color="b", zorder=2)
                plt.arrow(*q2, *.5*d2, width=.05, color="b", zorder=2)
                plt.annotate("%.1f"%barra.forca.q1, q1, textcoords='offset points', xytext=-15*d1, ha='center', va='center')
                plt.annotate("%.1f"%barra.forca.q2, q2, textcoords='offset points', xytext=-15*d2, ha='center', va='center')
        
        # deixa grafico num aspecto 1:1
        plt.axis("equal") 
        plt.show()            
        
            
    def exibirDeformado(self, escala=250):
        """
        Plota a deformada da estrutura
        
        Args:
            escala (float): Escala de distorção da deformada, padrão 250.

        Returns:
            None.

        """
        plt.title("Deformada da Estrutura")
        # plota os nós nas coordenadas deformadas
        # plt.scatter([x.coord[0] + escala*x.coordDesloc[0] for x in self.nos], [x.coord[1] + escala*x.coordDesloc[1] for x in self.nos], c='k', zorder=0)
        
        # plota as forças nos nós
        for no in self.nos:
            # plota apoios
            no.plotaApoio(plt)
               
            # plota reações
            nox, noy = no.coord
            forca = self.F[3*no.id:3*no.id+2]
            fcoord = no.coord - 1.5*norm(forca)
            
            if(no.apoio[0]):
                plt.arrow(fcoord[0], noy, norm(forca)[0], 0, width=.05, color="r", zorder=2)
                plt.annotate("%.2f"%forca[0], (fcoord[0], noy), 
                             textcoords='offset points', xytext=(-15*sgn(forca[0]),0), ha='center', va='center')
                
            if(no.apoio[1]):
                plt.arrow(nox, fcoord[1], 0, norm(forca)[1], width=.05, color="r", zorder=2)
                plt.annotate("%.2f"%forca[1], (nox, fcoord[1]), 
                             textcoords='offset points', xytext=(0, -15*sgn(forca[1])), ha='center', va='center')
                
            if(no.apoio[2]):
                momento = self.F[3*no.id+2]
                plt.plot(nox, noy, marker=r'$\circlearrow%s$'%('left' if momento > 0 else 'right'),ms=25, color="r", zorder=2)
                plt.annotate("%.2f"%momento, no.coord, textcoords='offset points', xytext=(15, 15), ha='center', va='center')
                    
                
                
        # plota as barras
        xi = np.linspace(-1, 1)
        for barra in self.elementos:
            plt.plot(barra.x, barra.y, 'k--', linewidth=1, zorder=0)
            ## constroi a deformada da barra
            
            # desconta a deformação axial da barra
            x = np.linspace(escala*barra._uel()[0], barra.L + escala*barra._uel()[3])

            # calcula a deformada ao longo da barra
            y = np.array([barra.N(i) @ barra._uel() for i in xi])
            
            # rotaciona e desloca a curva para a barra
            xy = deslocCurva(np.array([x, escala*y]), barra.no1.coord, barra.T)
            
            plt.plot(*xy, 'k', zorder=0)
        
        # deixa grafico num aspecto 1:1
        plt.axis("equal") 
        plt.show()
        
        
    def _sentido(self, t):
        """
        retorna +1 ou -1 para inverter ou nao o sentido do diagrama
        dependendo do angulo da barra, entre -45 e 135 graus nao inverte, restante
        inverte.

        Args:
            t (float): angulo da barra.

        Returns:
            x (int): +1 ou -1.

        """
        if t <= 135:
            return 1
        elif t <= 315:
            return -1
        else:
            return 1

        
    def diagramaMomento(self, escala=20):
        """
        Plota o diagrama de momento fletor da estrutura
        
        Args:
            escala (float): Escala de distorção do diagrama, padrão 20.

        Returns:
            None.

        """
        plt.title("Diagrama de Momento")
        # plota os nós nas coordenadas
        plt.scatter([x.coord[0] for x in self.nos], [x.coord[1] for x in self.nos], c='k', zorder=0)
        
        # plota as forças nos nós
        for no in self.nos:
            no.plotaApoio(plt)

                
        # plota as barras
        xi = np.linspace(-1, 1)
        for barra in self.elementos:
            plt.plot(barra.x, barra.y, 'k', linewidth=1, zorder=0)
            
            # monta o diagrama de momento, inclui os limitantes 2 vezes para fechar a curva
            x = np.concatenate([[0], np.linspace(0, barra.L), [barra.L]])
            
            # inclui zeros no final para "fechar" a curva
            y = np.array([0]+[barra.Momento(i) for i in xi]+[0])
                        
            # rotaciona e desloca a curva para a barra e inverte o sentido do diagrama de momento
            xy = deslocCurva(np.array([x, -y/escala], dtype=np.ndarray), barra.no1.coord, barra.T)
            plt.plot(*xy, 'b', zorder=0)
            
            plt.annotate("%.2f"%abs(y[1]), xy.transpose()[1], textcoords='offset points', xytext=(5, 5), ha='center', va='center')
            plt.annotate("%.2f"%abs(y[-2]), xy.transpose()[-2], textcoords='offset points', xytext=(5, 5), ha='center', va='center')
            
            
        
        # deixa grafico num aspecto 1:1
        plt.axis("equal") 
        plt.show()  
        
        
        
    def diagramaCortante(self, escala=20):
        """
        Plota o diagrama de Cortante da estrutura
        
        Args:
            escala (float): Escala de distorção do diagrama, padrão 20.

        Returns:
            None.

        """
        plt.title("Diagrama de Cortante")
        # plota os nós nas coordenadas
        plt.scatter([x.coord[0] for x in self.nos], [x.coord[1] for x in self.nos], c='k', zorder=0)
        
        # plota as forças nos nós
        for no in self.nos:
            no.plotaApoio(plt)

                
        # plota as barras
        xi = np.linspace(-1, 1)
        for barra in self.elementos:
            plt.plot(barra.x, barra.y, 'k', linewidth=1, zorder=0)
            
            # monta o diagrama de cortante, inclui os limitantes 2 vezes para fechar a curva
            x = np.concatenate([[0], np.linspace(0, barra.L), [barra.L]])
            
            # inclui zeros no final para "fechar" a curva
            y = np.array([0]+[barra.Cortante(i) for i in xi]+[0])
            
            # define o sentido da barra para inverter ou nao o sentido do diagrama
            invt = self._sentido(barra.theta)
            
            # rotaciona e desloca a curva para a barra, 
            xy = deslocCurva(np.array([x, invt*y/escala], dtype=np.ndarray), barra.no1.coord, barra.T)
            plt.plot(*xy, 'b', zorder=0)
            
            plt.annotate("%.2f"%y[25], xy.transpose()[25], textcoords='offset points', xytext=(10, 10), ha='center', va='center')

        # deixa grafico num aspecto 1:1
        plt.axis("equal") 
        plt.show()
        
    def diagramaNormal(self, escala=20):
        """
        Plota o diagrama de esforço Normal
        
        Args:
            escala (float): Escala de distorção do diagrama, padrão 20.

        Returns:
            None.

        """
        plt.title("Diagrama de Normal")
        # plota os nós nas coordenadas
        plt.scatter([x.coord[0] for x in self.nos], [x.coord[1] for x in self.nos], c='k', zorder=0)
        
        # plota as forças nos nós
        for no in self.nos:
            no.plotaApoio(plt)

                
        # plota as barras
        xi = np.linspace(-1, 1)
        for barra in self.elementos:
            plt.plot(barra.x, barra.y, 'k', linewidth=1, zorder=0)
            
            # monta o diagrama de cortante, inclui os limitantes 2 vezes para fechar a curva
            x = np.concatenate([[0], np.linspace(0, barra.L), [barra.L]])
            
            # inclui zeros no final para "fechar" a curva
            y = np.array([0]+[barra.Normal(i) for i in xi]+[0])
            
            # define o sentido da barra para inverter ou nao o sentido do diagrama
            invt = self._sentido(barra.theta)
            
            # rotaciona e desloca a curva para a barra, 
            xy = deslocCurva(np.array([x, invt*y/escala], dtype=np.ndarray), barra.no1.coord, barra.T)
            plt.plot(*xy, 'b', zorder=0)
            
            plt.annotate("%.2f"%y[25], xy.transpose()[25], textcoords='offset points', xytext=(10, 10), ha='center', va='center')

        # deixa grafico num aspecto 1:1
        plt.axis("equal") 
        plt.show()
            
        
        

        
            
        
        
    