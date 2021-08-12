# -*- coding: utf-8 -*-
"""
Módulo para definição da estrutura para método dos elementos finitos

@author: César Eduardo Petersen

@date: Tue Jul 20 13:25:50 2021

Fornece as definições para nós, barras e forças aplicadas e distribuidas

"""

import numpy as np
from math import pi, atan2
from .material import *
from .forcas import *

## FUNÇÕES DE APOIO
# Calcula a magnitude de um vetor.  
mag = lambda x: np.sqrt(x.dot(x))

# Normaliza o vetor
norm = lambda x: x/mag(x)

# função sinal
sgn = lambda x: 1 if x>0 else 0 if x==0 else -1

##### FORÇAS #####
        

  
##### GEOMETRIA #####      
  
class Node:
    __id = 0
    def __init__(self, x, y = None):
        """
        Define um nó

        Args:
            x (float): Vetor de coordenadas ou coordenada x.
            y (float, optional): coordenada y.

        Returns:
            None.

        """
        if isinstance(x, (tuple, list)):
            self.coord = np.array(x)
        else:
            self.coord = np.array([x, y])
            
        self.coordDesloc = np.array([0., 0., 0.])
        
        # inicializa os apoios como livres por padrão
        self.apoio =  np.array([False, False, False])
            
        self.forca = np.array([0., 0.])
        self.momento = 0.
        
        self.forcaEqv = np.array([0., 0., 0.])
        
        # incrementador de nó
        self.id = Node.__id
        Node.__id += 1
    
    def __iadd__(self, forca: ForcaNodal):
        """
        Adiciona uma força nodal no nó

        Args:
            forca (ForcaNodal): Força nodal adicionada.

        Returns:
            None.

        """
        if not isinstance(forca, ForcaNodal):
            raise TypeError("Precisa adicionar objeto do tipo ForcaNodal")
        else:            
            self.forca += forca.forca
            self.momento += forca.momento
            
        return self
    
    def plotaApoio(self, plt):
        """
        Plota os apoios no nó

        Args:
            plt (matplotlib.pyplot): gráfico.

        """
        
        # plota apoios
        if all(self.apoio): # engaste
            plt.plot(*self.coord, marker="s",ms=15, color="tab:blue", zorder=1)
        else:
            if self.apoio[0]:
                plt.plot(*self.coord, marker=6,ms=15, color="tab:blue", zorder=1)
            if self.apoio[1]:
                plt.plot(*self.coord, marker=5,ms=15, color="tab:blue", zorder=1)
        

def T_rot2(c, s):
    """
    Retorna a matriz de rotação

    Args:
        c (float): cosseno do angulo.
        s (float): seno do angulo.

    Returns:
        float[][]: Matriz de rotação do elemento.

    """
    return np.array([[ c, s, 0, 0, 0, 0],
                     [-s, c, 0, 0, 0, 0],
                     [ 0, 0, 1, 0, 0, 0],
                     [ 0, 0, 0, c, s, 0],
                     [ 0, 0, 0,-s, c, 0],
                     [ 0, 0, 0, 0, 0, 1]])



class Barra:
    """
    Subclasse para elemento tipo barra de portico
    """
    def __init__(self, no1: Node, no2: Node, mat: Material):
        self.no1 = no1
        self.no2 = no2
        self.forca = None
        
        # cria lista dos graus de liberade de cada nó
        self.gls = [3*no1.id + i for i in range(0,3)] + [3*no2.id + i for i in range(0,3)]
        
        # vetores de coordenadas para plotagem
        self.x = [no1.coord[0], no2.coord[0]]
        self.y = [no1.coord[1], no2.coord[1]]
        
        # definiçao do material
        self.mat = mat
        E = self.mat.E
        A = self.mat.A
        I = self.mat.I              
        
        # definição da matriz de rigidez local
        b = self.no2.coord - self.no1.coord
        L = mag(b)
        self.L = L 
        self.T = T_rot2(b[0]/L, b[1]/L) 
        self.theta = atan2(b[0]/L, b[1]/L)/(pi/180)
        self.K = E*np.array([[ A/L,     0    ,    0    ,-A/L,     0    ,    0    ],
                             [  0 , 12*I/L**3, 6*I/L**2,  0 ,-12*I/L**3, 6*I/L**2],
                             [  0 ,  6*I/L**2, 4*I/L   ,  0 , -6*I/L**2, 2*I/L   ],
                             [-A/L,     0    ,    0    , A/L,     0    ,    0    ],
                             [  0 ,-12*I/L**3,-6*I/L**2,  0 , 12*I/L**3,-6*I/L**2],
                             [  0 ,  6*I/L**2, 2*I/L   ,  0 , -6*I/L**2, 4*I/L   ]])
        
        # funções de forma
        self.N = lambda xi: np.array([
            0,
            1/2-3/4*xi+1/4*xi**3,
            L/8*(1-xi-xi**2+xi**3),
            0,
            1/2+3/4*xi-1/4*xi**3,  
            L/8*(-1-xi+xi**2+xi**3)
        ])
        
    def __iadd__(self, forca: ForcaDistribuida):
        """
        Adiciona uma força distribuida na barra

        Args:
            forca (ForcaDistribuida): Força distribuida.

        Returns:
            None.

        """
        if not isinstance(forca, ForcaDistribuida):
            raise TypeError("Precisa adicionar objeto do tipo ForcaNodal")
        else:            
            self.forca = forca
            fNodais = self.T.transpose() @ forca.nodaisEq(self.L)
            self.no1.forcaEqv += fNodais[0:3]
            self.no2.forcaEqv += fNodais[3:6]
            
        return self
    
    ##### ESFORÇOS #####
        
    def _uel(self):
        """
        Retorna os deslocamentos nodais no sistema local
        """
        return self.T @ np.concatenate((self.no1.coordDesloc, self.no2.coordDesloc))
        
    
    def Normal(self, xi):
        """
        Calcula o esforço normal no espaço ξ [-1,1]
        """
        # deslocamentos nodais locais
        return self.mat.EA*(self._uel()[3]-self._uel()[0])/self.L
    
    
    def Momento(self, xi):
        """
        Calcula o momento fletor no espaço ξ [-1,1]
        """
        # função de forma derivada 2x para momento
        B_xx = np.array([
            0, 
            3*xi/2, 
            self.L*(-2 +6*xi)/8,
            0,
            -3*xi/2,
            self.L*( 2 +6*xi)/8
        ])
                
        return 4*self.mat.EI*(B_xx @ self._uel())/self.L**2
            
           
    def Cortante(self, xi):
        """
        Calcula o esforço cortante no espaço ξ [-1,1]
        """
        # função de forma derivada 3x para cortante
        B_xxx = np.array([
            0, 
            3/2, 
            self.L*3/4,
            0,
            -3/2,
            self.L*3/4
        ])
        
            
        return 8*self.mat.EI*(B_xxx @ self._uel())/self.L**3
    
                    