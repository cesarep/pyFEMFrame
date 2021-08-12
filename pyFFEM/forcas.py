# -*- coding: utf-8 -*-
"""
Nome

@author: cesar

@date: Wed Aug 11 18:06:11 2021

Descrição

"""

import numpy as np

class ForcaNodal:
    def __init__(self, x, y = None, m = 0):
        """
        Define uma força no nó

        Args:
            x (float): Vetor de componentes ou componente x da força.
            y (float, optional): componente y da força.
            m (float, optional): momento aplicado. Defaults to 0.

        Returns:
            None.

        """
        if isinstance(x, (tuple, list)):
            self.forca = np.array(x)
        else:
            self.forca = np.array([x, y])
        
        self.momento = m
        

class ForcaDistribuida:
    def __init__(self, q1, q2 = None):
        """
        Define uma força distribuida retangular ou triangular

        Args:
            q1 (float): valor da força, ou altura inicial.
            q2 (float, optional): altura final da força. Defaults to None.
            
        Returns:
            None.

        """
        # define as cargas retangulares se os segundos parametros nao forem definidos
        self.q1 = q1      
        self.q2 = q1 if q2 is None else q2
        self.qm = (self.q2+self.q1)/2
        self.dq = self.q2-self.q1
        
        # função da carga no espaço ξ [-1,1]
        self.qx = lambda x: self.qm + self.dq*x/2
                
    def plotForca(self, barra):
        """
        retorna o conjunto de pontos para plotar a força

        Args:
            barra (Barra): Objeto de barra.

        Returns:
            float[]: matriz de coordenadas rotacionados.

        """
        L = barra.L
        
        # função da carga no espaço x [0,L] reescalada para a maior carga = 1  
        x = np.linspace(0, L)
        y = -self.qx(2*x/L-1)/max(abs(self.q1), abs(self.q2))
        
        return np.array([x, y])
                  
        
    def nodaisEq(self, L):
        """
        

        Args:
            L (float): comprimento da barra.

        Returns:
            float[]: vetor de forças nodais equivalentes.

        """
        qm = self.qm
        dq = self.dq
        return np.array([0, qm-dq/5, qm*L/6 - dq*L/60, 0, qm+dq/5, -qm*L/6 - dq*L/60])*L/2
        