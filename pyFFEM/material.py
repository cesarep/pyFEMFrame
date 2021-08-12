# -*- coding: utf-8 -*-
"""
Nome

@author: cesar

@date: Wed Aug 11 18:05:02 2021

Descrição

"""

class Material:
    def __init__(self, E, A, I):
        """
        Define um material elastico-linear

        Args:
            E (float): Módulo de Elasticidade.
            A (float): Área da seção transversal.
            I (float): Momento de Inercia da seção.

        Returns:
            None.

        """
        self.E = E
        self.A = A
        self.I = I
        
        self.EI = E*I
        self.EA = E*A