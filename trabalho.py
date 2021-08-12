# -*- coding: utf-8 -*-
"""
Trabalho MEF 1

@author: César Eduardo Petersen

@date: %(date)s

Utiliza o método dos elementos finitos para resolver um pórtico

esquema da estrutura:
      D-------E----G
     /|   b2  | b5
    / |       |b3
 b0/  |b1     F
  /   |       |b4
 /    |       |
A     B       C

"""

from pyFFEM.pyFFEM import *

# unidades: kN, m

# define material e seção
mat = Material(2e8, 4e-3, 1.5e-4)    

# define os nós
a=5.2
b=3.6

A = Node(0, 0)
B = Node(a, 0)
C = Node(2*a, 0)
D = Node(a, b)
E = Node(2*a, b)
F = Node(2*a, b/2)
G = Node(2.25*a, b)

# define os apoios
A.apoio = [True, True, False]
B.apoio = [True, True, True]
C.apoio = [True, True, False]

# aplica forças nodais
P = 36
H = 26

D += ForcaNodal(0, -P)
G += ForcaNodal(0, -P)
F += ForcaNodal(H, 0)

# define as barras
b0 = Barra(A, D, mat)
b1 = Barra(B, D, mat)
b2 = Barra(D, E, mat)
#b3 = Barra(F, E, mat)
#b4 = Barra(C, F, mat)
b3 = Barra(E, F, mat)
b4 = Barra(F, C, mat)
b5 = Barra(E, G, mat)

# aplica forças distribuidas
q = 16
b1 += ForcaDistribuida(q)
b2 += ForcaDistribuida(0, -q)

# monta a estrutura
estrutura = Estrutura([A, B, C, D, E, F, G], 
                      [b0, b1, b2, b3, b4, b5])

# exibe a estrutura
#estrutura.exibir()

# calcula
estrutura.calcular()

estrutura.exibirDeformado(750)

#estrutura.diagramaMomento()

#estrutura.diagramaCortante()

#estrutura.diagramaNormal()
