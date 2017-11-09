#!/usr/bin/env python3


import numpy as np

def write_Points(w , l , dx , d_or , cortes_y , i):
# Cria os vetores 'px', 'py' e 'd_or' a partir das entradas:
#     w - largura da partícula_i
#     l - comprimento da partícula_i
#     d - distância entre a partícula i e i-1
#     d_or - vetor com centros de partículas
#     i - índice da partícula atual.
##

    px = np.zeros(4)
    py = np.zeros(4)
# Descreve os pontos em X
    x1 = -w[i]/2
    x2 =  w[i]/2
    x3 =  x2
    x4 =  x1
# Descreve os pontos em Y, sem cortes.
    y1 =  l[i]/2
    y2 =  l[i]/2
    y3 = -l[i]/2
    y4 = -l[i]/2
# Concatena os valores em 'px' e 'py'
    px = [x1 , x2 , x3 , x4]
    py = [y1 , y2 , y3 , y4] - cortes_y

    if i==1:
        d_or[i,:]=[0 , 0 , 0] # primeira partícula na origem
    else:
        d_or[i,:]=[w[i-1]/2 + w[i]/2 + d_or [i-1,0] + dx[i] , 0 , 0] #
    # soma w/2 da particula anterior
    # com w/2 da particula atual
    # a distancia em x da origem acumulada
    # mais a distancia entre as duas partículas 'd'
    # OBS: Atualmente somente em X

    return px , py , d_or
