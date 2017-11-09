#!/usr/bin/env python3
import numpy as np

def compute_Happ(steps , s , Ms):
# Describes the applied field
# following the logic used by OOMMF
# N : é o numero de passos de simulação
#          | Inicial  |   Final   |Duracao|
#          |x_i y   z | x_f y   z | steps |
#          |1   2   3 | 4   5   6 |  7    |
# s=    [   0   0   0   1   0   0   221
#           1   0   0   1   1   0   221
#           1   1   0   0   1   0   221
#           0   1   0   0   0   0   221
#           0   0   0   0   0   0   221]
    if sum(s[ : , 6]) != steps:
        print("Total pulse duration != Total steps")

    h_app = np.zeros((steps,3))
    count = 0

    t2am = 1/(4*4*np.arctan(1)*10**-7)
    t2am = 1
    print(t2am)

    for line in range(len(s)): #iterate over lines of s
        for dimension in range(3):  # iterate over de dimension x, y, z
            slope = (s[line,dimension+3]-s[line, dimension])/(s[line,6]-1)
            for n in range(int(s[line,6])):  # iterate over duration steps of the ramp
                h_app[count+n,dimension] = (s[line,dimension]+n*slope)*t2am
        count = count + n+1

    return h_app
