#!/usr/bin/env python3

import numpy as np              # Lib for math computation and vectorial calculus
import matplotlib.pyplot as plt # Lib to plot results
import subprocess
import time
from rk4 import rk4               # Range-Kuta for solving sLLGS
from dm import dm                 # Function to calculate 'dm'
from write_Points import write_Points # Converts w l and d as px py and d_or
from compute_Happ import compute_Happ
from demag3D3_optimized import demag3D3


#=========================== Constants ================================
gama  = 2.211*10**5 			# Giromagnetic Ratio in 	[m/A.s]*
q     = 1.60217662*10**-19 	    # Electron charge 			[C]
pi    = 4*np.arctan(1)          # Pi 3.14
mu0   = 4*pi*10**-7             # [H/m] or 				    [T.m/A]^
kb    = 1.38064852*10**-23      # Boltzman"s Constant 		[m2.kg/s2.K]
hbar  = 1.0545718*10**-34 	    # Norm. Planks Const. 		[J.s/rad]
t2am  = 1/mu0	    		    # Convert Tesla [T] to A/m
j2ev  = 1/q				        # Convert Joule to eV
#* <cmath> already have a fun. called gamma with 2 M"s
#^ Pi value was calculated using atan(1)*4
#======================================================================



def main():
    # Simulation parameter
    n_of_particles = 2
    steps = 20
    nmc = 1000000   # number of samples used by the Monte Carlo method (demag3D3)

    # Particles Parameters
    Ms = 8E5
    m = np.zeros((steps,n_of_particles,3))
    s = np.zeros((n_of_particles,4,7))
    m[0 , : , : ] = [0 , 1 , 0]

    px = np.zeros((n_of_particles,4))
    py = np.zeros((n_of_particles,4))
    th = np.ones(n_of_particles)*10
    d_or = np.zeros((n_of_particles,3))

    s[:,:,:] = np.array([
        [0, 0, 0, 2, 0, 0, 5],
        [2, 0, 0, 1, 0, 0, 5],
        [1, 0, 0, 0, 0, 0, 5],
        [0, 0, 0, 0, 0, 0, 5]])
    #print(s[0,0])

    h_app = np.zeros((steps,n_of_particles,3))
    temp = np.zeros((steps,4))

    print(len(s[0,:,:]))
    h_app[:,0,:] = compute_Happ(steps, s[0,:,:], Ms)

    print(h_app[:,0,:])

    w = [50,50]
    l = [100,100]
    dx = [10,10]
    cortes_y = np.zeros(4)

    for j in range(n_of_particles):
        px[j,:], py[j,:], d_or = write_Points(w,l,dx,d_or,cortes_y,j)
    #p = subprocess.run(".\demag3D3.exe", stdout = subprocess.PIPE, shell = True, encoding="cp850").stdout
    #print(p)

    Nd = np.zeros(9)
    t = time.time()
    print("px: ", px[0,:], "py: ", py[0,:],"th: ", th[0],"\n")
    Nd = demag3D3(px[0,:],py[0,:],th[0], nmc)
    print("Elapsed time: ",time.time()-t,"s")

    # t = time.time()
    # print("px: ", px[0,:], "py: ", py[0,:],"th: ", th[0],"\n")
    # Nd = demag_2(px[0,:],py[0,:],th[0], nmc)
    # print("Elapsed time: ",time.time()-t,"s")
    print(Nd)

# Calls Main Function
if __name__ == "__main__":
    main()
