#!/usr/bin/env python3

# !
# !	PROGRAMA DE CALCULO DO TENSOR DESMAGNETIZACAO DE
# !	POLIGONOS IRREGULARES (SLANTED)
# !
# !	O arquivo de entrada (IN_demag3D.dat) deve conter:
# !	px(1),py(1): coordenadas do canto superior esquerdo (nm)
# !	px(2),py(2): coordenadas do canto superior direito  (nm)
# !	px(3),py(3): coordenadas do canto inferior direito  (nm)
# !	px(4),py(4): coordenadas do canto inferior esquerdo (nm)
# !	t	   : espessura (nm)
# !
# !	*************** RESTRICAO ****************
# !	-------> px(1) e px(4) DEVEM ser iguais
# !	-------> px(2) e px(3) DEVEM ser iguais
# !	******************************************
# !
# !	O arquivo de saida (OUT_demag3D.dat) contem:
# !	Linha 1: d1, d2, d3
# !	Linha 2: d4, d5, d6
# !	Linha 3: d7, d8, d9
# !
# !	OBS.: OS FATORES DE DESMAGNETIZACAO SAO OBTIDOS FAZENDO
# !
# !	      D11=d1/(4*Pi*V), D12=d2/(4*Pi*V), D13=d3/(4*Pi*V)
# !	      D21=d4/(4*Pi*V), D22=d5/(4*Pi*V), D23=d6/(4*Pi*V)
# !	      D31=d7/(4*Pi*V), D32=d8/(4*Pi*V), D33=d9/(4*Pi*V)
# !
# !	      ONDE V e DADO EM nm^3
# !
# !	A energia de desmagnetizacao (eV) e calculada como
# !
# !	Ud = K/2*(
# !	          (d1*COS(phi)^2 + d5*SIN(phi)^2 + (d2+d4)*SIN(phi)*COS(phi))*SIN(theta)^2 +
# !		  ((d3+d7)*COS(phi) + (d6+d8)*SIN(phi))*SIN(theta)*COS(theta) +
# !		  d9*COS(theta)^2
# !		 )
# !	ONDE
# !
# !	K = mu0*Ms*Ms*JtoeV*ten**(-27)/(four*Pi)
# !	Ms = 800 kA/m
# !	mu0 = 4*Pi*10^-7 H/m
# !	JtoeV = 6.242_lg*ten**(18)
# !
# !
# !
# !

import numpy as np
# import numba as nb

# # @nb.jit(nopython = True, nogil = True)
# def xyz(px,ad,bd,au,bu,th):
#     x_0 = np.empty((3),np.float64)
#     x_0[0] = np.random.uniform(px[0], px[1])
#     y_min_0 = ad*x_0[0]+bd
#     y_max_0 = au*x_0[0]+bu
#     Dy_0 = y_max_0 - y_min_0
#     x_0[1] = np.random.uniform(y_min_0, y_max_0)
#     x_0[2] = np.random.uniform(-0.5*th, 0.5*th)
#
#     return x_0, Dy_0

# @nb.jit(nopython = True, nogil = True)
def fu(x_0, xp_0, Dy_0):

    den = np.sum((x_0 - xp_0) ** 2, 0) ** (-3/2)

    return Dy_0*(xp_0 - x_0)*den

# @nb.jit(nopython = True, nogil = True)
def fafbfc(px, py, th, x_0, Dy_0, au, bu, ad, bd,nmc):
    xp_0 = np.empty((3,nmc),np.float64)
    int_0_d = np.empty((3,6),np.float64)
    # S1
    xp_0[0,:] = px[1]
    xp_0[1,:] = np.random.uniform(py[2], py[1],nmc)
    xp_0[2,:] = np.random.uniform(-0.5*th, 0.5*th,nmc)
    int_0 = fu(x_0, xp_0, Dy_0)
    int_0_d[:,0] = np.sum(int_0)/nmc
    # S2
    xp_0[0,:] = np.random.uniform(px[0], px[1],nmc)
    xp_0[1,:] = au*xp_0[0,:]+bu
    xp_0[2,:] = np.random.uniform(-0.5*th, 0.5*th,nmc)
    int_0 = fu(x_0, xp_0, Dy_0)
    int_0_d[:,1] = np.sum(int_0)/nmc
    # S3
    xp_0[0,:] = np.random.uniform(px[0], px[1],nmc)
    y_min_p = ad*xp_0[0,:] + bd
    y_max_p = au*xp_0[0,:] + bu
    rnd_0 = np.random.sample(nmc)
    xp_0[1,:] = y_min_p*(1.0 - rnd_0) + rnd_0*y_max_p
    xp_0[2,:] = 0.5*th
    int_0 = fu(x_0, xp_0, Dy_0)
    int_0_d[:,2] = np.sum(int_0*(y_max_p - y_min_p))/nmc
    # S4
    xp_0[0,:] = px[0]
    xp_0[1,:] = np.random.uniform(py[3], py[0],nmc)
    xp_0[2,:] = np.random.uniform(-0.5*th, 0.5*th,nmc)
    int_0 = fu(x_0, xp_0, Dy_0)
    int_0_d[:,3] = np.sum(int_0)/nmc
    # S5
    xp_0[0,:] = np.random.uniform(px[0], px[1],nmc)
    xp_0[1,:] = ad*xp_0[0,:]+bd
    xp_0[2,:] = np.random.uniform(-0.5*th, 0.5*th)
    int_0 = fu(x_0, xp_0, Dy_0)
    int_0_d[:,4] = np.sum(int_0)/nmc
    # S6
    xp_0[0,:] = np.random.uniform(px[0], px[1])
    y_min_p = ad*xp_0[0,:]+bd
    y_max_p = au*xp_0[0,:]+bu
    rnd_0 = np.random.sample(nmc)
    xp_0[1,:] = y_min_p*(1.0 - rnd_0) + rnd_0*y_max_p
    xp_0[2,:] = -0.5*th
    int_0 = fu(x_0, xp_0, Dy_0)
    int_0_d[:,5] = np.sum(int_0*(y_max_p - y_min_p))/nmc

    return int_0_d

# @nb.jit(nopython = True, nogil = True)
def MC_demag(px,py,th,w,ad,bd,au,bu,nmc):
    s = np.empty(6)
    int_d = np.empty((3,6),np.float64)
    intd = np.empty((3,6),np.float64)
    x_0 = np.empty((3,nmc),np.float64)

    s[0] = w*th*th*(py[1] - py[2])
    s[1] = w*w*th*th
    s[2] = w*w*th
    s[3] = w*th*th*(py[0] - py[3])
    s[4] = w*w*th*th
    s[5] = w*w*th

    fd = np.empty((3,6),np.float64)

    x_0[0,:] = np.random.uniform(px[0], px[1], nmc)
    y_min_0 = ad*x_0[0,:]+bd
    y_max_0 = au*x_0[0,:]+bu
    Dy_0 = y_max_0 - y_min_0
    x_0[1,:] = np.random.sample(nmc)
    x_0[1,:] = y_min_0*(1.0 - x_0[1,:]) + x_0[1,:]*y_max_0
    x_0[2,:] = np.random.uniform(-0.5*th, 0.5*th, nmc)
    intd = fafbfc(px, py, th, x_0, Dy_0, au, bu, ad, bd,nmc)

    int_d = s*intd

    return int_d

# @nb.jit(nopython = True, nogil = True)
def demag3D3(px, py, th, nmc):

    d = np.empty(9)
    w = px[1] - px[0]
    au = (py[1]-py[0])/w
    bu = py[0] - au*px[0]
    ad = (py[2]-py[3])/w
    bd = py[2] - ad*px[2]

    dk = MC_demag(px,py,th,w,ad,bd,au,bu,nmc)

    teta2 = np.arctan(au)
    alfa2 = -np.sin(teta2)
    beta2 = np.cos(teta2)
    teta5 = np.arctan(ad)
    alfa5 = np.sin(teta5)
    beta5 = -np.cos(teta5)

    d[0] = dk[0,0] + alfa2*dk[0,1] - dk[0,3] + alfa5*dk[0,4]
    d[1] = beta2*dk[0,1] + beta5*dk[0,4]
    d[2] = dk[0,2] - dk[0,5]

    d[3] = dk[1,0] + alfa2*dk[1,1] - dk[1,3] + alfa5*dk[1,4]
    d[4] = beta2*dk[1,1] + beta5*dk[1,4]
    d[5] = dk[1,2] - dk[1,5]

    d[6] = dk[2,0] + alfa2*dk[2,1] - dk[2,3] + alfa5*dk[2,4]
    d[7] = beta2*dk[2,1] + beta5*dk[2,4]
    d[8] = dk[2,2] - dk[2,5]

    return d
