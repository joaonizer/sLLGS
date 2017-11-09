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
import random

def frand(kmin, kmax):
    rnd_0 = random.uniform(0, 1)
    f_rand = kmin*(1.0 - rnd_0) + rnd_0*kmax

    return f_rand


def xyz(px,ad,bd,au,bu,th):
    x_0 = np.zeros((3),np.float64)
    x_0[0] = frand(px[0], px[1])
    y_min_0 = ad*x_0[0]+bd
    y_max_0 = au*x_0[0]+bu
    Dy_0 = y_max_0 - y_min_0
    x_0[1] = frand(y_min_0, y_max_0)
    x_0[2] = frand(-0.5*th, 0.5*th)

    return x_0, Dy_0

def fu(x_0, xp_0, Dy_0):
    int_0 = np.zeros((3),np.float64)

    den = (
    np.sqrt(
    (x_0[0] - xp_0[0])*(x_0[0] - xp_0[0]) +
    (x_0[1] - xp_0[1])*(x_0[1] - xp_0[1]) +
    (x_0[2] - xp_0[2])*(x_0[2] - xp_0[2])
    )
    )**(-3)

    int_0[0] = Dy_0*(xp_0[0] - x_0[0])*den
    int_0[1] = Dy_0*(xp_0[1] - x_0[1])*den
    int_0[2] = Dy_0*(xp_0[2] - x_0[2])*den

    return int_0

def fafbfc(px, py, th, x_0, Dy_0, au, bu, ad, bd):
    xp_0 = np.zeros((3),np.float64)
    int_0 = np.zeros((3),np.float64)
    int_0_a = np.zeros((6),np.float64)
    int_0_b = np.zeros((6),np.float64)
    int_0_c = np.zeros((6),np.float64)
    # S1
    xp_0[0] = px[1]
    xp_0[1] = frand(py[2], py[1])
    xp_0[2] = frand(-0.5*th, 0.5*th)
    int_0 = fu(x_0, xp_0, Dy_0)
    int_0_a[0] = int_0[0]
    int_0_b[0] = int_0[1]
    int_0_c[0] = int_0[2]
    # S2
    xp_0[0] = frand(px[0], px[1])
    xp_0[1] = au*xp_0[0]+bu
    xp_0[2] = frand(-0.5*th, 0.5*th)
    int_0 = fu(x_0, xp_0, Dy_0)
    int_0_a[1] = int_0[0]
    int_0_b[1] = int_0[1]
    int_0_c[1] = int_0[2]
    # S3
    xp_0[0] = frand(px[0], px[1])
    y_min_p = ad*xp_0[0] + bd
    y_max_p = au*xp_0[0] + bu
    xp_0[1] = frand(y_min_p, y_max_p)
    xp_0[2] = 0.5*th
    int_0 = fu(x_0, xp_0, Dy_0)
    int_0_a[2] = int_0[0]*(y_max_p - y_min_p)
    int_0_b[2] = int_0[1]*(y_max_p - y_min_p)
    int_0_c[2] = int_0[2]*(y_max_p - y_min_p)
    # S4
    xp_0[0] = px[0]
    xp_0[1] = frand(py[3], py[0])
    xp_0[2] = frand(-0.5*th, 0.5*th)
    int_0 = fu(x_0, xp_0, Dy_0)
    int_0_a[3] = int_0[0]
    int_0_b[3] = int_0[1]
    int_0_c[3] = int_0[2]
    # S5
    xp_0[0] = frand(px[0], px[1])
    xp_0[1] = ad*xp_0[0]+bd
    xp_0[2] = frand(-0.5*th, 0.5*th)
    int_0 = fu(x_0, xp_0, Dy_0)
    int_0_a[4] = int_0[0]
    int_0_b[4] = int_0[1]
    int_0_c[4] = int_0[2]
    # S6
    xp_0[0] = frand(px[0], px[1])
    y_min_p = ad*xp_0[0]+bd
    y_max_p = au*xp_0[0]+bu
    xp_0[1] = frand(y_min_p, y_max_p)
    xp_0[2] = -0.5*th
    int_0 = fu(x_0, xp_0, Dy_0)
    int_0_a[5] = int_0[0]*(y_max_p - y_min_p)
    int_0_b[5] = int_0[1]*(y_max_p - y_min_p)
    int_0_c[5] = int_0[2]*(y_max_p - y_min_p)

    return int_0_a, int_0_b, int_0_c

def MC_demag(px,py,th,w,ad,bd,au,bu,nmc):
    s = np.zeros(6)
    int_a = np.zeros((6),np.float64)
    int_b = np.zeros((6),np.float64)
    int_c = np.zeros((6),np.float64)
    inta = np.zeros((6),np.float64)
    intb = np.zeros((6),np.float64)
    intc = np.zeros((6),np.float64)

    s[0] = w*th*th*(py[1] - py[2])
    s[1] = w*w*th*th
    s[2] = w*w*th
    s[3] = w*th*th*(py[0] - py[3])
    s[4] = w*w*th*th
    s[5] = w*w*th

    fa = np.zeros((6),np.float64)
    fb = np.zeros((6),np.float64)
    fc = np.zeros((6),np.float64)


    for i in range(int(nmc)):
        x_0, Dy_0 = xyz(px,ad,bd,au,bu,th)
        inta, intb, intc = fafbfc(px, py, th, x_0, Dy_0, au, bu, ad, bd)
        fa = fa + inta
        fb = fb + intb
        fc = fc + intc

    #print("fa:\n",fa,"\nfb:\n",fb,"\nfc:\n",fc,"\n")
    int_a = s*fa/nmc
    int_b = s*fb/nmc
    int_c = s*fc/nmc
    #print("int_a:\n",int_a,"\nint_b:\n",int_b,"\nint_c:\n",int_c,"\n")
    return int_a, int_b, int_c

def demag3D3(px, py, th, nmc):

    #for i in range(10):
    #    print("RND",i,": ",random.uniform(0, 1),"\n")

    d = np.zeros(9)
    w = px[1] - px[0]
    #print("w: ",w)
    au = (py[1]-py[0])/w
    #print("au: ", au)
    bu = py[0] - au*px[0]
    #print("bu: ", bu)
    ad = (py[2]-py[3])/w
    #print("ad: ", ad)
    bd = py[2] - ad*px[2]
    #print("bd: ", bd)

    ak = np.zeros((6),np.float64)
    bk = np.zeros((6),np.float64)
    ck = np.zeros((6),np.float64)

    ak, bk, ck = MC_demag(px,py,th,w,ad,bd,au,bu,nmc)

    #print("int_a:\n",ak,"\nint_b:\n",bk,"\nint_c:\n",ck,"\n")

    teta2 = np.arctan(au)
    #print("teta2: ", teta2)
    alfa2 = -np.sin(teta2)
    #print("alfa2: ", alfa2)
    beta2 = np.cos(teta2)
    #print("beta2: ", beta2)

    teta5 = np.arctan(ad)
    #print("teta5: ", teta5)
    alfa5 = np.sin(teta5)
    #print("alfa5: ", alfa5)
    beta5 = -np.cos(teta5)
    #print("beta5: ", beta5)

    d[0] = ak[0] + alfa2*ak[1] - ak[3] + alfa5*ak[4]
    d[1] = beta2*ak[1] + beta5*ak[4]
    d[2] = ak[2] - ak[5]

    d[3] = bk[0] + alfa2*bk[1] - bk[3] + alfa5*bk[4]
    d[4] = beta2*bk[1] + beta5*bk[4]
    d[5] = bk[2] - bk[5]

    d[6] = ck[0] + alfa2*ck[1] -ck[3] +alfa5*ck[4]
    d[7] = beta2*ck[1] + beta5*ck[4]
    d[8] = ck[2] - ck[5]

    #for i in range(9):
    #    print(d[i])

    return d
