function [d] = demag3D3(px, py, th, nmc)
% !	PROGRAMA DE CALCULO DO TENSOR DESMAGNETIZACAO DE
% !	POLIGONOS IRREGULARES (SLANTED)
% !
% !	O arquivo de entrada (IN_demag3D.dat) deve conter:
% !	px(1),py(1) coordenadas do canto superior esquerdo (nm)
% !	px(2),py(2) coordenadas do canto superior direito  (nm)
% !	px(3),py(3) coordenadas do canto inferior direito  (nm)
% !	px(4),py(4) coordenadas do canto inferior esquerdo (nm)
% !	t	   : espessura (nm)
% !
% !	*************** RESTRICAO ****************
% !	-------> px(1) e px(4) DEVEM ser iguais
% !	-------> px(2) e px(3) DEVEM ser iguais
% !	******************************************
% !
% !	O arquivo de saida (OUT_demag3D.dat) contem:
% !	Linha 1: d1, d2, d3
% !	Linha 2: d4, d5, d6
% !	Linha 3: d7, d8, d9
% !
% !	OBS.: OS FATORES DE DESMAGNETIZACAO SAO OBTIDOS FAZENDO
% !
% !	      D11=d1/(4*Pi*V), D12=d2/(4*Pi*V), D13=d3/(4*Pi*V)
% !	      D21=d4/(4*Pi*V), D22=d5/(4*Pi*V), D23=d6/(4*Pi*V)
% !	      D31=d7/(4*Pi*V), D32=d8/(4*Pi*V), D33=d9/(4*Pi*V)
% !
% !	      ONDE V e DADO EM nm^3
% !
% !	A energia de desmagnetizacao (eV) e calculada como
% !
% !	Ud = K/2*(
% !	          (d1*COS(phi)^2 + d5*SIN(phi)^2 + (d2+d4)*SIN(phi)*COS(phi))*SIN(theta)^2 +
% !		  ((d3+d7)*COS(phi) + (d6+d8)*SIN(phi))*SIN(theta)*COS(theta) +
% !		  d9*COS(theta)^2
% !		 )
% !	ONDE
% !
% !	K = mu0*Ms*Ms*JtoeV*ten**(-27)/(four*Pi)
% !	Ms = 800 kA/m
% !	mu0 = 4*Pi*10^-7 H/m
% !	JtoeV = 6.242_lg*ten**(18)
% !

d = zeros(9,1);
w = px(2) - px(1);
au = (py(2)-py(1))/w;
bu = py(1) - au*px(1);
ad = (py(3)-py(4))/w;
bd = py(3) - ad*px(3);

[ak, bk, ck] = MC_demag(px,py,th,w,ad,bd,au,bu,nmc);

teta2 = atan(au);
alfa2 = -sin(teta2);
beta2 = cos(teta2);
teta5 = atan(ad);
alfa5 = sin(teta5);
beta5 = -cos(teta5);

d(1) = ak(1) + alfa2*ak(2) - ak(4) + alfa5*ak(5);
d(2) = beta2*ak(2) + beta5*ak(5);
d(3) = ak(3) - ak(6);

d(4) = bk(1) + alfa2*bk(2) - bk(4) + alfa5*bk(5);
d(5) = beta2*bk(2) + beta5*bk(5);
d(6) = bk(3) - bk(6);

d(7) = ck(1) + alfa2*ck(2) -ck(4) +alfa5*ck(5);
d(8) = beta2*ck(2) + beta5*ck(5);
d(9) = ck(3) - ck(6);

end

function f_rand = frand(kmin, kmax)
rnd_0 = rand();
f_rand = kmin*(1.0 - rnd_0) + rnd_0*kmax;
end


function [x_0, Dy_0] = xyz(px,ad,bd,au,bu,th)
x_0 = zeros(3,1);
x_0(1) = frand(px(1), px(2));
y_min_0 = ad*x_0(1)+bd;
y_max_0 = au*x_0(1)+bu;
Dy_0 = y_max_0 - y_min_0;
x_0(2) = frand(y_min_0, y_max_0);
x_0(3) = frand(-0.5*th, 0.5*th);
end

function int_0 = fu(x_0, xp_0, Dy_0)
int_0 = zeros(3,1);
den = (...
    sqrt(...
    (x_0(1) - xp_0(1))*(x_0(1) - xp_0(1)) +...
    (x_0(2) - xp_0(2))*(x_0(2) - xp_0(2)) +...
    (x_0(3) - xp_0(3))*(x_0(3) - xp_0(3))...
    )...
    )^(-3);
int_0(1) = Dy_0*(xp_0(1) - x_0(1))*den;
int_0(2) = Dy_0*(xp_0(2) - x_0(2))*den;
int_0(3) = Dy_0*(xp_0(3) - x_0(3))*den;
end

function [int_0_a, int_0_b, int_0_c] = fafbfc(px, py, th, x_0, Dy_0, au, bu, ad, bd)
xp_0 = zeros(3,1);
int_0_a = zeros(6,1);
int_0_b = zeros(6,1);
int_0_c = zeros(6,1);
% S1
xp_0(1) = px(2);
xp_0(2) = frand(py(3), py(2));
xp_0(3) = frand(-0.5*th, 0.5*th);
int_0 = fu(x_0, xp_0, Dy_0);
int_0_a(1) = int_0(1);
int_0_b(1) = int_0(2);
int_0_c(1) = int_0(3);
% S2
xp_0(1) = frand(px(1), px(2));
xp_0(2) = au*xp_0(1)+bu;
xp_0(3) = frand(-0.5*th, 0.5*th);
int_0 = fu(x_0, xp_0, Dy_0);
int_0_a(2) = int_0(1);
int_0_b(2) = int_0(2);
int_0_c(2) = int_0(3);
% S3
xp_0(1) = frand(px(1), px(2));
y_min_p = ad*xp_0(1) + bd;
y_max_p = au*xp_0(1) + bu;
xp_0(2) = frand(y_min_p, y_max_p);
xp_0(3) = 0.5*th;
int_0 = fu(x_0, xp_0, Dy_0);
int_0_a(3) = int_0(1)*(y_max_p - y_min_p);
int_0_b(3) = int_0(2)*(y_max_p - y_min_p);
int_0_c(3) = int_0(3)*(y_max_p - y_min_p);
% S4
xp_0(1) = px(1);
xp_0(2) = frand(py(4), py(1));
xp_0(3) = frand(-0.5*th, 0.5*th);
int_0 = fu(x_0, xp_0, Dy_0);
int_0_a(4) = int_0(1);
int_0_b(4) = int_0(2);
int_0_c(4) = int_0(3);
% S5
xp_0(1) = frand(px(1), px(2));
xp_0(2) = ad*xp_0(1)+bd;
xp_0(3) = frand(-0.5*th, 0.5*th);
int_0 = fu(x_0, xp_0, Dy_0);
int_0_a(5) = int_0(1);
int_0_b(5) = int_0(2);
int_0_c(5) = int_0(3);
% S6
xp_0(1) = frand(px(1), px(2));
y_min_p = ad*xp_0(1)+bd;
y_max_p = au*xp_0(1)+bu;
xp_0(2) = frand(y_min_p, y_max_p);
xp_0(3) = -0.5*th;
int_0 = fu(x_0, xp_0, Dy_0);
int_0_a(6) = int_0(1)*(y_max_p - y_min_p);
int_0_b(6) = int_0(2)*(y_max_p - y_min_p);
int_0_c(6) = int_0(3)*(y_max_p - y_min_p);
end

function [int_a, int_b, int_c] = MC_demag(px,py,th,w,ad,bd,au,bu,nmc)
s = zeros(6,1);
s(1) = w*th*th*(py(2) - py(3));
s(2) = w*w*th*th;
s(3) = w*w*th;
s(4) = w*th*th*(py(1) - py(4));
s(5) = w*w*th*th;
s(6) = w*w*th;
fa = zeros(6,1);
fb = zeros(6,1);
fc = zeros(6,1);



for i =1:nmc
    [x_0, Dy_0] = xyz(px,ad,bd,au,bu,th);
    [inta, intb, intc] = fafbfc(px, py, th, x_0, Dy_0, au, bu, ad, bd);
    fa = fa + inta;
    fb = fb + intb;
    fc = fc + intc;
end
int_a = s.*fa/nmc;
int_b = s.*fb/nmc;
int_c = s.*fc/nmc;
end
