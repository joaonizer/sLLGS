clc
clear all
bulk_sha = 0.4; % bul spin hall angle
th_shm = 1:0.1:20; % [nm] thickness of the spin hall material SHM
l_shm = 3.5%3.5; % [nm] SHM spin diffusion length
th=5;
Ap=50*th;
Am=180*th_shm;
theta_she=bulk_sha*(1-sech(th_shm/l_shm)); % Spin Hall Angle ()
J_shm=5*1.8e12; % Spin Hall current density (A/m2)

plot(th_shm,theta_she)
xlabel('Espessura do SHM (nm)')
ylabel('\theta_{SHE}')
sdf('P1')
%%
q = 1.60217662e-19; % carga do eletron C
mu0=4*pi*1e-7;      % N/A2 ->   m.kg.s-2
kb=1.38064852e-23;  % J/K  ->   m2.kg.s-2.K-1
A=13e-12;           % Exchange Stiffness J/m
t2am=1/mu0;       % converte T para A/m
gammamu0=1.760859644e11*mu0;   %m/(sA)
hbar=2.05457e-34; % J.s/rad  -> h/2pi
Ms=8e5; %A/m
th=1:0.2:25;
zeta=hbar*max(theta_she)*J_shm/2/q./th/1e-9/Ms;
plot(th,zeta)
xlabel('Espessura da Partícula (nm)')
ylabel('\zeta_{SHE}')
sdf('P1')

%%

clc
clear all
q = 1.60217662e-19; % carga do eletron C
mu0=4*pi*1e-7;      % N/A2 ->   m.kg.s-2
kb=1.38064852e-23;  % J/K  ->   m2.kg.s-2.K-1
A=13e-12;           % Exchange Stiffness J/m
t2am=1/mu0;       % converte T para A/m
gammamu0=1.760859644e11*mu0;   %m/(sA)
hbar=2.05457e-34; % J.s/rad  -> h/2pi
Ms=8e5; %A/m
bulk_sha = 0.4; % bul spin hall angle
th_shm = 1:1:25; % [nm] thickness of the spin hall material SHM
l_shm = 3.5%3.5; % [nm] SHM spin diffusion length
J_shm=1e12; % Spin Hall current density (A/m2)
w=40:1:60;
th=5;
for i=1:length(w)
    Ap=th*w(i);
    for j=1:length(th_shm)
        Am=(3*w(i)+30)*th_shm(j);
        theta_she=bulk_sha*(1-sech(th_shm(j)/l_shm)); % Spin Hall Angle ()
        zeta(i,j)=hbar*theta_she*J_shm/2/q./th/1e-9/Ms;
    end
end
figure;
surf(w,th_shm,zeta');
xlabel('Largura da Partícula (nm)')
ylabel('Espessura do SHM (nm)')
zlabel('\zeta_{SHE}')
%sdf('P1')
figure;
surf(w,th_shm,zeta');
xlabel('Largura da Partícula (nm)')
ylabel('Espessura do SHM (nm)')
zlabel('\zeta_{SHE}')
%sdf('P1')
figure;
semilogy(w,zeta(:,:));
xlabel('Largura da Partícula (nm)')
ylabel('\zeta_{SHE}')
%sdf('P1')

%% spin polarization

