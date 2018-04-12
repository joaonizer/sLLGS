%% Inicio
clc
close all
%more off
%clear all
%clear all
global alpha alpha_l Ms K1 HkMs sig kbT q time_step gammamu0 A n mu0;
%% Constantes
%Ku=0.26*1; %eV
q = -1.60217662e-19; % carga do eletron C
mu0=4*pi*1e-7;      % N/A2 ->   m.kg.s-2
kb=1.38064852e-23;  % J/K  ->   m2.kg.s-2.K-1
A=13e-12;           % Exchange Stiffness J/m
t2am=1/mu0;       % converte T para A/m
gammamu0=1.760859644e11*mu0;   %m/(sA)
hbar=2.05457e-34; % J.s/rad  -> h/2pi
%% Configuracoes do Algoritmo
if strcmp(computer,'GLNXA64')
    platform= 'lin';
else
    platform = 'win';
end
N = 40000;       % numero de passos
tempo_total=40e-9;% Tempo total de simulaÃ§Ã£o
alpha=0.05;%;0.054;
n = [0 1 0]/sqrt(1);
T=0;        % Kelvin
Ms=800e3;   % A/m
kbT=kb*T;   % J
ti = 0;     % instante inicial da variavel independente
[N,tempo_total,ti,tf,dt]=compute_Time(gammamu0,Ms,N,tempo_total,ti);
time_step=dt/(gammamu0*Ms); % segundos
if time_step>1e-12*alpha
    warning('Time_Step muito grande! Considere aumentar N!')
end
result_rk_up=zeros(N+1,3);
result_rk_down=result_rk_up;
count_up=0;
%% Configuracoes do Sistema
name=['./Results/2Particles/particles_she-' num2str(T) 'K-' num2str(N) 'steps-' num2str(tempo_total*1e9) 'ns-' num2str(alpha*100) 'alpha-force-module'];
grid=[
    1 0
    1 1
    ];
part_n=sum(sum(grid>0)); % quantidade de particulas

d_min=15; % distancia mÃ­nima entre as particulas


alpha_l=1/(1+alpha^2);
mi = [0 1 0]/sqrt(1); % valor inicial das variaveis dependente
m=zeros(N+1,3,part_n);
h_eff=zeros(N,3,part_n);
m(1,1,1)=0;
m(1,2,1)=-1;
m(1,3,1)=0;
for i=2:part_n % inicializa as partÃ­culas de forma antiferromagnetica
    m(:,:,i)=(-1)^(i-1)*m(:,:,1);
end
%% Dimensoes da Particula

w=ones(1,part_n)*50;  % width of particles

l=ones(1,part_n)*100; % length of particles

th=ones(1,part_n)*5;   %thickness of particles


px=zeros(part_n,4);
py=px;
d_or=zeros(part_n,3);
nn=size(grid,1);
mm=size(grid,2);
count=1;
for i=1:mm
    for j=1:nn
        if grid(j,i)==1 % normal retangular
            cortes_y(count,:)=[0 0 0 0];
            count=count+1;
        elseif grid(j,i)==2 %and
            cortes_y(count,:)=[0 0 0 -10];
            count=count+1;
        elseif grid(j,i)==3 %or
            cortes_y(count,:)=[30 0 0 0];
            count=count+1;
        end
    end
end
% Transforms w,l and d into px py and d_or
for i=1:length(w)
    [px(i,:),py(i,:)]=write_Pontos(w,l,cortes_y(i,:),i);
end

Ncs=zeros(30,9,9);

for distance=1:100
    dx=(w(1)+distance)*[0:50];
    dy=(l(1)+distance)*[0:50]; %% deslocamentos em y
    offset=0;
    count=1;
    for i=1:mm
        for j=1:nn
            if grid(j,i)
                d_or(count,:)=[dx(i) dy(nn-j+1) 0];
                count=count+1;
            end
        end
    end
    
    %d_or([3 4],1)=d_or([3 4],1)+25;
    
    %% Compute Tensores
    compute_NCND=1; % If TRUE computes the tensors again
    compute_PAR=0; % If TRUE uses paralel parfor to compute coupling tensor
    
    if compute_NCND
        Nd=sparse(3*part_n,3*part_n);
        Nc=sparse(3*part_n,3*part_n);
        [Nd, V]=compute_Nd(px,py,th,part_n,compute_PAR,platform);
        %[Nc,~,V]=compute_Tensores(px,py,d_or,th,part_n,platform,compute_PAR,Nc,Nd);
        %Nc_old=Nc;
        %Nd_old=Nd;
        %test_List
        radius=2;
        Nc=compute_Nc(px,py,th,d_or,radius,grid,part_n,compute_PAR,platform);
    else
        warning('Tensores nao foram recalculados!');
    end
    
    %% Nc's
    Ncs(distance,:,:)=Nc;
    %% Antiferro
    anti_ferro(distance,:)=[Nc(4,1), Nc(5,2), Nc(6,3)];
    %% Ferro
    ferro(distance,:)=[Nc(7,4), Nc(8,5), Nc(9,6)];
end

%% Plot
figure;
nc_code=['xyz'];
for i=1:3
    subplot(1,3,i)
    plot(abs(ferro(:,i)))
    hold on;
    plot(abs(anti_ferro(:,i)))
    xlabel('Distância (nm)')
    ylabel(['Nc_{' nc_code(i) nc_code(i) '}'])
    sdf('P1')
end


