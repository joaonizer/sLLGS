%% Inicio
clc
close all
%clear all
global alpha alpha_l Ms sig kbT q time_step gammamu0 mu0;
%% Constantes
q = 1.60217662e-19; % carga do eletron C
mu0=4*pi*1e-7;      % N/A2 ->   m.kg.s-2 A-2
kb=1.38064852e-23;  % J/K  ->   m2.kg.s-2.K-1
t2am=1/mu0;       % converte T para A/m
gammamu0=1.760859644e11*mu0;   %m/(sA)
hbar=2.05457e-34; % J.s/rad  -> h/2pi
%% Configuracoes do Algoritmo
if strcmp(computer,'GLNXA64')
    platform= 'lin';
else
    platform = 'win';
end

N = 60000;       % numero de passos
tempo_total=30e-9;% Tempo total de simulaÃ§Ã£o
alpha=0.05;%;0.054;
T=0;        % Kelvin
Ms=800e3;   % A/m
kbT=kb*T;   % J
ti = 0;     % instante inicial da variavel independente
[N,tempo_total,ti,tf,dt]=compute_Time(gammamu0,Ms,N,tempo_total,ti);
time_step=dt/(gammamu0*Ms); % segundos
if time_step>1e-12*alpha
    warning('Time_Step muito grande! Considere aumentar N!')
end
count_up=0;
%% Configuracoes do Sistema
name=['./Results/testNC/4_particles_she_down-' num2str(T) 'K-' num2str(N) 'steps-' num2str(tempo_total*1e9) 'ns-' num2str(alpha*100) 'alpha-force-module'];
grid=[
    0 1
    1 0
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
for i=2:part_n % inicializa as particulas de forma antiferromagnetica
    m(:,:,i)=(1)^(i-1)*m(:,:,1);
end
%for i=1:4
m(1,2,1)=1; % AND = 1

m(1,2,2)=1; %Porta A

% Dimensoes da Particula

w=ones(1,part_n)*50;  % width of particles

l=ones(1,part_n)*150; % length of particles

th=ones(1,part_n)*15;   %thickness of particles
%th(2)=15;

px=zeros(part_n,4);
py=px;
d_or=zeros(part_n,3);
nn=size(grid,1);
mm=size(grid,2);
count=1;
for i=1:mm
    for j=1:nn
        if grid(j,i)==1 % normal retangular
            cortes_y(count,:)=[0 0 0 0];%[25 25 -25 -25];%[35 0 0 -35]*1;
            count=count+1;
        elseif grid(j,i)==2 %and
            cortes_y(count,:)=[35 0 0 -35]*1;
            count=count+1;
        end
    end
end
% Transforms w,l and d into px py and d_or
for i=1:length(w)
    [px(i,:),py(i,:)]=write_Pontos(w,l,cortes_y(i,:),i);
end
dx=(w(1)+10)*[0:50];
dy=(l(1)+24)*[0:50]; %% deslocamentos em y
%dy=[l(1) 2*l(1)-10 3*l(1)-45 4*l(1)-80 5*l(1)-95]
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
%d_or(5,2)=d_or(5,2)+12.5;
%d_or(7,2)=d_or(7,2)+25;
%d_or(6,2)=d_or(6,2)-12.5;
%d_or(8,2)=d_or(8,2)-25;
%close all
%plot_Particles(px,py,d_or,dx,dy,cor,1,rows,cols,angles,name,eps);
%% Compute Tensores
compute_NCND=1; % If TRUE computes the tensors again
compute_PAR=0; % If TRUE uses paralel parfor to compute coupling tensor

if compute_NCND
    Nd=sparse(3*part_n,3*part_n);
    Nc=sparse(3*part_n,3*part_n);
    [Nd, V]=compute_Nd(px,py,th,part_n,compute_PAR,platform);
    %Nc=compute_Nc(px,py,th,d_or,radius,grid,part_n,compute_PAR,platform);
    Nc=compute_Nc2(px,py,th,d_or,part_n,compute_PAR,platform);
else
    warning('Tensores nao foram recalculados!');
end
%Nc=Nc2;
i=1;j=2;
ii=3*(i-1)+1:3*(i-1)+3;
ij=3*(j-1)+1:3*(j-1)+3;
fprintf('N: %3.4f\tC: %3.4f\t p: %3.4f\n',full(Nc2(ii(2),ij(2)))*1000,full(Nc(ii(2),ij(2)))*1000,full(Nc(ii(2),ij(2))/Nc2(ii(2),ij(2)))*100);
i=2;j=1;
ii=3*(i-1)+1:3*(i-1)+3;
ij=3*(j-1)+1:3*(j-1)+3;
fprintf('N: %3.4f\tC: %3.4f\t p: %3.4f\n',full(Nc2(ii(2),ij(2)))*1000,full(Nc(ii(2),ij(2)))*1000,full(Nc(ii(2),ij(2))/Nc2(ii(2),ij(2)))*100);
% Campo Aplicado
cor=zeros(part_n,3);
h_app=zeros(N+1,3,part_n);
a=1*150e-3; % T

colors = [
    0.4940    0.1840    0.5560  % Roxo
    0.3010    0.7450    0.9330  % Azul
    0.6350    0.0780    0.1840  % Vermelho
    0.4660    0.6740    0.1880  % Verde
    ];
for i=1:part_n
    phases=6;
    if sum(i==[1 2 6]) % X in
        cor(i,:)=colors(1,:);
        s=      [
            0   0   0   0   0   0   N/phases %1
            0   0   0   0   0   0   N/phases %2
            0   0   0   0   0   0   N/phases %3
            0   0   0   0   0   0   N/phases %4
            0   0   0   0   0   0   N/phases %5
            0   0   0   0   0   0   N/phases %6
            ];
    elseif (sum(i==[3 4 5]))
        cor(i,:)=colors(2,:);
        s=      [
            0   0   0   a   0   0   N/phases %1
            a   0   0   a   0   0   N/phases %2
            a   0   0   0   0   0   N/phases %3
            0   0   0   0   0   0   N/phases %4
            0   0   0   0   0   0   N/phases %5
            0   0   0   0   0   0   N/phases %6
            ];
    elseif (sum(i==[7 8 9]))
        cor(i,:)=colors(3,:);
        s=      [
            0   0   0   0   0   0   N/phases %1
            0   0   0   a   0   0   N/phases %2
            a   0   0   a   0   0   N/phases %3
            a   0   0   0   0   0   N/phases %4
            0   0   0   0   0   0   N/phases %5
            0   0   0   0   0   0   N/phases %6
            ];
    else
        cor(i,:)=colors(4,:);
        s=  [
            0   0   0   0   0   0   N/phases %1
            0   0   0   0   0   0   N/phases %2
            0   0   0   a   0   0   N/phases %3
            a   0   0   a   0   0   N/phases %4
            a   0   0   0   0   0   N/phases %5
            0   0   0   0   0   0   N/phases %6
            ];
    end
    h_app(:,:,i)=compute_Happ(N,s); % aplicado
end
close all
rows=1;
cols=1;
angles=[90 90];
eps=0;
name=['teste'];
plot_Particles(px,py,d_or,dx,dy,cor,1,rows,cols,angles,name,eps);
