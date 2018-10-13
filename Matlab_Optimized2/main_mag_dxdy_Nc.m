%% Inicio
% Neste experimento ser�o variadas as dist�ncias entre duas part;iculas
% dx e dy
% para verificar qual valor de threshold � considerado ideal para zerar
% os termos do tensor Nc
% Nd ser� mantido constante para todos as rodadas,e vitando que o campo
% de desmagnetiza��o influencie na press�o da magnetiza��o entre uma
% rodada e outa.
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

N = 2500;       % numero de passos
tempo_total=2.5e-9;% Tempo total de simulaÃ§Ã£o
alpha=1;%;0.054;
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
name=['./Results/testNC/4_particles_she_down-' num2str(T) 'K-' num2str(N) 'steps-' num2str(tempo_total*1e9) 'ns-' num2str(alpha*100) 'alpha-force-module'];
grid=[
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
    m(1,:,i)=[0.9999 -0.0141 0];%(-1)^(i-1)*m(:,:,1);
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
            cortes_y(count,:)=[0 0 0 -25];
            count=count+1;
        elseif grid(j,i)==3 %or
            cortes_y(count,:)=[25 0 0 0];
            count=count+1;
        elseif grid(j,i)==4 %or
            cortes_y(count,:)=[0 25 -25 0];
            count=count+1;
        end
    end
end
% Transforms w,l and d into px py and d_or
for i=1:length(w)
    [px(i,:),py(i,:)]=write_Pontos(w,l,cortes_y(i,:),i);
end

    m1=zeros(N+1,3,2,50,2);
    Nc1=zeros(6,6,50,2);

for dxx=10:10
dx=(w(1)+dxx)*[0:3];
dy=(l(1)+24)*[0:3]; %% deslocamentos em y
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
for jj=1:2
%% Compute Tensores
compute_NCND=1; % If TRUE computes the tensors again
compute_PAR=0; % If TRUE uses paralel parfor to compute coupling tensor

if compute_NCND
    %Nd=sparse(3*part_n,3*part_n);
    %[Nd, V]=compute_Nd(px,py,th,part_n,compute_PAR,platform);
    %[Nc,~,V]=compute_Tensores(px,py,d_or,th,part_n,platform,compute_PAR,Nc,Nd);
    %radius=2;
    %Nc=compute_Nc(px,py,th,d_or,radius,grid,part_n,compute_PAR,platform);
    if jj==1
        Nc=zeros(3*part_n,3*part_n);
        Nc=compute_Nc2(px,py,th,d_or,part_n,compute_PAR,platform);
    else
        Nc=zeros(3*part_n,3*part_n);
    end
else
    warning('Tensores nao foram recalculados!');
end
%% Campo Aplicado
cor=zeros(part_n,3);
    h_app=zeros(N+1,3,part_n);
    a=1*150e-3; % T
    
    % EXP tem todas as combinaÃ§Ãµes de entradas
    %    X  Y
    experiment=0.1*[1  1
        1 -1
        -1  1
        -1 -1
        ];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ex=experiment(jj,:); % testa a primeira combinaÃ§Ã£o de 8 possibilidades
    %%%%%% ^ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    colors = [
        0.4940    0.1840    0.5560  % Roxo
        0.3010    0.7450    0.9330  % Azul
        0.6350    0.0780    0.1840  % Vermelho
        0.4660    0.6740    0.1880  % Verde
        ];
    for i=1:part_n
        phases=4;
        if (i==1) % X in
            cor(i,:)=colors(1,:);
            s=      [
                0   0   0   0   0   0   N/phases %2
                0   0   0   0   0   0   N/phases %2
                0   0   0   0   0   0   N/phases %2
                0   0   0   0   0   0   N/phases %2
                ];
        elseif (sum(i==[2]))
            cor(i,:)=colors(2,:);
            s=      [
                0   0   0   a   0   0   N/phases %1
                a   0   0   a   0   0   N/phases %2
                a   0   0   0   0   0   N/phases %3
                0   0   0   0   0   0   N/phases %4
                ];
        elseif (sum(i==[3]))
            cor(i,:)=colors(3,:);
            s=      [
                0   0   0   a   0   0   N/phases %1
                a   0   0   a   0   0   N/phases %2
                a   0   0   0   0   0   N/phases %3
                0   0   0   0   0   0   N/phases %4
                0   0   0   0   0   0   N/phases %5
                ];
        else
            cor(i,:)=colors(4,:);
            s=  [
                0   0   0   0   0   0   N/phases %1
                0   0   0   a   0   0   N/phases %2
                a   0   0   a   0   0   N/phases %3
                a   0   0   0   0   0   N/phases %4
                0   0   0   0   0   0   N/phases %5
                ];
        end
        h_app(:,:,i)=compute_Happ(N,s); % aplicado
    end
%% Corrente de Spin
% Define a curva da corrente de spin aplicada
bulk_sha = 0.4; % bul spin hall angle
th_shm = 5; % [nm] thickness of the spin hall material SHM
l_shm = 3.5; % [nm] SHM spin diffusion length
theta_she=bulk_sha*(1-sech(th_shm/l_shm)); % Spin Hall Angle ()
J_shm=0*1.8e12; % Spin Hall current density (A/m2)

zeta=hbar*theta_she*J_shm/2/q./th/1e-9/Ms;
%Ns = 2*Ms*V/gammamu0/hbar;
%is=I_s./(q*gammamu0*mu0*Ms*Ns); % magnitude normalizada da corrente de spin
i_s=ones(N+1,3,part_n);
i_s=h_app/a;
h_app=zeros(N+1,3,part_n);
for i=1:part_n
    i_s(:,:,i)=squeeze(i_s(:,:,i)).*zeta(i);
end
    %% Campo TÃ©rmico
    K1=0;
    HkMs=2*K1/Ms/mu0/Ms;
    % dt (adimensional) --> dt_real = dt/gammamu0*Ms
    %sig=sqrt(2*alpha*kbT/mu0/V(1)/dt)/Ms;
    
    %sig=sqrt(2*alpha*kb*T/gammamu0/mu0/Ms/V(1)/time_step)/Ms; % old version
    %sig=sqrt(2*alpha*kb*T/mu0/Ms/Ms/V(1))*sqrt(dt); % new
    sig=sqrt(2*alpha*kb*T/gammamu0/Ms/mu0/V(1)/time_step)/Ms*sqrt(dt);
    %version
    dW=zeros(N+1,3,part_n);
    hT=dW;
    v = zeros(3,part_n);
    for j=1:part_n
        rng(jj+1);
        dW(2:end,:,j)=(randn(N,3));
        hT(:,:,j)=sig*dW(:,:,j)*sqrt(V(1)/V(j));
        v(:,j)=[sig sig sig]*sqrt(V(1)/V(j));
    end
    %% Metodo para solucao numerica
    clc
%     fprintf('\n------------------------------------\n');
%     fprintf('            Range-Kutta            \n');
%     fprintf('------------------------------------\n');
%     fprintf('Plataform:\t\t%s\n',platform);
%     fprintf('dt real:\t\t%3.3e s\n',time_step);
%     fprintf('dt line:\t\t%3.3e\n',dt);
%     fprintf('total_time:\t\t%3.3e s\n',tempo_total);
%     fprintf('N: \t\t\t%d\n',N);
%     fprintf('T: \t\t\t%d Kelvin\n',T);
%     fprintf('alpha: \t\t\t%.3f\n',alpha);
%     fprintf('Particles: \t\t%d\n',part_n);
%     fprintf('Experiment No: \t%d\n',rodada);
    %fprintf('------------------------------------\n');
    %dispstat('','init'); % One time only initialization
    %dispstat(sprintf('Begining the process...'),'keepthis','timestamp');
    tic;
    hc=zeros(N+1,3,part_n); % inicializa o campo tÃ©rmico
    hd=hc;
    n_part_n=ones(part_n,1)*n;
    fprintf('Experiment No: %d/%d\n',dxx,jj);
    for i = 1 :N
        %progress=i/N*100;
        %dispstat(sprintf('Progress %3.2f %%',progress),'timestamp');
        
        temp1=reshape(m(i,:,:),1,3*part_n);
        hd(i,:,:) = -reshape(temp1*Nd,3,part_n);
        hc(i,:,:) = -reshape(temp1*Nc,3,part_n);
        %+squeeze(hT(i,:,:))...                  % Campo termico (adimensional)
        h_eff(i,:,:) = ...
            +squeeze(h_app(i,:,:)) ...           % Campo externo aplicado (adimensional)
            +transpose(HkMs*dot(n_part_n,squeeze(m(i,:,:))',2)*n) ...	% Anistropia magnetocristalina (adimensional)
            +squeeze(hd(i,:,:)) ...   
        +squeeze(hc(i,:,:));                    % Campo de Acoplamento (adimensional)
        % Range-Kuta-4
        % m(i+1,:,:)=rk4(squeeze(m(i,:,:)),squeeze(h_eff(i,:,:)),squeeze(hT(i,:,:)),squeeze(i_s(i,:,:)),dt);
        % RK_SDE
        m(i+1,:,:)=rk_sde(squeeze(m(i,:,:)),squeeze(h_eff(i,:,:)),squeeze(i_s(i,:,:)), v, dt,squeeze(dW(i,:,:)));
        m(i+1,:,:)=m(i+1,:,:)./sqrt(sum(m(i+1,:,:).^2));
    end
    m1(:,:,:,(dxx+1)/2,jj)=m;
    Nc1(:,:,(dxx+1)/2,jj)=Nc;

end
 
end