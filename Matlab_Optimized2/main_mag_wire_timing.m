%% Inicio
clc
close all
%more off
clear all
%clear all
global alpha alpha_l Ms sig kbT q time_step gammamu0 mu0;
%% Constantes
%Ku=0.26*1; %eV
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
runs=30;
N = 4000;       % numero de passos
tempo_total=2e-9;% Tempo total de simulaÃ§Ã£o
alpha=0.05;%;0.054;
n = [0 1 0]/sqrt(1);
T=10;        % Kelvin
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
tempo=zeros(5,runs);
for grid_size=1:5
    grid=ones(2^(grid_size));
    part_n=sum(sum(grid>0)); % quantidade de particulas
    
    d_min=15; % distancia mÃ­nima entre as particulas
    
    
    alpha_l=1/(1+alpha^2);
    mi = [0 1 0]/sqrt(1); % valor inicial das variaveis dependente
    m=zeros(N+1,3,part_n);
    h_eff=zeros(N,3,part_n);
    m(1,1,1)=0.98;
    m(1,2,:)=0.199;
    m(1,3,1)=0;
    for i=2:part_n % inicializa as partÃ­culas de forma antiferromagnetica
        m(:,:,i)=(1)^(i-1)*m(:,:,1);
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
    dx=(w(1)+10)*[0:130];
    dy=(l(1)+24)*[0:130]; %% deslocamentos em y
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
        tic;
        [Nd, V]=compute_Nd(px,py,th,part_n,compute_PAR,platform);
        tempo_Nd(grid_size)=toc;
        %[Nc,~,V]=compute_Tensores(px,py,d_or,th,part_n,platform,compute_PAR,Nc,Nd);
        %radius=2;
        %Nc=compute_Nc(px,py,th,d_or,radius,grid,part_n,compute_PAR,platform);
        tic;
        Nc=compute_Nc(px,py,th,d_or,2,grid,part_n,compute_PAR,platform);
        tempo_Nc(grid_size)=toc;
    else
        warning('Tensores nao foram recalculados!');
    end
for run=1:runs
    %Nc=0.5*Nc;
    % Nc(:,:,15,[9,11])=zeros(3,3,1,2);
    % Nc(:,:,[9,11],15)=zeros(3,3,2,1);
    % Nc(:,:,16,[12,14])=zeros(3,3,1,2);
    % Nc(:,:,[12,14],16)=zeros(3,3,2,1);
    % Nc(:,:,24,[19,21])=zeros(3,3,1,2);
    % Nc(:,:,[19,21],24)=zeros(3,3,2,1);
    % Nc(:,:,16,22)=zeros(3,3,1,1);
    % Nc(:,:,22,16)=zeros(3,3,1,1);
    % Nc(:,:,15,18)=zeros(3,3,1,1);
    % Nc(:,:,18,15)=zeros(3,3,1,1);
    %% Campo Aplicado
    h_app=zeros(N+1,3,part_n);
    i_s=h_app;
    i_s(:,1,:)=ones(N+1,1,part_n)*0.15;
    sig=sqrt(2*alpha*kb*T/mu0/Ms/Ms/V(1)); % new
    dW=zeros(N+1,3,part_n);
    hT=dW;
    v = zeros(3,part_n);
    for j=1:part_n
        %rng(jj+1);
        dW(2:end,:,j)=(randn(N,3))*sqrt(dt);
        hT(:,:,j)=sig*dW(:,:,j)*sqrt(V(1)/V(j));
        v(:,j)=[sig sig sig]*sqrt(V(1)/V(j));
    end
    %% Metodo para solucao numerica

        clc
        fprintf('\n------------------------------------\n');
        fprintf('            Range-Kutta            \n');
        fprintf('------------------------------------\n');
        fprintf('Plataform:\t\t%s\n',platform);
        fprintf('dt real:\t\t%3.3e s\n',time_step);
        fprintf('dt line:\t\t%3.3e\n',dt);
        fprintf('total_time:\t\t%3.3e s\n',tempo_total);
        fprintf('N: \t\t\t%d\n',N);
        fprintf('T: \t\t\t%d Kelvin\n',T);
        fprintf('alpha: \t\t\t%.3f\n',alpha);
        fprintf('Particles: \t\t%d\n',part_n);
        fprintf('Experiment No: \t%d Grid %d\n',run,grid_size);
        fprintf('------------------------------------\n');
        dispstat('','init'); % One time only initialization
        dispstat(sprintf('Begining the process...'),'keepthis','timestamp');
        hc=zeros(N+1,3,part_n); % inicializa o campo tÃ©rmico
        hd=hc;
        n_part_n=ones(part_n,1)*n;
        tic;
        fid=fopen('results.txt','w');
        for i = 1 :N
            %progress=i/N*100;
            %dispstat(sprintf('Progress %3.2f %%',progress),'timestamp');
            
            temp1=reshape(m(i,:,:),1,3*part_n);
            hd(i,:,:) = -reshape(temp1*Nd,3,part_n);
            hc(i,:,:) = -reshape(temp1*Nc,3,part_n);
            %+squeeze(hT(i,:,:))...                  % Campo termico (adimensional)
            h_eff(i,:,:) = ...
                +squeeze(h_app(i,:,:)) ...           % Campo externo aplicado (adimensional)
                +squeeze(hd(i,:,:)) ...
                +squeeze(hc(i,:,:));                    % Campo de Acoplamento (adimensional)
            m(i+1,:,:)=rk_sde(squeeze(m(i,:,:)),squeeze(h_eff(i,:,:)),squeeze(i_s(i,:,:)), v, dt,squeeze(dW(i,:,:)));
            m(i+1,:,:)=m(i+1,:,:)./sqrt(sum(m(i+1,:,:).^2));
            if ~mod(i,10)
                fprintf(fid,'%.4f ',m(i+1,:,:));
            end
        end
        fclose all;
        tempo(grid_size,run)=toc;
        pause(1)
        delete results.txt
end
%save tempo tempo
%save tempo_Nc tempo_Nc
%save tempo_Nd tempo_Nd
end
