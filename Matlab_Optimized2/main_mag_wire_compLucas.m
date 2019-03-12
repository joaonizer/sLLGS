%% Inicio
clc
close all
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

N = 40000;       % numero de passos
tempo_total=20e-9;% Tempo total de simulaÃ§Ã£o
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
    1 1
    ];
part_n=sum(sum(grid>0)); % quantidade de particulas

d_min=10; % distancia mÃ­nima entre as particulas


alpha_l=1/(1+alpha^2);
mi = [0 1 0]/sqrt(1); % valor inicial das variaveis dependente
m=zeros(N+1,3,part_n);
h_eff=zeros(N,3,part_n);
m(1,1,1)=0.1411;
m(1,2,1)=-0.99;
m(1,3,1)=0;
for i=2:part_n % inicializa as partÃ­culas de forma antiferromagnetica
    m(:,:,i)=(1)^(i-1)*m(:,:,1);
end
%% Dimensoes da Particula

w=ones(1,part_n)*50;  % width of particles

l=ones(1,part_n)*100; % length of particles

th=ones(1,part_n)*15;   %thickness of particles


px=zeros(part_n,4);
py=px;
d_or=zeros(part_n,3);
nn=size(grid,1);
mm=size(grid,2);
count=1;
for i=1:mm
    for j=1:nn
        if grid(j,i)==1 % normal retangular
            cortes_y(count,:)=[0 20 -20 0]*0;
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
dx=(w(1)+1000)*[0:50]; %% deslocamentos em x
dy=(l(1)+24)*[0:50]; %% deslocamentos em y
count=1;
for i=1:mm
    for j=1:nn
        if grid(j,i)
            d_or(count,:)=[dx(i) dy(nn-j+1) 0];
            count=count+1;
        end
    end
end

for teste=0:35
    %% Read Lucas Data
    name=['test' num2str(teste)];
    filename1=['C:\Users\joaon\Desktop\Teste3\' name '.csv'];
    filename2=['C:\Users\joaon\Desktop\Teste3\testdesc.csv'];
    test0 = importfileData(filename1);
    test0desc = importfileDesc(filename2);
    %%
    % M0
    px(1,:)=test0desc(1,1:4);
    py(1,:)=test0desc(1,5:8);
    th(1)=test0desc(1,9);
    d_or(1,:)=[test0desc(3,1:2) 0];
    m(1,:,1)=test0desc(3,3:5);
    % M1
    px(2,:)=test0desc(2,1:4);
    py(2,:)=test0desc(2,5:8);
    py(2,4)=py(2,4)+teste;
    th(2)=test0desc(2,9);
    d_or(2,:)=[test0desc(4,1:2) 0];
    m(1,:,2)=test0desc(4,3:5);
    %% Compute Tensores
    compute_NCND=1; % If TRUE computes the tensors again
    compute_PAR=0; % If TRUE uses paralel parfor to compute coupling tensor
    Nd=sparse(3*part_n,3*part_n);
    Nc=sparse(3*part_n,3*part_n);
    if compute_NCND
        [Nd, V]=compute_Nd2(px,py,th,part_n,compute_PAR,platform);
        Nc=compute_Nc2(px,py,th,d_or,part_n,compute_PAR,platform);
    else
        warning('Tensores nao foram recalculados!');
    end
    
    %% Modelagem do campo de Clock
    cor=zeros(part_n,3);
    h_app=zeros(N+1,3,part_n);
    a=1*150e-3; % Intesidade arbitraria
    
    colors = [
        0.4940    0.1840    0.5560  % Roxo
        0.3010    0.7450    0.9330  % Azul
        0.6350    0.0780    0.1840  % Vermelho
        0.4660    0.6740    0.1880  % Verde
        ];
    for i=1:part_n
        phases=4;
        if (sum(i==[1])) % X in
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
                0   0   0   0   0   0   N/phases %5
                0   0   0   a   0   0   N/phases %1
                a   0   0   a   0   0   N/phases %2
                a   0   0   0   0   0   N/phases %3
                0   0   0   0   0   0   N/phases %4
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
    bulk_sha = 0.4; % angulo de spin
    l_shm = 3.5; % [nm] comprimento de difusao de spin
    th_shm = test0desc(5,1); % [nm] espessura do material pesado
    J_shm=test0desc(5,2)*1e12; % densidade de corrente de spin (A/m2)
    theta_she=bulk_sha*(1-sech(th_shm/l_shm)); % Spin Hall Angle ()
    zeta=-hbar*theta_she*J_shm/2/q./th/1e-9/Ms;
    i_s=ones(N+1,3,part_n);
    i_s=h_app/a;
    h_app=zeros(N+1,3,part_n);
    for i=1:part_n
        i_s(:,:,i)=squeeze(i_s(:,:,i)).*zeta(i);
    end
    %% Campo Termico
    sig=sqrt(2*alpha*kb*T/mu0/Ms/Ms/V(1)); % desvio padrao do ruido termico
    dW=zeros(N+1,3,part_n);
    v = zeros(3,part_n);
    for j=1:part_n
        rng(j+1);
        dW(2:end,:,j)=(randn(N,3))*sqrt(dt);
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
    fprintf('Particles: \t\t%d\n',teste);
    fprintf('------------------------------------\n');
    dispstat('','init'); % One time only initialization
    dispstat(sprintf('Begining the process...'),'keepthis','timestamp');
    tic;
    hc=zeros(N+1,3,part_n); % inicializa o campo tÃ©rmico
    hd=hc;
    for i = 1 :N
        temp1=reshape(m(i,:,:),1,3*part_n);
        hd(i,:,:) = -reshape(temp1*Nd,3,part_n);
        hc(i,:,:) = -reshape(temp1*Nc,3,part_n);
        h_eff(i,:,:) = ...
            +squeeze(h_app(i,:,:)) ... % Campo externo aplicado (adimensional)
            +squeeze(hd(i,:,:)) ...    % Campo de desmagnetizacao
            +squeeze(hc(i,:,:));       % Campo de Acoplamento (adimensional)
        % RK_SDE
        m(i+1,:,:)=rk_sde(squeeze(m(i,:,:)),squeeze(h_eff(i,:,:)),squeeze(i_s(i,:,:)), v, dt,squeeze(dW(i,:,:)));
        module=sqrt(sum(m(i+1,:,:).^2));
        %module(module<1)=1;
        m(i+1,:,:)=m(i+1,:,:)./module; % Reprojecao
    end
    toc;
    dispstat('Finished.','keepprev');
    %% Plot Some Results
    close all
    t=0:1:N; % vetor tempo t
    t=t*time_step/1e-9; %transforma em ns
    angles=m(end,2,:)*90; % orientacao final de cada part�cula
    cols=mm; %numero de colunas no plot
    rows=nn; % numero de linhas
    eps=0; % 0 - nao plota resultado em .eps 1 - plota
    for i=1:part_n
        h_app(:,:,i)=squeeze(i_s(:,:,i))./zeta(i); % apenas converte o sinal em i_S para ser plotado na funcao seguinte
    end
    %plot_M_and_H(m,h_app,t,part_n,1,1,cols,rows,cor,grid,name,eps);
    %plot_Particles(px,py,d_or,dx,dy,cor,1,rows,cols,angles,name,eps);
    plot_ComparacaoLucas
    close all
end