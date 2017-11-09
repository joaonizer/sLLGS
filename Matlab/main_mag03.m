%% Inicio
clc
close all
clear all
%clear all
global alpha alpha_l Ms K1 HkMs sig kbT q time_step gammamu0 A T n;
%% Constantes
%Ku=0.26*1; %eV
q = 1.60217662e-19; % carga do eletron C
mu0=4*pi*1e-7;      % H/m ou T.m/A
kb=1.38064852e-23;  % m2.kg.s-2.K-1
A=13e-12;           % Exchange Stiffness J/m
t2am=7.943e5;       % converte T para A/m
gammamu0=2.211e5;   % (sA/m)-1
hbar=1.054571800e-34; % J.s/rad  -> h/2pi
%% Configuracoes do Algoritmo
N = 5000;       % numero de passos
tempo_total=100e-9;% Tempo total de simulação
alpha=1.0;
n = [0 1 0];
T=0;        % Kelvin
Ms=800e3;   % A/m
kbT=kb*T;   % J
ti = 0;     % instante inicial da variavel independente
[N,tempo_total,ti,tf,dt]=compute_Time(gammamu0,Ms,N,tempo_total,ti);
time_step=dt/(gammamu0*Ms); % segundos
if time_step>1e-12*alpha
    warning('Time_Step muito grande! Considere aumentar N!')
end
%% Configuracoes do Sistema


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part_n=23; % quantidade de particulas
%%%%%% ^ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ^ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_min=10; % distancia mínima entre as partículas


alpha_l=1/(1+alpha^2);
mi = [1 0 0]/sqrt(1); % valor inicial das variaveis dependente
m=zeros(N+1,3,part_n);
h_eff=zeros(N,3);
m(1,1,1)=0;
m(1,2,1)=1;
m(1,3,1)=0;
for i=2:part_n % inicializa as partículas de forma antiferromagnetica
    m(:,:,i)=(-1)^(i-1)*m(:,:,1);
end
%% Dimensoes da Particula
platform='lin';

w=ones(1,part_n)*50;  % width of particles

l=ones(1,part_n)*100; % length of particles

th=ones(1,part_n)*10;   %thickness of particles

d=[0,ones(1,part_n-1)]*d_min; % distance among particles

px=zeros(part_n,4);
py=px;
d_or=zeros(part_n,3);
cortes_y=zeros(part_n,4);

cortes_y=[
    0   0   0   0
    0   0   0   0
    0   0   0   0
    0   0   0   0
    0   0   0   0
    0   0   0   0
    0   0   0   -20
    0   0   0   0
    0   0   0   0
    
    0   0   0   0
    0   0   0   0
    0   0   0   0
    0   0   0   0
    0   0   0   0
    0   0   0   0
    0   0   0   -20
    0   0   0   0
    0   0   0   0
    
    0   0   0   0
    0   0   0   0
    20   0   0   0
    0   0   0   0
    0   0   0   0
    ];

% Transforms w,l and d into px py and d_or
for i=1:length(w)
    [px(i,:),py(i,:),d_or]=write_Pontos(w,l,d,d_or,cortes_y(i,:),i);
end

% Todas as partículas tem mesmo tamanho w l e th
% Abaixo exemplo pra montar uma majority com os deslocamentos
% d_or=[0     0     0
%     55     110     0
%     55     -110     0
%     55      0     0
%     110     0     0
%     165     0     0
%     220     0     0];

d_or=[
    0   0   0
    60  0   0
    120 0   0
    0   220 0
    60  220 0
    120 220 0
    120 110 0
    180 110 0
    240 110 0
    
    0   440   0
    60  440   0
    120 440   0
    0   660 0
    60  660 0
    120 660 0
    120 550 0
    180 550 0
    240 550 0
    
    240 220 0
    240 440 0
    240 330 0
    300 330 0
    360 330 0
    ];



%% Compute Tensores
compute_NCND=1; % If TRUE computes the tensors again
compute_PAR=1; % If TRUE uses paralel parfor to compute coupling tensor

if compute_NCND
    [Nc,Nd,V]=compute_Tensores(px,py,d_or,th,part_n,platform,compute_PAR);
else
    warning('Tensores não foram recalculados!');
end
Nd(2,1,:) = Nd(1,2,:);
Nd(3,1,:) = Nd(1,3,:);
Nd(3,2,:) = Nd(2,3,:);


Nc(:,:,3,5:7)=zeros(3,3,1,3);
Nc(:,:,2,5:7)=zeros(3,3,1,3);
Nc(:,:,5:7,2)=zeros(3,3,1,3);
Nc(:,:,5:7,3)=zeros(3,3,1,3);
%% Corrente de Spin
% Define a curva da corrente de spin aplicada
% Utiliza a mesma estrutura que do campo aplicado +info: help campute_Happ2
I_s=0e-3;%A
% corrente de spin normalizada i_s = I_s./(q*gammamu0*mu0*Ms*Ns)
%
Ns = 2*Ms*V/gammamu0/hbar;
is=I_s./(q*gammamu0*mu0*Ms*Ns); % magnitude normalizada da corrente de spin
i_s=zeros(N+1,3,part_n);

%% Campo Aplicado
h_app=zeros(N+1,3,part_n);
a=150e-3; % T

% EXP tem todas as combinações de entradas
exp=[1 1 1
    1 1 -1
    1 -1 1
    1 -1 -1
    -1 1 1
    -1 1 -1
    -1 -1 1
    -1 -1 -1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ex=exp(5,:); % testa a primeira combinação de 8 possibilidades
%%%%%% ^ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:part_n
    if i==1
        s=      [   0   0   0   0   0   0   N/10 %1
            0   ex(1)*a   0   0   ex(1)*a   0   N/10 %2
            0   0   0   0   0   0   N/10 %3
            0   0   0   0   0   0   N/10 %4
            0   0   0   0   0   0   N/10 %5
            0   0   0   0   0   0   N/10 %6
            0   0   0   0   0   0   N/10 %7
            0   0   0   0   0   0   N/10 %8
            0   0   0   0   0   0   N/10 %9
            0   0   0   0   0   0   N/10 %10
            ];
    elseif i==2
        s=      [  0   0   0   0   0   0   N/10 %1
            0   ex(2)*a   0   0   ex(2)*a   0   N/10 %2
            0   0   0   0   0   0   N/10 %3
            0   0   0   0   0   0   N/10 %4
            0   0   0   0   0   0   N/10 %5
            0   0   0   0   0   0   N/10 %6
            0   0   0   0   0   0   N/10 %7
            0   0   0   0   0   0   N/10 %8
            0   0   0   0   0   0   N/10 %9
            0   0   0   0   0   0   N/10 %10
            ];
    elseif i==3
        s=      [  0   0   0   0   0   0   N/10 %1
            0   ex(3)*a   0   0   ex(3)*a   0   N/10 %2
            0   0   0   0   0   0   N/10 %3
            0   0   0   0   0   0   N/10 %4
            0   0   0   0   0   0   N/10 %5
            0   0   0   0   0   0   N/10 %6
            0   0   0   0   0   0   N/10 %7
            0   0   0   0   0   0   N/10 %8
            0   0   0   0   0   0   N/10 %9
            0   0   0   0   0   0   N/10 %10
            ];
    elseif i==4
        s=      [  0   0   0   a   0   0   N/10 %1
            a   0   0   a   0   0   N/10 %2
            a   0   0   0   0   0   N/10 %3
            0   0   0   0   0   0   N/10 %4
            0   0   0   0   0   0   N/10 %5
            0   0   0   0   0   0   N/10 %6
            0   0   0   0   0   0   N/10 %7
            0   0   0   0   0   0   N/10 %8
            0   0   0   0   0   0   N/10 %9
            0   0   0   0   0   0   N/10 %10
            ];
    elseif i<=6
        s=      [  0   0   0   a   0   0   N/10 %1
            a   0   0   a   0   0   N/10 %2
            a   0   0   0   0   0   N/10 %3
            0   0   0   0   0   0   N/10 %4
            0   0   0   0   0   0   N/10 %5
            0   0   0   0   0   0   N/10 %6
            0   0   0   0   0   0   N/10 %7
            0   0   0   0   0   0   N/10 %8
            0   0   0   0   0   0   N/10 %9
            0   0   0   0   0   0   N/10 %10
            ];
    else
        s=      [  0   0   0   a   0   0   N/10 %1
            a   0   0   a   0   0   N/10 %2
            a   0   0   0   0   0   N/10 %3
            0   0   0   0   0   0   N/10 %4
            0   0   0   0   0   0   N/10 %5
            0   0   0   0   0   0   N/10 %6
            0   0   0   0   0   0   N/10 %7
            0   0   0   0   0   0   N/10 %8
            0   0   0   0   0   0   N/10 %9
            0   0   0   0   0   0   N/10 %10
            ];
    end
    h_app(:,:,i)=compute_Happ(N,s); % aplicado
end
%% Campo Térmico
K1=500;
HkMs=2*K1/Ms/mu0/Ms;
% dt (adimensional) --> dt_real = dt/gammamu0*Ms
%sig=sqrt(2*alpha*kbT/mu0/V(1)/dt)/Ms;

sig=sqrt(2*alpha*kbT/gammamu0/mu0/Ms/V(1)/time_step)/Ms;

hT=zeros(N+1,3,part_n);
for j=1:part_n
    hhT=randn(N+1,3);
    %     mod_hT=sqrt(sum(hhT'.^2))';
    %     hhT(:,1)=hhT(:,1)./mod_hT;
    %     hhT(:,2)=hhT(:,2)./mod_hT;
    %     hhT(:,3)=hhT(:,3)./mod_hT;
    %     for i = 1:length(hhT)
    %         hhT(i,:)=mod_hT(i)*hhT(i,:);
    %     end
    hT(:,:,j)=sig*hhT;
end
%% Metodo para solucao numerica
fprintf('\n------------------------------------\n');
fprintf('            Range-Kutta            \n');
fprintf('------------------------------------\n');
fprintf('Plataforma:%s\n',platform);
fprintf('dt real:\t\t%3.3e s\n',time_step);
fprintf('dt linha:\t\t%3.3e\n',dt);
fprintf('tempo_total:\t\t%3.3e s\n',tempo_total);
fprintf('N: \t\t\t%d\n',N);
fprintf('------------------------------------\n');

tic;
hc=zeros(N+1,3,part_n); % inicializa o campo térmico
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the process...'),'keepthis','timestamp');
for i = 1 :N
    progress=i/N*100;
    dispstat(sprintf('Progress %3.2f %%',progress),'timestamp');
    
    %RK4-Method
    for j=1:part_n
        for k=1:part_n
            
            hc(i,:,j)=hc(i,:,j)-squeeze(m(i,:,k))*squeeze(Nc(:,:,j,k));
            
        end
        
        h_eff(i,:,j)=compute_Heff(squeeze(m(i,:,j)), squeeze(h_app(i,:,j)), squeeze(hT(i,:,j))*sqrt(V(1)/V(j)), squeeze(hc(i,:,j)),squeeze(Nd(:,:,j)));
        
        m(i+1,:,j)=rk4(squeeze(m(i,:,j)),squeeze(h_eff(i,:,j)),squeeze(i_s(i,:,j)),dt); % Range-Kuta-4
        
    end
end
dispstat('Finished.','keepprev');
toc;
%% Plot Some Results
close all
t=0:1:N;
t=t*time_step/1e-9; %transforma em ns


%% Cores diferentes para particulas que estao em zonas de clock diferentes
cor=[
    0.4940    0.1840    0.5560  % P1
    0.3010    0.7450    0.9330  % P2
    0.3010    0.7450    0.9330  % P3
    0.4940    0.1840    0.5560  % P4
    0.3010    0.7450    0.9330  % P5
    0.3010    0.7450    0.9330  % P6
    0.6350    0.0780    0.1840  % P7
    0.6350    0.0780    0.1840  % P8
    0.6350    0.0780    0.1840  % P9
    0.4940    0.1840    0.5560  % P10
    0.3010    0.7450    0.9330  % P11
    0.3010    0.7450    0.9330  % P12
    0.4940    0.1840    0.5560  % P13
    0.3010    0.7450    0.9330  % P14
    0.3010    0.7450    0.9330  % P15
    0.6350    0.0780    0.1840  % P16
    0.6350    0.0780    0.1840  % P17
    0.6350    0.0780    0.1840  % P18
    0.6350    0.0780    0.1840  % P19
    0.6350    0.0780    0.1840  % P20
    0.4660    0.6740    0.1880  % P21
    0.4660    0.6740    0.1880  % P22
    0.4660    0.6740    0.1880  % P23
    0.3010    0.7450    0.9330  % P
    0.3010    0.7450    0.9330  % P
    0.3010    0.7450    0.9330  % P
    0.3010    0.7450    0.9330  % P
    0.6350    0.0780    0.1840  % P
    0.6350    0.0780    0.1840  % P
    0.6350    0.0780    0.1840  % P
    0.6350    0.0780    0.1840  % P
    0.6350    0.0780    0.1840  % P
    0.4660    0.6740    0.1880  % P
    0.4660    0.6740    0.1880  % P
    0.4660    0.6740    0.1880  % P
    0.4660    0.6740    0.1880
    0.4940    0.1840    0.5560];

% Plota as magnetizacoes e campos aplicados
figure('Position',[100 100 600 600], ...
    'Name','Magnetizações e Campos Aplicados');

cols=2; %numero de colunas no plot
rows=ceil(part_n/cols); % numero de linhas


for j=1:part_n
    subplot(rows,cols,j);
    plot(t,squeeze(m(:,1:2,j))); % Plota a Magnetização
    hold on
    plot(t,squeeze(h_app(:,1:2,j))/a,'--'); % Plota o Campo aplicado em X (1) normalizado por a
    ylim([-1.2 1.2]);
    if j==part_n
        title([' Output P_{' num2str(j) '}'], 'color', cor(1,:)); % cosiderando como ultima particula
    else
        title(['P_{' num2str(j) '}'], 'color', cor(1,:));
    end
    hl=legend('m_x','m_y','h_{app_x}','h_{app_y}');
    set(hl,'Orientation','Horizontal','Location','Best')
end


%% Plota as particulas
figure('Position',[800 100 600 600], ...
    'Name','Alocação das Partículas');

espaco1=1:2:2*part_n-1;
espaco2=0:1:part_n-1;
p_a=espaco1*25+espaco2*d_min;

for i=1:part_n
    fill(px(i,:)+d_or(i,1)+px(1,2),py(i,:)+py(1,1)+d_or(i,2),cor(i,:),'EdgeColor','none','LineStyle','none')
    hold on
    text(d_or(i,1)+sum(abs(px(i,1:2)))/2,d_or(i,2)+py(i,2),['P_{' num2str(i) '}'],'fontsize',10,'FontWeight','bold','Color','black','HorizontalAlignment','Center');
end


daspect([1 1 1])
title(['Distâncias ' num2str(d_min) 'nm'])
ylim([min(d_or(:,2))-25 max(d_or(:,2))+125])
xlim([min(d_or(:,1))-25 max(d_or(:,1))+75])
xlabel('nm')
ylabel('nm')

% figure
% stem(1:part_n,round(squeeze(m(5001,2,:))))