%% Inicio
clc
close all
%more off
%clear all
%clear all
global alpha alpha_l Ms K1 HkMs sig kbT q time_step gammamu0 A n mu0;
%% Constantes
%Ku=0.26*1; %eV
q = 1.60217662e-19; % carga do eletron C
mu0=4*pi*1e-7;      % H/m ou T.m/A
kb=1.38064852e-23;  % m2.kg.s-2.K-1
A=13e-12;           % Exchange Stiffness J/m
t2am=7.943e5;       % converte T para A/m
gammamu0=2.211e5;   % (sA/m)-1
hbar=2.05457e-34; % J.s/rad  -> h/2pi
%% Configuracoes do Algoritmo
if strcmp(computer,'GLNXA64')
    platform= 'lin';
else
    platform = 'win';
end

N = 10000;       % numero de passos
tempo_total=1e-9;% Tempo total de simulaÃ§Ã£o
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
name=['comp_slanted_up_down-' num2str(T) 'K-' num2str(N) 'steps-' num2str(tempo_total*1e9) 'ns-' num2str(alpha*100) 'alpha'];
grid=[
    1 0 0
    0 2 0
    0 0 3
    %               0   0   1   1 1
    %               0   0   1   0 0
    %               1   1   1   0 0
    %               0   0   1   0 0
    %               0   0   1   1 1
    
    %     1   1   1 0  0   0   0
    %     0   0   2 1  1   0   0
    %     1   1   1 0  1   0   0
    %     0   0   0 0  1   0   0
    %     0   0   0 0  3   1   1
    %     0   0   0 0  1   0   0
    %     1   1   1 0  1   0   0
    %     0   0   2 1  1   0   0
    %     1   1   1 0  0   0   0
    ];
part_n=sum(sum(grid>0)); % quantidade de particulas

d_min=15; % distancia mÃ­nima entre as particulas


alpha_l=1/(1+alpha^2);
mi = [0 1 0]/sqrt(1); % valor inicial das variaveis dependente
m=zeros(N+1,3,part_n);
h_eff=zeros(N,3,part_n);
m(1,1,1)=sqrt(1-2*0.1^2);
m(1,2,1)=0.1;
m(1,3,1)=0.1;
for i=2:part_n % inicializa as partÃ­culas de forma antiferromagnetica
    m(:,:,i)=(-1)^(i-1)*m(:,:,1);
end
%% Dimensoes da Particula

w=ones(1,part_n)*50;  % width of particles

l=ones(1,part_n)*100; % length of particles

th=ones(1,part_n)*10;   %thickness of particles


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
            cortes_y(count,:)=[0 20 0 0];
            count=count+1;
        elseif grid(j,i)==3 %or
            cortes_y(count,:)=[0 20 -20 0];
            count=count+1;
        end
    end
end
% Transforms w,l and d into px py and d_or
for i=1:length(w)
    [px(i,:),py(i,:)]=write_Pontos(w,l,cortes_y(i,:),i);
end
dx=(w(1)+100)*[0:50];
dy=(l(1)+250)*[0:50]; %% deslocamentos em y
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
NMC=[0.5 1 5 10]*1e6;
%% Compute Tensores
Nd=zeros(length(NMC),1000,3*part_n,3*part_n);
compute_PAR=1;
for j=1:length(NMC)
    for i=1:1000
        fprintf('%d/100 %d',i,j)
        [Nd(j,i,:,:), V]=compute_Nd_NMC(px,py,th,part_n,compute_PAR,platform,NMC(j));
        pause(1)
    end
end

%% Plots
cor=[0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    ];
figure('Position',[0 0 1920 1080])
for j=1:9
    subplot(3,3,j)
    for i=1:4
        %subplot(2,2,i)
        %pd(j,i,:) = fitdist(squeeze(Nd(i,:,j,j))','Normal');
        %x_values=pd(j,i,:).mu-4*pd(j,i,:).sigma:0.0001:pd(j,i,:).mu+4*pd(j,i,:).sigma;
        %y=pdf(pd(j,i,:),x_values);
        if ~j==1
            h=histfit(squeeze(Nd(i,abs(Nd(i,:,j,j))<.5,j,j)));
        else
            h=histfit(squeeze(Nd(i,abs(Nd(i,:,j,j))<1,j,j)));
        end
        delete(h(1))
        set(h(2),'color',cor(i,:))
        %plot(x_values,y)
        hold on
        %title(['Coeficiente zz - NMC=' num2str(NMC(i)/1e6) '\times 10^6'])
        %xlim([min(min(squeeze(Nd(:,:,dim,dim)))),max(max(squeeze(Nd(:,:,dim,dim))))])
        %ylim([0 70])
        
        
        
    end
    
    if sum(j==[1 2 3])
        particle='Retangular';
    elseif sum(j==[4 5 6])
        particle='Corte Sup 20nm';
    else
        particle='Corte Sup/Inf 20nm';
    end
    
    if sum(j==[1 4 7])
        index='_{xx}';
    elseif sum(j==[2 5 8])
        index='_{yy}';
    else
        index='_{zz}';
    end
    lh=legend('0.5M ', '1M', '5M', '10M');
    set(lh,'Location','Best')
    title([particle ' Nd' index])
    %sdf('P1')
    set(lh,'FontSize',12)
    %set(gca,'LooseInset',get(gca,'TightInset')+0.01);
end
%

%% Plots Box plot
cor=[0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    ];
figure('Position',[0 0 1920 1080])
Nd_box=Nd;
Nd_box(Nd_box>1)=NaN;
for j=1:9
    
    if sum(j==[1 2 3])
        particle='Retangular 50 \times 100';
    elseif sum(j==[4 5 6])
        particle='Corte Sup 20nm';
    else
        particle='Corte Sup/Inf 20nm';
    end
    
    samples={'0.5M', '1M', '5M', '10M'};
    
    if sum(j==[1 4 7])
        index='_{xx}';
    elseif sum(j==[2 5 8])
        index='_{yy}';
    else
        index='_{zz}';
    end
    
    subplot(3,3,j)
    boxplot(squeeze(Nd_box(:,:,j,j))',samples,'Notch','on');
    if j>6
    xlabel('Número de Amostras')
    end
    ylabel(['Nd' index],'interpreter','tex')
    title(particle);    
        
    
    
    
    
    %set(lh,'Location','Best')
    %title([particle ' Nd' index])
    sdf('P1')
    %set(lh,'FontSize',12)
    set(gca,'LooseInset',get(gca,'TightInset')+0.01);
end
%
