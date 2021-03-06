%% Inicio
clc
close all
%clear all
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
N = 10000;       % numero de passos
tempo_total=200e-9;% Tempo total de simulação
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
part_n=26; % quantidade de particulas
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
    0     0     0     0
    0     0     0     0
    0     0     0     0
    0     0     0     0
    0     0     0     0
    0     0     0     0
    0     0     0     0
    0     0     0     0
    0     0     0     0
    0     0     0     0
    0    0     0     0
    0     0     0     0
    0    0     0     0
    0     0     0     0
    0    0     0     0
    0     0     0     0
    0    0     0     0
    0     0     0     0
    0    0     0     0
    0     0     0     0
    0     0     0   -10
    10     0     0    0
    0     0     0     0
    0     0     0     0
    0     0     0     0
    0     0     0     0
    ];
%reverse =[1  5  9  2  6 10 13 15 17  3  7 11  4  8 12 14 16 18 19 20 21 22 23];
%cortes_y(reverse,:)=cortes_y(1:23,:);

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

dy=124*[1:6]; %% deslocamentos em y
dx=60*[0:8]; %% deslocamento em x
offset=0;
d_or=[
    dx(1)     0     0 % P1
    dx(1)   dy(2)     0 % P2
    dx(1)   dy(4)+offset     0 % P3
    dx(1)   dy(6)+offset     0 % P4
    
    dx(2)     0     0 % P1
    dx(2)   dy(2)     0 % P2
    dx(2)   dy(4)+offset     0 % P3
    dx(2)   dy(6)+offset     0 % P4
    
    dx(3)     0     0 % P1
    dx(3)   dy(2)     0 % P2
    dx(3)   dy(4)+offset     0 % P3
    dx(3)   dy(6)+offset     0 % P4
    
    dx(4)     0     0 % P1
    dx(4)   dy(2)     0 % P2
    dx(4)   dy(4)+offset     0 % P3
    dx(4)   dy(6)+offset     0 % P4
    
    dx(5)     0     0 % P1
    dx(5)   dy(2)     0 % P2
    dx(5)   dy(4)+offset     0 % P3
    dx(5)   dy(6)+offset     0 % P4
    
    
    dx(5)   dy(1)     0 % P13
    dx(5)   dy(5)+offset     0 % P14
    dx(6)   dy(1)     0 % P15
    dx(6)   dy(5)+offset     0 % P16
    dx(7)   dy(1)     0 % P17
    dx(7)   dy(5)+offset     0 % P18
    dx(7)   dy(2)     0 % P19
    dx(7)   dy(4)     0 % P20
    dx(7)   dy(3)     0 % P21
    dx(8)   dy(3)     0 % P22
    dx(9)   dy(3)     0 % P23
    ];

%d_or(13:end,1)=d_or(13:end,1)+25;
%d_or(21:end,1)=d_or(21:end,1)+25;


%% Compute Tensores
compute_NCND=1; % If TRUE computes the tensors again
compute_PAR=1; % If TRUE uses paralel parfor to compute coupling tensor
imprimir=1;
if compute_NCND
    [Nc,Nd,V]=compute_Tensores(px,py,d_or,th,part_n,platform,compute_PAR);
    Nc_old = Nc;
else
    warning('Tensores não foram recalculados!');
end

results=[1 -1 -1 1 -1 1 1 -1 1 -1 -1 1 -1 -1 1 1 -1 -1 -1 -1 -1 1 -1
    -1 -1 1 1 1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 -1 1 -1 1 1 -1 1
    1 1 -1 -1 -1 -1 1 1 1 1 -1 -1 1 -1 -1 1 1 -1 1 -1 1 -1 1
    -1 1 1 -1 1 -1 -1 1 -1 1 1 -1 -1 -1 1 1 -1 -1 -1 -1 -1 1 -1];

results_and_or = [
%    1	1	1	1	-1	-1	-1	-1	1	1	1	1	1	1	-1	-1	1	1
%    -1	1	-1	1	1	-1	1	-1	-1	1	-1	1	-1	1	1	-1	-1	1
%    1	-1	1	-1	-1	1	-1	1	1	-1	1	-1	-1	1	1	-1	-1	1
%    -1	-1	-1	-1	1	1	1	1	-1	-1	-1	-1	-1	-1	1	1	-1	-1
1	1	1	1	-1	-1	-1	-1	1	1	1	1	-1	-1	-1	-1	1	1	1	1	1	1	-1	-1	1	1
-1	1	-1	1	1	-1	1	-1	-1	1	-1	1	1	-1	1	-1	-1	1	-1	1	-1	1	1	-1	-1	1
1	-1	1	-1	-1	1	-1	1	1	-1	1	-1	-1	1	-1	1	1	-1	1	-1	-1	1	1	-1	-1	1
-1	-1	-1	-1	1	1	1	1	-1	-1	-1	-1	1	1	1	1	-1	-1	-1	-1	-1	-1	1	1	-1	-1
];


analysis=zeros(10,4);
compare=zeros(10,4,part_n);
for iiii=5:5
    reduce=1-0.1*(iiii-1); % percentual
    %reduce=0.4;
    Nc(:,:,23,17:18)=Nc_old(:,:,23,17:18)*reduce;
    Nc(:,:,17:18,23)=Nc_old(:,:,17:18,23)*reduce;
    Nc(:,:,24,19:20)=Nc_old(:,:,24,19:20)*reduce;
    Nc(:,:,19:20,24)=Nc_old(:,:,19:20,24)*reduce;
    
    Nc(:,:,21,13:14)=Nc_old(:,:,21,13:14)*reduce;
    Nc(:,:,13:14,21)=Nc_old(:,:,13:14,21)*reduce;
    Nc(:,:,22,15:16)=Nc_old(:,:,22,15:16)*reduce;
    Nc(:,:,15:16,22)=Nc_old(:,:,15:16,22)*reduce;
    %Nc(:,:,22,19:20)=Nc_old(:,:,22,19:20)*reduce;
    %Nc(:,:,19:20,22)=Nc_old(:,:,19:20,22)*reduce;
    %Nc(:,:,15,19)=Nc_old(:,:,15,19)*reduce;
    %Nc(:,:,19:15)=Nc_old(:,:,19:15)*reduce;
    %Nc(:,:,16,20)=Nc_old(:,:,16,20)*reduce;
    %Nc(:,:,20,16)=Nc_old(:,:,20,16)*reduce;
    
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
    
    
    for jj=1:4
        fprintf('%d/%d %d/%d',iiii,10,jj,4);
        h_app=zeros(N+1,3,part_n);
        a=200e-3; % T
        
        
        % EXP tem todas as combinações de entradas
        %    X  Y
        exp=[1  1
            1 -1
            -1  1
            -1 -1
            ];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ex=exp(jj,:); % testa a primeira combinação de 8 possibilidades
        %%%%%% ^ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        colors = [
            0.4940    0.1840    0.5560  % Roxo
            0.3010    0.7450    0.9330  % Azul
            0.6350    0.0780    0.1840  % Vermelho
            0.4660    0.6740    0.1880  % Verde
            0.4940    0.1840    0.5560  % Roxo
            ];
        for i=1:part_n
            if i==1
                cor(i,:)=colors(1,:);
                s=      [
                    0   ex(2)*a   0   0   ex(2)*a   0   N/10 %2
                    0   0   0   0   0   0   N/10 %3
                    0   0   0   0   0   0   N/10 %4
                    0   0   0   0   0   0   N/10 %5
                    0   0   0   0   0   0   N/10 %6
                    0   0   0   0   0   0   N/10 %7
                    0   0   0   0   0   0   N/10 %8
                    0   0   0   0   0   0   N/10 %9
                    0   0   0   0   0   0   N/10 %10
                    0   0   0   0   0   0   N/10 %1
                    ];
            elseif i==2
                cor(i,:)=colors(1,:);
                s=      [
                    0   ex(1)*a   0   0   ex(1)*a   0   N/10 %2
                    0   0   0   0   0   0   N/10 %3
                    0   0   0   0   0   0   N/10 %4
                    0   0   0   0   0   0   N/10 %5
                    0   0   0   0   0   0   N/10 %6
                    0   0   0   0   0   0   N/10 %7
                    0   0   0   0   0   0   N/10 %8
                    0   0   0   0   0   0   N/10 %9
                    0   0   0   0   0   0   N/10 %10
                    0   0   0   0   0   0   N/10 %1
                    ];
            elseif i==3
                cor(i,:)=colors(1,:);
                s=      [
                    0   ex(2)*a   0   0   ex(2)*a   0   N/10 %2
                    0   0   0   0   0   0   N/10 %3
                    0   0   0   0   0   0   N/10 %4
                    0   0   0   0   0   0   N/10 %5
                    0   0   0   0   0   0   N/10 %6
                    0   0   0   0   0   0   N/10 %7
                    0   0   0   0   0   0   N/10 %8
                    0   0   0   0   0   0   N/10 %9
                    0   0   0   0   0   0   N/10 %10
                    0   0   0   0   0   0   N/10 %1
                    ];
            elseif i==4
                cor(i,:)=colors(1,:);
                s=      [
                    0   ex(1)*a   0   0   ex(1)*a   0   N/10 %2
                    0   0   0   0   0   0   N/10 %3
                    0   0   0   0   0   0   N/10 %4
                    0   0   0   0   0   0   N/10 %5
                    0   0   0   0   0   0   N/10 %6
                    0   0   0   0   0   0   N/10 %7
                    0   0   0   0   0   0   N/10 %8
                    0   0   0   0   0   0   N/10 %9
                    0   0   0   0   0   0   N/10 %10
                    0   0   0   0   0   0   N/10 %1
                    ];
            elseif i<=16
                cor(i,:)=colors(2,:);
                s=      [
                    0   0   0   a   0   0   N/10 %1
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
            elseif i<=26
                cor(i,:)=colors(3,:);
                s=      [
                    0   0   0   0   0   0   N/10 %10
                    0   0   0   a   0   0   N/10 %1
                    a   0   0   a   0   0   N/10 %2
                    a   0   0   0   0   0   N/10 %3
                    0   0   0   0   0   0   N/10 %4
                    0   0   0   0   0   0   N/10 %5
                    0   0   0   0   0   0   N/10 %6
                    0   0   0   0   0   0   N/10 %7
                    0   0   0   0   0   0   N/10 %8
                    0   0   0   0   0   0   N/10 %9
                    
                    ];
            elseif i<=21
                cor(i,:)=colors(4,:);
                s=  [
                    0   0   0   0   0   0   N/10 %9
                    0   0   0   0   0   0   N/10 %10
                    0   0   0   a   0   0   N/10 %1
                    a   0   0   a   0   0   N/10 %2
                    a   0   0   0   0   0   N/10 %3
                    0   0   0   0   0   0   N/10 %4
                    0   0   0   0   0   0   N/10 %5
                    0   0   0   0   0   0   N/10 %6
                    0   0   0   0   0   0   N/10 %7
                    0   0   0   0   0   0   N/10 %8
                    ];
            else
                cor(i,:)=colors(5,:);
                s=  [
                    0   0   0   0   0   0   N/10 %9
                    0   0   0   0   0   0   N/10 %10
                    0   0   0   0   0   0   N/10 %1
                    0   0   0   a   0   0   N/10 %2
                    a   0   0   a   0   0   N/10 %3
                    a   0   0   0   0   0   N/10 %4
                    0   0   0   0   0   0   N/10 %5
                    0   0   0   0   0   0   N/10 %6
                    0   0   0   0   0   0   N/10 %7
                    0   0   0   0   0   0   N/10 %8
                    ];
            end
            h_app(:,:,i)=compute_Happ(N,s); % aplicado
        end
        % Campo Térmico
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
        %     %% Plot Some Results
        close all
        t=0:1:N;
        t=t*time_step/1e-9; %transforma em ns
        set(groot, 'defaultAxesTickLabelInterpreter','latex');
        set(groot, 'defaultLegendInterpreter','latex');
        
        %% Plota as magnetizacoes e campos aplicados
        if imprimir
            figure('Position',[0 0 1400 1400], ...
                'Name','Magnetizações e Campos Aplicados');
            
            cols=7; %numero de colunas no plot
            rows=7;%ceil(part_n/cols); % numero de linhas
            plot_place=[
                43, 29, 15,  1, ...
                44, 30, 16,  2, ...
                45, 31, 17,  3, ...
                46, 32, 18,  4, ...
                47, 33, 19,  5, ...
                40, 12, 41,13, 42, 14];
            
            for j=1:part_n
                subplot(rows,cols,plot_place(j));
                plot(t,squeeze(m(:,1:2,j))); % Plota a Magnetização
                hold on
                plot(t,squeeze(h_app(:,1:2,j))/a,'--'); % Plota o Campo aplicado em X (1) normalizado por a
                ylim([-1.2 1.2]);
                if j==part_n
                    title(['$Output\ P_{' num2str(j) '}$'], 'color', cor(j,:),'Interpreter','latex'); % cosiderando como ultima particula
                elseif j==1
                    title(['$P_{' num2str(j) '} - Y$'], 'color', cor(j,:),'Interpreter','latex');
                elseif j==2
                    title(['$P_{' num2str(j) '} - {X}$'], 'color', cor(j,:),'Interpreter','latex');
                elseif j==3
                    title(['$P_{' num2str(j) '} - {Y}$'], 'color', cor(j,:),'Interpreter','latex');
                elseif j==4
                    title(['$P_{' num2str(j) '} - X$'], 'color', cor(j,:),'Interpreter','latex');
                else
                    title(['$P_{' num2str(j) '}$'], 'color', cor(j,:),'Interpreter','latex');
                end
                %hl=legend('m_x','m_y','h_{app_x}','h_{app_y}');
                %set(hl,'Orientation','Horizontal','Location','Best')
            end
            sdf('P1');
            print('-dbmp','-r150',['XOR_' num2str(jj) '.bmp'])
            close all
        end
        m_response=round(squeeze(m(end,2,:)));
        compare(iiii,jj,:)=m_response==results_and_or(jj,1:part_n)';
        if sum(m_response==results_and_or(jj,1:part_n)')==part_n
            analysis(iiii,jj)=1
        else
            analysis(iiii,jj)=0
        end
    end
    
    
end

%% Plota as particulas
if imprimir
    figure('Position',[800 100 300 600], ...
        'Name','Alocação das Partículas');
    
    espaco1=1:2:2*part_n-1;
    espaco2=0:1:part_n-1;
    p_a=espaco1*25+espaco2*d_min;
    
    for i=1:part_n
        fill(px(i,:)+px(1,2)+d_or(i,1),py(i,:)+py(1,2)+d_or(i,2),cor(i,:),'EdgeColor','none','LineStyle','none')
        hold on
        text(d_or(i,1)+sum(abs(px(i,1:2)))/2,d_or(i,2)+py(i,2),['$P_{' num2str(i) '}$'],'Interpreter','latex','fontsize',9,'FontWeight','bold','Color','black','HorizontalAlignment','Center');
    end
    
    
    daspect([1 1 1])
    title(['XOR Architecture'],'Interpreter','latex')
    ylim([min(d_or(:,2))-25 max(d_or(:,2))+125])
    xlim([min(d_or(:,1))-25 max(d_or(:,1))+75])
    xlabel('$nm$','Interpreter','latex')
    ylabel('$nm$','Interpreter','latex')
    sdf('P1');
    print('-dbmp','-r150',['XOR_Architecture.bmp'])
end
% % figure
% % stem(1:part_n,round(squeeze(m(5001,2,:))))