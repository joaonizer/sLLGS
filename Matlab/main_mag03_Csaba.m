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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part_n=15; % quantidade de particulas
%%%%%% ^ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ^ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm_time=[1 10 200 2000]; %ns
n_runs = 50; % numero de simulações em cada tipo de rampa
comp=zeros(part_n,30,4);
for tt=1:length(rm_time)
    for jj =1:n_runs
        
        
        N=round((rm_time(tt)+7.5)/10e-3);
        %N = 5000;       % numero de passos
        tempo_total=N*10e-12;% Tempo total de simulação
        alpha=1.0;
        n = [0 1 0];
        T=300;        % Kelvin
        Ms=800e3;   % A/m
        kbT=kb*T;   % J
        ti = 0;     % instante inicial da variavel independente
        %[N,tempo_total,ti,tf,dt]=compute_Time(gammamu0,Ms,N,tempo_total,ti);
        dt=tempo_total/N*gammamu0*Ms;
        time_step=dt/(gammamu0*Ms); % segundos
        if time_step>1e-12*alpha
            warning('Time_Step muito grande! Considere aumentar N!')
        end
        %% Configuracoes do Sistema
        
        
        
        d_min=10; % distancia mínima entre as partículas
        
        
        alpha_l=1/(1+alpha^2);
        mi = [0 1 0]/sqrt(1); % valor inicial das variaveis dependente
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
        
        % Transforms w,l and d into px py and d_or
        for i=1:length(w)
            [px(i,:),py(i,:),d_or]=write_Pontos(w,l,d,d_or,cortes_y(i,:),i);
        end
        
        
        
        %% Compute Tensores
        if (tt==1 && jj==1)
            compute_NCND=1; % If TRUE computes the tensors again
        else
            compute_NCND=0;
        end
        compute_PAR=0; % If TRUE uses paralel parfor to compute coupling tensor
        
        if compute_NCND
            [Nc,Nd,V]=compute_Tensores(px,py,d_or,th,part_n,platform,compute_PAR);
        else
            warning('Tensores não foram recalculados!');
        end
        
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
        
        for i=2:part_n
            if i>1
                s=      [
                    a   0   0   a   0   0   500 %2
                    a   0   0   0   0   0   rm_time(tt)/10e-3 %3
                    0   0   0   0   0   0   250 %4
                    
                    ];
            end
            h_app(:,:,i)=compute_Happ(N,s); % aplicado
        end
        %% Campo Térmico
        K1=500;
        HkMs=2*K1/Ms/mu0/Ms;
        % dt (adimensional) --> dt_real = dt/gammamu0*Ms
        %sig=sqrt(2*alpha*kbT/mu0/V(1)/dt)/Ms;
        seed=length(rm_time)*tt + jj+50;
        rng(seed);
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
        clc
        fprintf('\n--------  %d/4 -- %d/30  ---------\n',tt,jj);
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
        
        mask = [1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1];% -1 1 -1 1 -1 1 -1 1 -1 1];
        comp(:,jj,tt) = round(squeeze(m(N,2,:)))'==mask;
    end
    %% Final analisys
error_count=zeros(part_n,4);
for tt=1:length(rm_time)
    for jj=1:n_runs
        error_place=find(comp(:,jj,tt)==0);
        if ~isempty(error_place)
            error_count(error_place(1):end,tt)=error_count(error_place(1):end,tt)+1;
        end
    end
end
figure('Position',[100 100 1000 500], ...
    'Name','Alocação das Partículas');
plot(1:15,error_count(:,:)/100,'o-');

lh=legend('1 ns', '10 ns', '200 ns', '2 $\mu$s')
set(lh,'Location','Best');
%xlabel('Particle')
title('100 Runs')
ylabel('P_{error}')
grid on
sdf('P1')
print('-djpeg','-r300',['P_error.jpeg'])
end

