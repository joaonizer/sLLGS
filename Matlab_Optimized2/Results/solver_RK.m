function [m]=solver_RK(platform, time_step,dt,tempo_total,N,T,alpha,...
                        part_n, n,jj,Nd,Nc,h_eff,h_app,...
                        i_s,v,dW,HkMs,m)    
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
    fprintf('Experiment No: \t%d\n',jj);
    fprintf('------------------------------------\n');
    dispstat('','init'); % One time only initialization
    dispstat(sprintf('Begining the process...'),'keepthis','timestamp');
    tic;
    hc=zeros(N+1,3,part_n); % inicializa o campo tÃ©rmico
    hd=hc;
    n_part_n=ones(part_n,1)*n;
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
    toc;
    dispstat('Finished.','keepprev');