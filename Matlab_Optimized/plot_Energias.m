function plot_Energias(m, hT, h_app, hc, V,t,data,cell_size,im_id,Wi,Li)
% Calcula e plota as energias:
%   -Térmica
%   -Dissipada
%   -Desmagnetização
%   -Acoplamento

fprintf('Imprimindo resultados:')
global Ms q time_step alpha gammamu0 K1 T;
mu0=4*pi*1e-7;
%% Energia Térmica
% Eth(t) = -mu0*V*Ms*(mx(t)*Hthx(t) + my(t)*Hthy(t) + mz(t)*Hthz(t))*6.242e+18
Eth=-mu0*V*Ms^2*m.*hT/q;
Eth= sum(Eth,2);
%% Energia Zeeman
% E0(t)  = -mu0*V*Ms*(mx(t)*Hx(t))*6.242e+18
E0=-mu0*V*Ms^2*m(:,1).*h_app(:,1)/q;
%% Energia Desmagnetização
% Ed(t)  = -mu0*V*Ms*(mx(t)*Hdx(t) + my(t)*Hdy(t) + mz(t)*Hdz(t))*6.242e+18,
% Onde:
% Hdx(t) = -Ms(Nxx*mx(t) + Nxy*my(t) + Nxz*mz(t))
% Hdy(t) = -Ms(Nyx*mx(t) + Nyy*my(t) + Nyz*mz(t))
% Hdz(t) = -Ms(Nzx*mx(t) + Nzy*my(t) + Nzz*mz(t)),
% onde Nxx, Nxy etc. são as componentes do teNcor desmagnetização
hd(:,1)=Nc(1,1)*m(:,1)+Nc(1,2)*m(:,2)+Nc(1,3)*m(:,3);
hd(:,2)=Nc(2,1)*m(:,1)+Nc(2,2)*m(:,2)+Nc(2,3)*m(:,3);
hd(:,3)=Nc(3,1)*m(:,1)+Nc(3,2)*m(:,2)+Nc(3,3)*m(:,3);
Ed=mu0*V*Ms^2*hd.*m/q/2;
Ed=sum(Ed,2);
%% Energia de Acoplamento
% Ec(t)  = -mu0*V*Ms*(mx(t)*Hcx + my(t)*Hcy + mz(t)*Hcz)*6.242e+18
Ec=-mu0*V*Ms^2*m.*hc/q;
Ec= sum(Ec,2);


%% Plots

% Plot the Results
figure('Visible','off',...
    'PaperPosition',[0 0 11.692 8.267]*2,...
    'PaperSize',[11.692 8.267]*2)
subplot(3,2,1);
plot(data(:,1)*time_step/1e-9,data(:,2:4),'LineWidth',2);
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,m(:,:),'--','LineWidth',2);
grid on;
title({['\alpha=' num2str(alpha) '|\gamma=' num2str(gammamu0)...
    '|dt=' num2str(time_step) 's'];['K1 = ' num2str(K1) '|' num2str(T) 'K'...
    '|CS=' num2str(cell_size(1)/1e-9) 'x' num2str(cell_size(2)/1e-9) 'x' num2str(cell_size(3)/1e-9) 'nm3']})
%ylim([-1.8 1.2]);
ylabel('m');
xlabel('Tempo (Nc)')
xlim([0 max(t)]);
ylim([-1.6 1.2]);
set(gca,'FontSize',24);
% Legendas auxiliares
h = zeros(5, 1);
ax.ColorOrderIndex = 1;
h(1) = plot(NaN,NaN,'LineWidth',2);
h(2) = plot(NaN,NaN,'LineWidth',2);
h(3) = plot(NaN,NaN,'LineWidth',2);
ax.ColorOrderIndex = 1;
h(4) = plot(NaN,NaN,'LineWidth',2);
ax.ColorOrderIndex = 1;
h(5) = plot(NaN,NaN,'--','LineWidth',2);
lh=legend(h, 'Mx','My','Mz','OOMMF','RK');
set(lh,'Orientation','horizontal','Location','south','FontSize',12);

% Plota a Particula com iNcet
axes('Position',[.05 .8 .05 .05])
box on
imshow(['./OOMMF_sim/' im_id '.gif'])
xlabel(['W = ' num2str(Wi) 'nm']);
ylabel(['L = ' num2str(Li) 'nm']);
set(gca,'FontSize',10);
%% Energia Térmica
subplot(3,2,6);
plot(t,Eth,'LineWidth',1);
title('Energia Térmica')
ylabel('eV');
xlabel('Tempo (Nc)')
set(gca,'FontSize',24);
%% Desmagnetizacao Alternativo
% subplot(3,2,4)
% plot(t,Ud_oommf,'LineWidth',2);
% ax = gca;
% ax.ColorOrderIndex = 3;
% hold on
% plot(t,Ud,'LineWidth',2);
% legend('OOMMF','RK');
% title('Energia Desmagnetização')
% ylabel('eV');
% xlabel('Tempo (Nc)')
% ylim([-5 max(max(Ed),max(data(:,8)/q))+5]);
% set(gca,'FontSize',24);
% subplot(3,2,6)
% plot(theta2*180/pi,Ud_oommf,'*');
% hold on
% ax = gca;
% ax.ColorOrderIndex = 3;
% plot(theta1*180/pi,Ud,'*');
% legend('OOMMF','RK');
% xlim([0 360]);
% title('Energia Desmagnetização')
% ylabel('eV');
% xlabel('\theta (Graus)')
% set(gca,'FontSize',24);
%% Demagnetização
subplot(3,2,3);
%plot(data(:,1)*time_step/1e-9,data(:,8)/q,'LineWidth',2);
%hold on
plot(t,data(:,8)/q,'LineWidth',2);
hold on
%plot(data(:,1)*time_step/1e-9,Ed_oommf,'LineWidth',2);
plot(t,Ed,'LineWidth',2);
title('Energia Desmagnetização')
ylabel('eV');
xlabel('Tempo (Nc)')
ylim([-5 max(max(Ed),max(data(:,8)/q))+5]);
set(gca,'FontSize',24);
%lh=legend('OOMMF','OOMMF w Eq','RK');
lh=legend('OOMMF','RK');
set(lh,'Location','South','Orientation','horizontal','FontSize',12);
%% Zeeman Energy
subplot(3,2,4);
plot(t,data(:,10)/q,'LineWidth',2);
hold on
plot(t,E0,'--','LineWidth',2);
title('Energia Zeeman')
ylabel('eV');
xlabel('Tempo (Nc)')
set(gca,'FontSize',24);
lh=legend('OOMMF','RK');
set(lh,'FontSize',12);
%% Desmagnetização Theta
subplot(3,2,5);
theta_RK=cart2sph(m(:,1),m(:,2),m(:,3))*180/pi;
theta_OOMMF=cart2sph(m_oommf(:,1),m_oommf(:,2),m_oommf(:,3))*180/pi;
theta_RK(theta_RK<0)=theta_RK(theta_RK<0)+360;
theta_OOMMF(theta_OOMMF<0)=theta_OOMMF(theta_OOMMF<0)+360;
plot(theta_OOMMF,data(:,8)/q,'*','LineWidth',2);
hold on
%plot(theta_OOMMF,Ed_oommf,'*','LineWidth',2);
plot(theta_RK,Ed,'*','LineWidth',2);
title('Energia Desmagnetização')
ylabel('eV');
xlabel('\theta (Graus)')
ylim([-5 max(max(Ed),max(data(:,8)/q))+5]);
set(gca,'FontSize',24);
%lh=legend('OOMMF','OOMMF--RK','RK');
lh=legend('OOMMF','RK');
set(lh,'Orientation','horizontal','Location','south','FontSize',12);
%% Acoplamento
subplot(3,2,2);
plot(t,Ec,'LineWidth',2);
title('Energia Acoplamento')
ylabel('eV');
xlabel('Tempo (Nc)')
set(gca,'FontSize',24);
%% Salva o PDF
print(gcf,'-dpdf','-r300',['./img/HxM_' im_id '.pdf'])
fprintf('\tOK!\n');
end