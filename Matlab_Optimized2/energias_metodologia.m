%%
% linkado com main mag wire
figure;
plot(t,squeeze(m(:,:,1)));
hold on
set(gca,'ColorOrderIndex',1);
plot(t,squeeze(m(:,:,2)),'--');
xlabel('Tempo (ns)')
ylabel('m')
legend('m_{x_1}','m_{y_1}','m_{z_1}','m_{x_2}','m_{y_2}','m_{z_2}');
ylim([-1.2 1.2]);
sdf('P1');
%%
%subplot(1,3,1)
Ud1=-0.5*mu0*V(1)*Ms^2*sum(squeeze(m(1:end-1,:,1)).*squeeze(hd(1:end-1,:,1)),2)/abs(q);
Ud2=-0.5*mu0*V(2)*Ms^2*sum(squeeze(m(1:end-1,:,2)).*squeeze(hd(1:end-1,:,2)),2)/abs(q);
plot(t(1:end-1),Ud1);
hold on
plot(t(1:end-1),Ud2);
xlabel('Tempo (ns)')
ylabel('U_d (eV)')
legend('U_{d_1}','U_{d_2}');
sdf('P1');
%%
%subplot(1,3,2)
Uc1=-0.5*mu0*V(1)*Ms^2*sum(squeeze(m(1:end-1,:,1)).*squeeze(hc(1:end-1,:,1)),2)/abs(q);
Uc2=-0.5*mu0*V(2)*Ms^2*sum(squeeze(m(1:end-1,:,2)).*squeeze(hc(1:end-1,:,2)),2)/abs(q);
plot(t(1:end-1),Uc1);
hold on
plot(t(1:end-1),Uc2);
xlabel('Tempo (ns)')
ylabel('U_c (eV)')
legend('U_{c_1}','U_{c_2}');
sdf('P1');
%%
%subplot(1,3,3)
Ur=Ud1+Ud2+0.5*Uc1+0.5*Uc2;
plot(t(1:end-1),Ur);
xlabel('Tempo (ns)')
ylabel('U_r (eV)')
sdf('P1');