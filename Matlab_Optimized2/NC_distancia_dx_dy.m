figure;
i=1;
j=2;
t=[0:N]*time_step;
plot(t/1e-9,squeeze(m1(:,:,1,i,j)));
hold on
set(gca,'ColorOrderIndex',1)
plot(t/1e-9,squeeze(m1(:,:,2,i,j)),'--');
xlabel('Tempo (ns)')
ylabel('m')
legend('m_{x_1}', 'm_{y_1}', 'm_{z_1}','m_{x_2}', 'm_{y_2}', 'm_{z_2}')
title(['dy=' num2str(i+(i-1)*2) ' nm'])
ylim([-1.1 1.1]);
sdf('P1')
print('-depsc',['dy=' num2str(i+(i-1)*2) num2str(j) 'nm.eps'])
%%
i=1;
figure;
t=[0:N]*time_step;
subplot(2,2,1)
plot(t/1e-9,squeeze(m1(:,:,1,i,1)));
xlabel('Tempo (ns)')
ylabel('m')
legend('m_x', 'm_y', 'm_z')
ylim([-1.1 1.1]);
subplot(2,2,2)
plot(t/1e-9,squeeze(m1(:,:,2,i,1)));
xlabel('Tempo (ns)')
ylabel('m')
legend('m_x', 'm_y', 'm_z')
ylim([-1.1 1.1]);
subplot(2,2,3)
plot(t/1e-9,squeeze(m1(:,:,1,i,2)));
xlabel('Tempo (ns)')
ylabel('m')
ylim([-1.1 1.1]);
legend('m_x', 'm_y', 'm_z')
subplot(2,2,4)
plot(t/1e-9,squeeze(m1(:,:,2,i,2)));
xlabel('Tempo (ns)')
ylabel('m')
ylim([-1.1 1.1]);
legend('m_x', 'm_y', 'm_z')
sdf('P1')
%%
erro=squeeze(sum(sum(sqrt(((m1(:,:,:,:,1)-m1(:,:,:,:,2))/(N+1)).^2))));

fid=fopen('resultado_NC_distancia_dy.txt','w');
fprintf(fid,'dx & \\multicolumn{3}{c}{erro} & \\multicolumn{3}{c}{Nc} \\\\\n');
fprintf(fid,' & erro1 & erro2 & Nc_{xx} & Nc_{yy} & Nc_{zz} \\\\\n');
for i =1:50
fprintf(fid,'%d %.2f & %.2f & %.6f & %.6f & %.6f \\\\\n',i+(i-1)*1,erro(1,i),erro(2,i),Nc1(1,4,i,1),Nc1(2,5,i,1),Nc1(3,6,i,1));
end
fclose(fid)

