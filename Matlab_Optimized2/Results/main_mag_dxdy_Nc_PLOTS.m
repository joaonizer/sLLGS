%%
t=[0:N]*time_step/1e-9;
figure;
subplot(1,2,1)
plot(t,squeeze(m1(:,:,1,2,1)));
hold on
set(gca,'ColorOrderIndex',1);
plot(t,squeeze(m1(:,1,1,2,2)),'--o','MarkerIndices',1:350:N);
plot(t,squeeze(m1(:,2,1,2,2)),'--o','MarkerIndices',1:400:N);
plot(t,squeeze(m1(:,3,1,2,2)),'--o','MarkerIndices',1:450:N);
title('P_1 | dy=11nm')
ylim([-1.2 1.2])
xlim([0 4]);
xlabel('Tempo (ns)')
ylabel('m')
legend('m_x','m_y','m_z','m_{x0}','m_{y0}','m_{z0}');
subplot(1,2,2)
plot(t,squeeze(m1(:,:,2,2,1)));
hold on
set(gca,'ColorOrderIndex',1);
plot(t,squeeze(m1(:,1,2,1,2)),'--o','MarkerIndices',1:350:N);
plot(t,squeeze(m1(:,2,2,1,2)),'--o','MarkerIndices',1:400:N);
plot(t,squeeze(m1(:,3,2,1,2)),'--o','MarkerIndices',1:450:N);
title('P_2 | dy=11nm')
ylim([-1.2 1.2])
xlim([0 4]);
xlabel('Tempo (ns)')
ylabel('m')
legend('m_x','m_y','m_z','m_{x0}','m_{y0}','m_{z0}');
sdf('P1')
%%
figure;
subplot(1,2,1)
plot(t,squeeze(m1(:,:,1,23,1)));
hold on
set(gca,'ColorOrderIndex',1);
plot(t,squeeze(m1(:,1,1,23,2)),'--o','MarkerIndices',1:350:N);
plot(t,squeeze(m1(:,2,1,23,2)),'--o','MarkerIndices',1:400:N);
plot(t,squeeze(m1(:,3,1,23,2)),'--o','MarkerIndices',1:450:N);
title('P_1 | dy=221nm')
ylim([-1.2 1.2])
xlim([0 4]);
xlabel('Tempo (ns)')
ylabel('m')
legend('m_x','m_y','m_z','m_{x0}','m_{y0}','m_{z0}');

subplot(1,2,2)
plot(t,squeeze(m1(:,:,2,23,1)));
hold on
set(gca,'ColorOrderIndex',1);
plot(t,squeeze(m1(:,1,2,23,2)),'--o','MarkerIndices',1:350:N);
plot(t,squeeze(m1(:,2,2,23,2)),'--o','MarkerIndices',1:400:N);
plot(t,squeeze(m1(:,3,2,23,2)),'--o','MarkerIndices',1:450:N);
title('P_2 | dy=221nm')
ylim([-1.2 1.2])
xlim([0 4]);
xlabel('Tempo (ns)')
ylabel('m')
legend('m_x','m_y','m_z','m_{x0}','m_{y0}','m_{z0}');
sdf('P1')