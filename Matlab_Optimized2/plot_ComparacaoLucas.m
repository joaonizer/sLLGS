figure('Position',[0,0,1900/2,1080/3]);
subplot(1,2,1);
plot(t(1:2:end),squeeze(m(1:2:end,:,1)));
hold on
set(gca,'colororderindex',1)
plot(test0.t,[test0.x0 test0.y0 test0.z0],'--', 'MarkerIndices',1:100:20000);
title('M0')
xlabel('Tempo (ns)')
ylabel('m')
%lh=legend('m_x J', 'm_y J','m_z J','m_x L','m_y L','m_z L')
%set(lh,'Location','northeastoutside')
xlim([0 20]);
ylim([-1.2 1.2]);
subplot(1,2,2);
plot(t(1:2:end),squeeze(m(1:2:end,:,2)));
hold on
set(gca,'colororderindex',1)
plot(test0.t,[test0.x1 test0.y1 test0.z1],'--', 'MarkerIndices',1:100:20000);
title('M1')
xlabel('Tempo (ns)')
ylabel('m')
%legend('m_x J', 'm_y J','m_z J','m_x L','m_y L','m_z L')
xlim([0 20]);
ylim([-1.2 1.2]);
sdf('P2')
%print('-dpng','-r300',[name '_3.png'])
print('-djpeg',[name '_3.jpeg'])
close all

