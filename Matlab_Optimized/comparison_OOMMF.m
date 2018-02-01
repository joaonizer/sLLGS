r0K=read_ODT('/home/joaonizer/Desktop/Comparison/results300K.odt');


figure;
plot(t,squeeze(m(:,:,1)),'-','MarkerIndices',1:200:length(t),'MarkerFaceColor','w')
%plot(r0K1(:,1)*1e9,r0K1(:,2:4),'-','MarkerIndices',1:50:length(r0K(:,1)),'MarkerFaceColor','w')
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(r0K(:,1)*1e9,r0K(:,2:4),'o-','MarkerIndices',1:200:length(r0K(:,1)),'MarkerFaceColor','w')
legend('m_x sRK', 'm_y sRK', 'm_z sRK', 'm_x OOMMF', 'm_y OOMMF', 'm_z OOMMF')
xlim([0 1])
ylim([-1.2 1.2])
title('T=300K | \alpha=1 | Ms = 8e5 | OOMMF Cell Size 5x5x5')
xlabel('Tempo (ns)')
ylabel('M/M_s')
sdf('P1')

figure;
plot(t,1-sqrt(sum(squeeze(m(:,:,1)).^2,2)),'-','MarkerIndices',1:200:length(t),'MarkerFaceColor','w')
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(r0K(:,1)*1e9,1-sqrt(sum(r0K(:,2:4).^2,2)),'o-','MarkerIndices',1:200:length(t),'MarkerFaceColor','w')
xlim([0 1])
legend('sRK','OOMMF')
xlabel('Tempo (ns)')
ylabel('1-|M/M_s|')
sdf('P1')