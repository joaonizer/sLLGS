%r0K=read_ODT('C:\Users\joaon\Desktop\Comparison\rectangular_50_100\rectangular_50_100nm.odt');
%r0K=read_ODT('C:\Users\joaon\Desktop\Comparison\slanted_up\slanted_u_20nm.odt');
r0K=read_ODT('C:\Users\joaon\Desktop\Comparison\slanted_up_down\slanted_ud_20nm.odt');
%figure('Position',[0 0 850 600]);
 figure;
 plot(NaN,NaN,'Linewidth',2);
 hold on
 plot(NaN,NaN,'Linewidth',2);
  plot(NaN,NaN,'Linewidth',2);
  ax = gca;
ax.ColorOrderIndex = 1;
plot(t(1:10:end),squeeze(m(1:10:end,:,1)),'--','MarkerIndices',1:200:length(t),'MarkerFaceColor','w')
%plot(r0K1(:,1)*1e9,r0K1(:,2:4),'-','MarkerIndices',1:50:length(r0K(:,1)),'MarkerFaceColor','w')
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(r0K(1:10:end,1)*1e9,r0K(1:10:end,2:4),'-','MarkerIndices',1:200:length(r0K(:,1)),'MarkerFaceColor','w')
%lh=legend('m_x sRK', 'm_y sRK', 'm_z sRK', 'm_x OOMMF', 'm_y OOMMF', 'm_z OOMMF');
lh=legend('m_x', 'm_y', 'm_z');

xlim([0 1])
ylim([-1.2 1.2])
%title('T=0K | \alpha=0.05 | Ms = 8e5 | OOMMF Cell Size 5x5x10')
xlabel('Tempo (ns)')
ylabel('m')
sdf('P1')
set(lh,'FontSize',18)
set(gca,'LooseInset',get(gca,'TightInset')+0.01);

%set(gcf,'color','none')

%figure('Position',[0 0 800 600]);
 figure;
plot(t,log10(abs(1-sqrt(sum(squeeze(m(:,:,1)).^2,2)))),'--','MarkerIndices',1:200:length(t),'MarkerFaceColor','w')
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(r0K(:,1)*1e9,log10(abs(1-sqrt(sum(r0K(:,2:4).^2,2)))),'-','MarkerIndices',1:200:length(t),'MarkerFaceColor','w')
xlim([0 1])
lh=legend('sRK','OOMMF');
xlabel('Tempo (ns)')
ylabel('log_{10}(1-|m|)')
sdf('P1')
set(lh,'FontSize',12)
set(gca,'LooseInset',get(gca,'TightInset')+0.01);
fclose all