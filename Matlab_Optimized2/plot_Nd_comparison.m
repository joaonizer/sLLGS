colormap(lines(3))
plot(t,squeeze(m_statistic(:,:,1)),'-');
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,squeeze(m_statistic(:,:,2)),'--');
ax.ColorOrderIndex = 1;
plot(t,squeeze(m_statistic(:,:,3)),'-.');
ylim([-1.2 1.2])
xlim([0 1])
xlabel('Tempo (ns)')
ylabel('M/M_s')
ax.ColorOrderIndex = 1;
h(1) = plot(NaN,NaN,'LineWidth',2);
h(2) = plot(NaN,NaN,'LineWidth',2);
h(3) = plot(NaN,NaN,'LineWidth',2);
ax.ColorOrderIndex = 1;
h(4) = plot(NaN,NaN,'LineWidth',2);
ax.ColorOrderIndex = 1;
h(5) = plot(NaN,NaN,'--','LineWidth',2);
ax.ColorOrderIndex = 1;
h(6) = plot(NaN,NaN,'-.','LineWidth',2);
lh=legend(h,'$m_x$','$m_y$','$m_z$','$\overline{N_d}$','$\overline{N_d}+2\sigma$','$\overline{N_d}-2\sigma$')
sdf('P1')
set(lh,'interpreter','latex','fontsize',10)