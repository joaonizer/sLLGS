plot([0:N/8]*time_step*1e9,squeeze(hT(1:N/8+1,1,1)))
hold on;
plot([0:N/8]*time_step*1e9,squeeze(hT1(1:N/8+1,1,1)))
legend('T=300 K', 'T=72 K')
xlabel('Tempo (ns)')
ylabel('\nu dW_x')
%%
plot(squeeze(h_app(:,1,1)))
set(gca,'xtick',[0 2500 5000 7500 10000])
xlim([0 10000])
xlabel('Passo')
ylabel('Clock')
sdf('P1')