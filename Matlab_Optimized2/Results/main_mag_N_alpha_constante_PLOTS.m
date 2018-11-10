%%
figure;
for i=1:length(alphas)
    subplot(2,ceil(length(alphas)/2),i)
plot(t,squeeze(m_alphas(i,:,:)))
title(['\alpha=' num2str(alphas(i))])
ylim([-1.2 1.2])
xlim([0 max(t)])
xlabel('Tempo (ns)')
ylabel('m')
end
legend('m_x','m_y','m_z')
sdf('P1')
%%
figure;
for i=1:length(alphas)
    subplot(2,ceil(length(alphas)/2),i)
plot([0:1:alphas(i)]*4/(alphas(i)),squeeze(m_alphas(i,1:alphas(i)+1,:)))
title(['N=' num2str(alphas(i))])
ylim([-1.2 1.2])
xlim([0 max(t)])
xlabel('Tempo (ns)')
ylabel('m')
end
legend('m_x','m_y','m_z')
sdf('P1')