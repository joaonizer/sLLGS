%%
figure;
for i=1:length(temperatures)
    for j=1:length(alphas)
        subplot(length(temperatures),length(alphas),(i-1)*3+j)
        plot(t,squeeze(m_alphas(i,j,:,:)))
        title(['\alpha=' num2str(alphas(j)) ' T=' num2str(temperatures(i)) 'K'])
        ylim([-1.2 1.2])
        xlim([0 max(t)])
        xlabel('Tempo (ns)')
        ylabel('m')
    end
end
legend('m_x','m_y','m_z')
sdf('P1')