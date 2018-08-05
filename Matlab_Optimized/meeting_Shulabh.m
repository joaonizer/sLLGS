%%
plot(t,m(:,:,2))
hold on;
plot(t,i_s(:,1,2)./max(max(max(abs(i_s)))))
legend('m_x', 'm_y', 'm_z', 'is_x')
xlabel('Time (ns)')
ylabel('M/M_s, I_s/|I_s|')
sdf('P1')