%%
tt0=0;
tt1=60;
subplot(3,1,1)
plot(t(t<tt1&t>tt0),squeeze(m(t<tt1&t>tt0,2,1:3)));
hold on
plot(t(t<tt1&t>tt0),squeeze(i_s(t<tt1&t>tt0,1,1:3))/0.25);
legend('P1','P2','P3'); grid on
xlim([tt0 tt1]);ylim([-1.2 1.2]);
subplot(3,1,2)
plot(t(t<tt1&t>tt0),squeeze(m(t<tt1&t>tt0,2,4:6)));
legend('P4','P5','P6'); grid on
xlim([tt0 tt1]);ylim([-1.2 1.2]);
subplot(3,1,3)
plot(t(t<tt1&t>tt0),squeeze(m(t<tt1&t>tt0,2,7:9)));
legend('P7','P8','P9'); grid on
xlim([tt0 tt1]);ylim([-1.2 1.2]);
sdf('P1')
