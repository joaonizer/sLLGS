m_result=round(m_test);
mask=[1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1];
p=zeros(2,13);
for k=1:2
    for i=1:100
        erro =0;
        j=1;
        if i<51
            mult=1;
        else
            mult=-1;
        end
        while erro==0 && j<14
            if m_result(k,i,2,j)==mult*mask(j)
                j=j+1;
                erro=0;
            else
                erro=1;
                p(k,j:end)=p(k,j:end)+1;
            end
        end
    end
end

plot(p'/100,'-o');
ylim([0 1])
xlim([1 13])
xlabel('Índice da Partícula')
set(gca,'xtick',[1 2 3 4 5 6 7 8 9 10 11 12 13]);
set(gca,'xticklabel',{'In'; '1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12'})
ylabel('Percentual de Erro')
sdf('P1')
%%
figure('rend','painters','pos',[10 10 1200 300]);
plot(t(1:100:end),squeeze(i_s(1:100:end,1,2))/zeta(2)*J_shm/1e12,'--o','Color',colors(2,:),'MarkerIndices',[1:15:length(t)/100])
hold on
plot(t(1:100:end),squeeze(i_s(1:100:end,1,5))/zeta(5)*J_shm/1e12,'--^','Color',colors(3,:),'MarkerIndices',[1:20:length(t)/100])
plot(t(1:100:end),squeeze(i_s(1:100:end,1,8))/zeta(8)*J_shm/1e12,'--s','Color',colors(4,:),'MarkerIndices',[1:25:length(t)/100])
plot(t(1:100:end),squeeze(i_s(1:100:end,1,11))/zeta(11)*J_shm/1e12,'--x','Color',colors(5,:),'MarkerIndices',[1:30:length(t)/100])
xlim([0 35])
xlabel('Tempo (ns)')

legend('Clock_1','Clock_2','Clock_3','Clock_4')
sdf('P1')

%%
figure('rend','painters','pos',[10 10 1200 350]);% IN P1 P2 P3
plot(t(1:100:end),squeeze(m(1:100:end,2,1)),'-^','Color',colors(1,:),'MarkerIndices',[1:10:length(t)/100])
hold on;
plot(t(1:100:end),squeeze(m(1:100:end,2,2)),'-','Color',colors(2,:),'MarkerIndices',[1:10:length(t)/100])
plot(t(1:100:end),squeeze(m(1:100:end,2,3)),'--x','Color',colors(2,:),'MarkerIndices',[1:15:length(t)/100])
plot(t(1:100:end),squeeze(m(1:100:end,2,4)),'-.o','Color',colors(2,:),'MarkerIndices',[1:20:length(t)/100])
% title(['\downarrow_{In}\uparrow_{P_1}\downarrow_{P_2}\uparrow_{P_3}\color{black}',...
%         '\downarrow_{P_4}\uparrow_{P_5}\color{black}\downarrow_{P_6}',...
%         '\downarrow_{P_7}\uparrow_{P_8}\color{black}\downarrow_{P_9}',...
%     '\downarrow_{P_{10}}\uparrow_{P_{11}}\color{black}\downarrow_{P_{12}}']);
%     
xlabel('Tempo (ns)')
xlim([0 35])
ylabel('m_y')
ylim([-1.2 1.2])
yyaxis right
plot(t(1:100:end),squeeze(i_s(1:100:end,1,2))/zeta(2)*J_shm/1e12,'--','Color','k','MarkerIndices',[1:15:length(t)/100])
ylabel('J_{MHS_y}\muA/nm^2')
ylim([-25-25*.2 25+25*.2])
set(gca,'YColor',[0 0 0]);
legend('Input','P_1','P_2','P_3','Clock_1','P_4','P_5','P_6','P_7','P_8','P_9','P_{10}','P_{11}','P_{12}')
sdf('P1')
%
figure('rend','painters','pos',[10 10 1200 350]);% IN P1 P2 P3% P4 P5 P6
plot(t(1:100:end),squeeze(m(1:100:end,2,5)),'-','Color',colors(3,:),'MarkerIndices',[1:10:length(t)/100])
hold on;
plot(t(1:100:end),squeeze(m(1:100:end,2,6)),'--x','Color',colors(3,:),'MarkerIndices',[1:10:length(t)/100])
plot(t(1:100:end),squeeze(m(1:100:end,2,7)),'-.o','Color',colors(3,:),'MarkerIndices',[1:10:length(t)/100])
% title(['\downarrow_{In}\uparrow_{P_1}\downarrow_{P_2}\uparrow_{P_3}\color{black}',...
%         '\downarrow_{P_4}\uparrow_{P_5}\color{black}\downarrow_{P_6}',...
%         '\downarrow_{P_7}\uparrow_{P_8}\color{black}\downarrow_{P_9}',...
%     '\downarrow_{P_{10}}\uparrow_{P_{11}}\color{black}\downarrow_{P_{12}}']);
%     
xlabel('Tempo (ns)')
xlim([0 35])
ylabel('m_y')
ylim([-1.2 1.2])
yyaxis right
plot(t(1:100:end),squeeze(i_s(1:100:end,1,5))/zeta(5)*J_shm/1e12,'--','Color','k','MarkerIndices',[1:15:length(t)/100])
ylabel('J_{MHS_y}\muA/nm^2')
ylim([-25-25*.2 25+25*.2])
set(gca,'YColor',[0 0 0]);
legend('P_4','P_5','P_6','Clock_2','P_7','P_8','P_9','P_{10}','P_{11}','P_{12}')
sdf('P1')

figure('rend','painters','pos',[10 10 1200 350]);% IN P1 P2 P3;% P7 P8 P9
plot(t(1:100:end),squeeze(m(1:100:end,2,8)),'-','Color',colors(4,:),'MarkerIndices',[1:10:length(t)/100])
hold on;
plot(t(1:100:end),squeeze(m(1:100:end,2,9)),'--x','Color',colors(4,:),'MarkerIndices',[1:10:length(t)/100])
plot(t(1:100:end),squeeze(m(1:100:end,2,10)),'-.o','Color',colors(4,:),'MarkerIndices',[1:10:length(t)/100])
% title(['\downarrow_{In}\uparrow_{P_1}\downarrow_{P_2}\uparrow_{P_3}\color{black}',...
%         '\downarrow_{P_4}\uparrow_{P_5}\color{black}\downarrow_{P_6}',...
%         '\downarrow_{P_7}\uparrow_{P_8}\color{black}\downarrow_{P_9}',...
%     '\downarrow_{P_{10}}\uparrow_{P_{11}}\color{black}\downarrow_{P_{12}}']);
%     
xlabel('Tempo (ns)')
xlim([0 35])
ylabel('m_y')
ylim([-1.2 1.2])
yyaxis right
plot(t(1:100:end),squeeze(i_s(1:100:end,1,8))/zeta(8)*J_shm/1e12,'--','Color','k','MarkerIndices',[1:15:length(t)/100])
ylabel('J_{MHS_y}\muA/nm^2')
ylim([-25-25*.2 25+25*.2])
set(gca,'YColor',[0 0 0]);
legend('P_7','P_8','P_9','Clock_3','P_{10}','P_{11}','P_{12}')
sdf('P1')
figure('rend','painters','pos',[10 10 1200 350]);% IN P1 P2 P3; % P10 P11 P12
plot(t(1:100:end),squeeze(m(1:100:end,2,11)),'-','Color',colors(5,:),'MarkerIndices',[1:10:length(t)/100])
hold on;
plot(t(1:100:end),squeeze(m(1:100:end,2,12)),'--x','Color',colors(5,:),'MarkerIndices',[1:10:length(t)/100])
plot(t(1:100:end),squeeze(m(1:100:end,2,13)),'-.o','Color',colors(5,:),'MarkerIndices',[1:10:length(t)/100])
% title(['\downarrow_{In}\uparrow_{P_1}\downarrow_{P_2}\uparrow_{P_3}\color{black}',...
%         '\downarrow_{P_4}\uparrow_{P_5}\color{black}\downarrow_{P_6}',...
%         '\downarrow_{P_7}\uparrow_{P_8}\color{black}\downarrow_{P_9}',...
%     '\downarrow_{P_{10}}\uparrow_{P_{11}}\color{black}\downarrow_{P_{12}}']);
%     
xlabel('Tempo (ns)')
xlim([0 35])
ylabel('m_y')
ylim([-1.2 1.2])
yyaxis right
plot(t(1:100:end),squeeze(i_s(1:100:end,1,11))/zeta(11)*J_shm/1e12,'--','Color','k','MarkerIndices',[1:15:length(t)/100])
ylabel('J_{MHS_y}\muA/nm^2')
ylim([-25-25*.2 25+25*.2])
set(gca,'YColor',[0 0 0]);
legend('P_{10}','P_{11}','P_{12}','Clock_4')
sdf('P1')