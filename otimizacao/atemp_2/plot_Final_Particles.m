function plot_Final_Particles(x)
part_n=11;

platform='lin';
%x=[50 100 10 0 0 0 0 ...
%     50 100 24 -10 0 0 0 ];

w = [50 50 50 50 50 50 50 50 50 50 50];
l = [100 100 100 100 100 100 100 100 100 100 100];
th = [ 10 10 10 10 10 10 10 10 10 10 10];

cortes_y = -1*[
    0 x(1) -x(1) 0
    x(2) 0 0 -x(2)
    0 0 0 10 % AND GATE
    0 x(5) -x(5) 0
    x(6) 0 0 -x(6)
    0 0 0 0
    0 x(5) -x(5) 0
    x(6) 0 0 -x(6)
    0 0 0 0
    0 x(3) -x(3) 0
    x(4) 0 0 -x(4)
    ];

dx=[10 10 10 10 10 10 10 10 10 10 10];
d_or=[0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0];

for i=1:part_n
    [px(i,:),py(i,:),d_or]=write_Pontos(w,l,dx,d_or,cortes_y(i,:),i);
end
d_or=[
    0 0 0
    60 0 0
    60 -124 0
    120 -124 0
    180 -124 0
    240 -124 0
    300 -124 0
    360 -124 0
    420 -124 0
    420 0 0
    480 0 0];
%%
d_min=10;
cor = [
        0.4940    0.1840    0.5560  % Roxo
        0.4940    0.1840    0.5560  % Roxo
        0.3010    0.7450    0.9330  % Azul
        0.6350    0.0780    0.1840  % Vermelho
        0.6350    0.0780    0.1840  % Vermelho
        0.6350    0.0780    0.1840  % Vermelho
        0.4660    0.6740    0.1880  % Verde
        0.4940    0.1840    0.5560  % Roxo
        0.4940    0.1840    0.5560  % Roxo
        0.6350    0.0780    0.1840  % Vermelho
        0.4660    0.6740    0.1880  % Verde
        ];

%figure('Position',[800 100 300 600], ...
%    'Name','Alocação das Partículas');

espaco1=1:2:2*part_n-1;
espaco2=0:1:part_n-1;
p_a=espaco1*25+espaco2*d_min;

for i=1:part_n
    fill(px(i,:)+px(1,2)+d_or(i,1),py(i,:)+py(1,2)+d_or(i,2),cor(i,:),'EdgeColor','none','LineStyle','none')
    hold on
    text(d_or(i,1)+sum(abs(px(i,1:2)))/2,d_or(i,2)+py(i,2),['$P_{' num2str(i) '}$'],'Interpreter','latex','fontsize',9,'FontWeight','bold','Color','black','HorizontalAlignment','Center');
end


daspect([1 1 1])
title(['XOR Architecture'],'Interpreter','latex')
ylim([min(d_or(:,2))-25 max(d_or(:,2))+125])
xlim([min(d_or(:,1))-25 max(d_or(:,1))+75])
xlabel('$nm$','Interpreter','latex')
ylabel('$nm$','Interpreter','latex')

