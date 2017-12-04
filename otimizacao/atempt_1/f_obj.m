function [obj]=f_obj(x)
%% Constantes
% Ms=800e3; % A/m
% Ku=0.26*1; %eV
% q = 1.60217662e-19; % carga do eletron C
% mu0=4*pi*1e-7; % H/m ou T.m/A
% T=100; % Kelvin
% kb=1.38064852e-23; %m2.kg.s-2.K-1
% kbT=kb*T;% J
% A=13e-12; % Exchange Stiffness J/m
% t2am=7.943e5; % converte T para A/m
% gammamu0=2.211e5; % (sA/m)-1

part_n=4;

platform='win';
%x=[50 100 10 0 0 0 0 ...
%     50 100 24 -10 0 0 0 ];

w = [50 x(1) x(8) x(8)];
l = [100 x(2) x(9) x(9)];
th = [ 10 10 10 10];

dx=[0 x(3) 0 0];
dy=[0 0 x(10) x(10)];
cortes_y = [0     0     0   -10 % AND GATE
    x(4:7)
    x(11:14)
    -x(14:-1:11)
    ];


d_or=[0 0 0
    0 0 0
    0 0 0
    0 0 0];

for i=1:part_n
    [px(i,:),py(i,:),d_or]=write_Pontos(w,l,dx,d_or,cortes_y(i,:),i);
end
d_or(3,1) = 0;
d_or(3,2) = 100 + dy(3);
d_or(4,1) = 0;
d_or(4,2) = -100 - dy(3);
%% Computa os tensores
Nc = zeros(3,3,part_n,part_n);
for i=2:part_n-1
    %fprintf('%d -> ',i)
    for j=i:part_n
        if j~=i
            %fprintf('%d ',j)
            [Nc(:,:,i,j), V2]=write_FileDipolar3D_noDAT([px(i,:);px(j,:)], ...
                [py(i,:);py(j,:)], [th(i);th(j)],...
                [d_or(j,:);d_or(i,:)], platform);
            K = 4*pi*V2/1e-27;
            Nc(:,:,j,i)=Nc(:,:,i,j)*K(1)/K(2);
        end
    end
    %fprintf('\n')
end

obj = abs( Nc(2,2,3,2) ) + abs( Nc(2,2,2,3) ) + abs( Nc(2,2,4,2) ) + abs( Nc(2,2,2,4) );
%%
% d_min=10;
% cor = [
%         0.4940    0.1840    0.5560  % Roxo
%         0.3010    0.7450    0.9330  % Azul
%         0.6350    0.0780    0.1840  % Vermelho
%         0.4660    0.6740    0.1880  % Verde
%         ];
% 
% figure('Position',[800 100 300 600], ...
%     'Name','Alocação das Partículas');
% 
% espaco1=1:2:2*part_n-1;
% espaco2=0:1:part_n-1;
% p_a=espaco1*25+espaco2*d_min;
% 
% for i=1:part_n
%     fill(px(i,:)+px(1,2)+d_or(i,1),py(i,:)+py(1,2)+d_or(i,2),cor(i,:),'EdgeColor','none','LineStyle','none')
%     hold on
%     text(d_or(i,1)+sum(abs(px(i,1:2)))/2,d_or(i,2)+py(i,2),['$P_{' num2str(i) '}$'],'Interpreter','latex','fontsize',9,'FontWeight','bold','Color','black','HorizontalAlignment','Center');
% end
% 
% 
% daspect([1 1 1])
% title(['XOR Architecture'],'Interpreter','latex')
% ylim([min(d_or(:,2))-25 max(d_or(:,2))+125])
% xlim([min(d_or(:,1))-25 max(d_or(:,1))+75])
% xlabel('$nm$','Interpreter','latex')
% ylabel('$nm$','Interpreter','latex')
% sdf('P1');

end