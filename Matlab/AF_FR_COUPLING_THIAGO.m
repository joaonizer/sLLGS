%% Algitmo para Calcular os Tensores pro Thiago
%
%
%
clear all
clc
colors = [
    0.4940    0.1840    0.5560  % Roxo
    0.3010    0.7450    0.9330  % Azul
    0.6350    0.0780    0.1840  % Vermelho
    0.4660    0.6740    0.1880  % Verde
    0.4940    0.1840    0.5560  % Roxo
    0.3010    0.7450    0.9330  % Azul
    0.6350    0.0780    0.1840  % Vermelho
    0.4660    0.6740    0.1880  % Verde
    0.4940    0.1840    0.5560  % Roxo
    0.3010    0.7450    0.9330  % Azul
    0.6350    0.0780    0.1840  % Vermelho
    0.4660    0.6740    0.1880  % Verde
    ];
%% Dimensoes da Particula
part_n = 2;
d_min_x=10;
d_min_y=29;
w0=50;
l0=100;
th0=10;
platform='lin';

w=ones(1,part_n)*w0;  % width of particles

l=ones(1,part_n)*l0; % length of particles

th=ones(1,part_n)*th0;   %thickness of particles

d=[0,ones(1,part_n-1)]*d_min_x; % distance among particles

px=zeros(part_n,4);
py=px;
d_or=zeros(part_n,3);
cortes_y=zeros(part_n,4);

% Transforms w,l and d into px py and d_or
for i=1:length(w)
    [px(i,:),py(i,:),d_or]=write_Pontos(w,l,d,d_or,cortes_y(i,:),i);
end

dist = 1:1:50;
data = zeros(2,length(dist),4);
for i=1:2
    for j=1:length(dist)
        clc;
        fprintf('%d/%d\t\t%d/%d',i,2,j,length(dist));
        cor(1:3,:)=colors(1:3,:);
        if i==1
            d_or = [0 0 0
                w0+dist(j) 0 0
                ];
        else
            d_or = [0 0 0
                0 l0+dist(j) 0
                ];
        end
        
        
        
        
        %% Compute Tensores
        compute_NCND=1; % If TRUE computes the tensors again
        compute_PAR=1; % If TRUE uses paralel parfor to compute coupling tensor
        
        if compute_NCND
            [Nc,Nd,V]=compute_Tensores(px,py,d_or,th,part_n,platform,compute_PAR);
        else
            warning('Tensores não foram recalculados!');
        end
        
        data(i,j,:) = [i dist(j) Nc(1,1,1,2) Nc(2,2,1,2)];
    end
end
%% Fin X_dist and Y-dist for min
y_values=squeeze(data(1,:,4));
x_values = abs(squeeze(data(2,:,4)));
x=dist;

y_xx = pchip(x_values,x,y_values);
for i=1:length(dist)
fprintf('%3.2f %3.2f\n',dist(i), y_xx(i));
end
% Plota Resultados
figure('Position',[100 400 1600 400]);
subplot(1,3,1)
plot(squeeze(data(1,:,2)),squeeze(data(1,:,4)))
hold on
plot(squeeze(data(1,:,2)),abs(squeeze(data(2,:,4))))
xlabel('Distance (nm)');
lh=legend('Nc_{xx} Variation', '|Nc_{yy}| Variation');
set(lh,'Location','Best','Interpreter','tex')
subplot(1,3,2)
plot(squeeze(data(1,:,2)),abs(squeeze(data(2,:,4)))./squeeze(data(1,:,4)))
xlabel('Distance (nm)');
ylabel('|Nc_{yy}| / Nc_{xx}')
subplot(1,3,3)
plot(y_xx./dist);
xlabel('Distance (nm)')
ylabel('dY / dX')
sdf('P1');
%% Plota as particulas
% figure('Position',[800 100 300 600], ...
%     'Name','Alocação das Partículas');
% 
% espaco1=1:2:2*part_n-1;
% espaco2=0:1:part_n-1;
% p_a=espaco1*25+espaco2*d_min_x;
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