%% Inicio
clc
close all
%clear all

%% Configuracoes do Sistema

part_n=6; % quantidade de particulas

dx=10; % distancia x entre as particulas
dy=20; % distancia y entre as particulas

platform='lin';

w=50;  % width of particles

l=100; % length of particles
r=w/2; % radius in nnm

px=zeros(part_n,4);
py=px;
d_or=zeros(part_n,3);

%% Especificações do Grid
grid = [
    1 2 1
    2 1 2
    ];
clock = [
    1 0 2
    0 3 4
    ];
xshift = [
    0 0 0
    0 0 0
    ];
yshift = [
    0 0 0
    0 0 0
    ];
%% Data Verification
data.name={'tipo';'xshift';'yshift';'clock'};
data.data={grid;xshift;yshift;clock};
disp(data.data{1,:});

%% Criação dos pontos

[px,py,d_or]=write_Pontos_New(w,l,r,dx,dy,grid,xshift,yshift);


%% Cores diferentes para particulas que estao em zonas de clock diferentes
cor=[204,85,90
    74,180,195
    192,94,175
    109,164,81
    118,117,203
    192,137,61]/255;

%% Plota as particulas
figure('Position',[0 0 600 600], ...
    'Visible','on', ...
    'Name','Alocação das Partículas');
d_min=dx;
espaco1=1:2:2*part_n-1;
espaco2=0:1:part_n-1;
p_a=espaco1*25+espaco2*d_min;

for i=1:size(grid,1) % iterate over row
    for j=1:size(grid,2) % iterate over columns
        index = (i-1)*size(grid,2) + j;
        switch grid(i,j)
            case 1
                fill(px(index,:)+d_or(index,1)+px(1,2),py(index,:)+py(1,1)+d_or(index,2),cor(index,:),'EdgeColor','none','LineStyle','none');
                hold on
            case 2
                t=[0:0.1:360]*pi/180;
                x=r*cos(t)+d_or(index,1);
                y=r*sin(t)+d_or(index,2);
                fill(x,y,cor(index,:),'EdgeColor','none','LineStyle','none');
            otherwise
                fill(px(index,:)+d_or(index,1)+px(1,2),py(index,:)+py(1,1)+d_or(index,2),cor(index,:),'EdgeColor','none','LineStyle','none');
                hold on
                
        end
    end
end


daspect([1 1 1])
%title(['Distâncias ' num2str(d_min) 'nm'])
ylim([0 max(d_or(:,2))+l/2])
xlim([0 max(d_or(:,1))+w/2])
view([0 -90])
%xlabel('nm')
%ylabel('nm')
box off
axis off
% Usado para Cortar os Espaços em Branco
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%

% Pega a imagem da do PLOT com FRAME
image=getframe;
[image_gif,map_gif] = rgb2ind(image.cdata,256);
% Plota a imagem como gif
imwrite(image_gif,map_gif,'teste.gif','gif');