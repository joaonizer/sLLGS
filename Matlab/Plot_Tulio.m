%% Inicio
clc
close all
%clear all

%% Configuracoes do Sistema

dx=20; % distancia x entre as particulas
dy=40; % distancia y entre as particulas

w=60;  % width of particles
l=120; % length of particles

r=40; % radius in nnm



%%

% Especificação do Grid
% 0 - sem partícula
% 1 - retangular
% 2 - circular
grid = [
    0 1 0 1 0
    1 1 0 1 1
    1 1 0 0 1
    ];

% Especificação do CLOCK
% 0 - INPUT
% 1 - Rede 1
% ...
% 4 - Rede 4
% 5 - Cruzador campo em y
clock = [
    0 1 2 3 4
    4 3 2 1 0
    4 3 2 1 0
    ];
% Deslocamento em nm de cada particula
xshift = [
    0 0 0 0 0
    0 10 0 0 0
    0 0 0 0 0
    ];
yshift = [
    0 -10 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    ];

%% Data Verification
data.name={'tipo';'xshift';'yshift';'clock'};
data.data={grid;xshift;yshift;clock};
disp(data.data{1,:});

%% Criação dos pontos

[px,py,d_or]=write_Pontos_New(w,l,r,dx,dy,grid,xshift,yshift);


%% Cores diferentes para particulas que estao em zonas de clock diferentes
cor=[204,85,90 % Input
    192,94,175 % Clock 1 
    109,164,81 % Clock 2
    118,117,203 % Clock 3
    192,137,61  % Clock 4
    74,180,195 % Cruzador
    ]/255;

%% Plota as particulas
figure('Position',[150 150 600 600], ...
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
                fill(px(index,:)+d_or(index,1)+px(1,2),py(index,:)+py(1,1)+d_or(index,2),cor(clock(i,j)+1,:),'EdgeColor','none','LineStyle','none');
                hold on
            case 2
                t=[0:0.1:360]*pi/180;
                x=r*cos(t)+d_or(index,1);
                y=r*sin(t)+d_or(index,2);
                fill(x,y,cor(clock(i,j)+1,:),'EdgeColor','none','LineStyle','none');
            otherwise
                fill(px(index,:)+d_or(index,1)+px(1,2),py(index,:)+py(1,1)+d_or(index,2),cor(clock(i,j)+1,:),'EdgeColor','none','LineStyle','none');
                hold on
                
        end
    end
end


daspect([1 1 1])
%title(['Distâncias ' num2str(d_min) 'nm'])
ylim([0 size(grid,1)*(l+dy)-dy])
xlim([0 size(grid,2)*(w+dx)-dx])
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