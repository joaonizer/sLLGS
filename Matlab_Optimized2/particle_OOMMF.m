function [Wi,Li,im_id]=particle_OOMMF(px,py,th,res_f,dx,dy)
% Define a imagem da partícula em formato .gif
% OBRIGATORIO:
% px=[-30 30 30 -30]; %nm Pontos para Largura X
% py=[30 45 -45 -45]; %nm Pontos para Comprimento Y
% OPCIONAIS:
% res_f = 4; %nm Fator de resolucao da Imagem
% dx = 10; %nm Desloca a particula em X
% dy = 10; %nm Desloca a particula em Y

if nargin<5
    res_f=4; %fator de resolucao 1 pixel/nm
    dx=10; % Desloca a particula em dx pixeis (nm)
    dy=10; % Desloca a particula em dy pixeis (nm)
end
%% Calcula Largura e Comprimento da Imagem
Wp=(max(px)-min(px)); % Largura da Particula
Lp=(max(py)-min(py)); % Comprimento da Particula
Wi=Wp+2*dx; % Largura da imagem
Li=Lp+2*dy; % Comprimento da imagem
if Wi<Wp+dx || Li<Lp+dy
    warning('Particle is out of image boundaries!')
end
%% Desloca para (0,0) a Particula
px=px+abs(min(px));
py=py+abs(min(py));
%% Altera a Resolucao da Imagem
px=px*res_f;
py=py*res_f;
m=Wi*res_f;
n=Li*res_f;
dx=dx*res_f;
dy=dy*res_f;
%% Gera a Imagem
Iu=zeros(m,n);
Id=zeros(m,n);
s=min(px):1:max(px); % Vetor de passos na Largura da Particula
[yup, ydown]=f_reta(s,px,py); % Calcula os valores de y sup e inf;
% Computa a imagem UP
for i=1+dx:m;
    for j=1+dy:n;
        if i>=min(px)+dx && i<=max(px)+dx && j<=yup(i-dx)+dy
            Iu(i,j)=1;
        else
            Iu(i,j)=0;
        end
    end
end
% Computa a imagem DOWN
for i=1+dx:m
    for j=1+dy:n
        if i>=min(px)+dx && i<=max(px)+dx && j<=ydown(i-dx)+dy
            Id(i,j)=0;
        else
            Id(i,j)=1;
        end
    end
end
%% Imagem final e Escrita em .gif
I=~flipud((Iu.*Id)'); % Faz o produto entre as imagens e rotaciona
i=0;
im_id=['particle' num2str(Wp) 'x' num2str(Lp) 'x' num2str(th)];
while exist(['./OOMMF_sim/' im_id '.gif'], 'file')
    i=i+1;
    im_id=['particle' num2str(Wp) 'x' num2str(Lp) 'x' num2str(th) '_' num2str(i)];
end
imwrite(I,['./OOMMF_sim/' im_id '.gif']); % salva a imagem