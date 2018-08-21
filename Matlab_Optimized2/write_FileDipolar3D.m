function [I, V]=write_FileDipolar3D(px, py, th, d_or, platform)
% Calcula o tensor Ns de anisotropia de forma
% Retorna o tensor Ns ja normalizado e V em m3
% px - coordendas x dos pontos nm
% py - coordenada y dos pontos nm
% th - espessura da particula nm
% d_or - distancia da origem de cada particula [x y z]
% platform - 'lin'ou 'win'
%Ex:
% px=[-30 30 30 -30;-30 30 30 -30];
% py=[45 45 -45 -45;45 45 -45 -45];
% th=[15;15];
% d_or=[0 0 0;65 0 0];
%system('chmod -rwx demag3D2.o');
%% Verifica a restricoes em px
if (px(1,1)~=px(1,4) | px(1,2)~=px(:,3) | px(2,1)~=px(2,4) | px(2,2)~=px(2,3)) % Verifica a restricao ditada pelo .90f
    error('Verifique as cordenadas px da geometria!\n-------> px(1) e px(4) DEVEM ser iguais!\n-------> px(2) e px(3) DEVEM ser iguais! %s\n',char(' '));
end
%% Gera o arquivo de entrada
if (size(px,1)~=size(py,1) | size(px,2)~=size(py,2)) % verifica se as matrizes possuem tamanhos =
    error('Dimensões Inconsistentes!');
end
% Verifica se existe interseção entre as particulas
% OBS mesmo que a particula seja cortada sera considerado como um retangulo
w=max(px(1,2)-px(1,1),px(1,4)-px(1,3));
h=max(py(1,1)-py(1,4),py(1,2)-py(1,3));
A=[d_or(1,1)-w/2,d_or(1,2)-h/2,w,h];

w=max(px(2,2)-px(2,1),px(2,4)-px(2,3));
h=max(py(2,1)-py(2,4),py(2,2)-py(2,3));
B=[d_or(2,1)-w/2,d_or(2,2)-h/2,w,h];
if rectint(A,B) % rectint calcula a intersecao entre os retangulo A e B
    error('Existe Interseçao entre as partículas!');
end


fid = fopen('IN_dipolar3D.dat','w');
for i=1:2
    for j=1:4 % Imprime as coordenadas
        fprintf(fid,'%f\t%f\n',px(i,j),py(i,j));
    end
    fprintf(fid,'%f\n',th(1)); % Imprime a espessura
    fprintf(fid,'%f %f %f\n',d_or(i,:)); % Imprime a espessura
end
fclose(fid); % fecha o arquivo
%% Calcula os parametros pelo .exe
if strcmp(platform,'lin')
    system('./dipolar3D'); % roda o executavel para gerar o OUT_demag3D2
elseif strcmp(platform,'win')
    system('dipolar3D.exe');
else
    error('Verifique a Plataforma');
end
%% Leitura e Normalização do Tensor
% Calculo do Volume
for i=1:2
    V(i)=abs(sum(trapz(px(i,:),py(i,:))))*th(i);
end
fileID = fopen('OUT_dipolar3D.dat','r');
delimiter = ' ';
formatSpec = '%f%f%f%[^\n\r]';
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,0.0, 'ReturnOnError', false);
fclose(fileID);
I = [dataArray{1:end-1}];
V=V*1e-27; % retorna em m3
end