function [Nd, V]=write_FileDemag3D_rec_NMC(px, py, th, platform,NMC)
% Calcula o tensor Nd de anisotropia de forma
% Retorna o tensor Nd ja normalizado e V em m3
% px - coordendas x dos pontos nm
% py - coordenada y dos pontos nm
% thick - espessura da particula nm
% platform - 'lin'ou 'win'
%Ex:
%px=[-25 25 25 -25];
%py=[0 50 -50 0];
%thick=[10];
% Utiliza a versão sem arquivos de texto do algoritmo DEMAG3D3 ->
%                                                   DEMAG3d3_noDAT
% OBS: Esse código em particlar calcula apenas para particulas retangulares
% sem corte
%% Verifica a restricoes em px
if (px(1)~=px(4) || px(2)~=px(3)) % Verifica a restri��o ditada pelo .90f
    error('Verifique as cordenadas px da geometria!\n-------> px(1) e px(4) DEVEM ser iguais!\n-------> px(2) e px(3) DEVEM ser iguais! %s\n',char(' '));
end
%% Gera o arquivo de entrada
if (length(px)~=length(py)) % verifica se as matrizes possuem tamanhos =
    error('Dimensões Inconsistentes!');
end

%% Calcula os parametros pelo .exe
data_p=[
    num2str(px(2)-px(1)) ' '...
    num2str(py(2)-py(3)) ' '...
    num2str(th) ' ' num2str(NMC)...
    ];
if strcmp(platform,'lin')
    [flag,Nstr]=system(['echo '' ' data_p ' '' | ./Cpp/demag_rec_NMC']); % roda o executavel para gerar o OUT_demag3D2
else
    [flag,Nstr]=system(['echo ' data_p ' | .\Cpp\demag_rec_NMC.exe']); % no single quotes for windows
end
%% Leitura e Normaliza��o do Tensor
% Calculo do Volume
V=abs(sum(trapz(px,py)))*th;
Nd_temp=str2num(Nstr)/4/pi/V;
Nd=diag(Nd_temp);
%Ns=table2array(readtable('OUT_demag3D3.dat', 'Format', '%f %f %f'));
V=V*1e-27; % retorna em m3
end