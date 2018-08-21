function [Ns, V]=write_FileDemag3D(px, py, th, platform)
% Calcula o tensor Ns de anisotropia de forma
% Retorna o tensor Ns ja normalizado e V em m3
% px - coordendas x dos pontos nm
% py - coordenada y dos pontos nm
% thick - espessura da particula nm
% platform - 'lin'ou 'win'
%Ex:
%px=[-25 25 25 -25];
%py=[0 50 -50 0];
%thick=[10];
%system('chmod -rwx demag3D2.o');
%% Verifica a restricoes em px
if (px(1)~=px(4) || px(2)~=px(3)) % Verifica a restri��o ditada pelo .90f
    error('Verifique as cordenadas px da geometria!\n-------> px(1) e px(4) DEVEM ser iguais!\n-------> px(2) e px(3) DEVEM ser iguais! %s\n',char(' '));
end
%% Gera o arquivo de entrada
if (length(px)~=length(py)) % verifica se as matrizes possuem tamanhos =
    error('Dimensões Inconsistentes!');
end

fid = fopen('IN_demag3D3.dat','w');
for i=1:length(px) % Imprime as coordenadas
    fprintf(fid,'%f\t%f\n',px(i),py(i));
end
fprintf(fid,'%f',th); % Imprime a espessura
fclose(fid); % fecha o arquivo
%% Calcula os parametros pelo .exe
if strcmp(platform,'lin')
    system('./demag3D3'); % roda o executavel para gerar o OUT_demag3D2
elseif strcmp(platform,'win')
    system('demag3D3.exe');
else
    error('Verifique a Plataforma');
end
%% Leitura e Normaliza��o do Tensor
% Calculo do Volume
V=abs(sum(trapz(px,py)))*th;
Ns=table2array(readtable('OUT_demag3D3.dat', 'Format', '%f %f %f'))/4/pi/V;
%Ns=table2array(readtable('OUT_demag3D3.dat', 'Format', '%f %f %f'));
V=V*1e-27; % retorna em m3
%% Plota a Geometria da Particula
% figure('Visible','off',...
%       'PaperPosition',[0 0 7.5 10],...
%       'PaperSize',[7.5 10])
%   plot(0,0,'ro','LineWidth',2);
%   hold on
%   plot([px px(1)]',[py py(1)]','LineWidth',2);
% set(gca,'BoxStyle','full','Box','on')
% title(['Geometria da Partícula - Thickness = ' num2str(th) ' nm']);
% xlabel('x (nm)');
% ylabel('y (nm)');
% xlim([min(px)-.1*max(px) max(px)+.1*max(px)]);
% ylim([min(py)-.1*max(px) max(py)+.1*max(px)]);
% grid on
% set(gcf,'Position',[0 0 1920 1080])
% print(gcf,'-dpdf','particle.pdf')

end