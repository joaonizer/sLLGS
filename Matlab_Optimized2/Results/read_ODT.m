function [data]=read_ODT(f_id)
% Função para leitura dos dados exportados pelo OOMMF
% f_id - Arquivo no formato .ODT
% data - Resultados da leitura
%       Stage,Mx,My,Mz,Bx,By,Bz, Demag Energy, Total Energy
% maxM - Valor máximo para normalização da magnetização
%
% OBS: Pode-se alterar a função para exportar mais dados (see .ODT)
fid=fopen(f_id,'r');
x = importdata(f_id); %// read lines
data = zeros(numel(x)-6,25); %// preallocate with 26 cols (acccording to the file)
zz = sscanf(fgetl(fid),'%f'); %% Pula as cinco primeiras linhas
zz = sscanf(fgetl(fid),'%f');
zz = sscanf(fgetl(fid),'%f');
zz = sscanf(fgetl(fid),'%f');
zz = sscanf(fgetl(fid),'%f');
for n = 1:numel(x)-6
    %     zz = sscanf(fgetl(fid),'%f');
    %    size(zz)
    %     n
    %     zz
    data(n,:)=sscanf(fgetl(fid),'%f');
    
    %y(n,:) = regexp(x{n}, '(\s\s+)|\t', 'split'); %// split each line into
    %// columns using as separator either more than one space or a tab
    %//(according to your file)
end
%y=y(6:numel(x),2:25); %% Exclui Col 1 e 26 (sao vazias)
%data(:,:)=reshape(str2num(strjoin(y)),numel(x)-5,24);
% for i=1:numel(x)-5 % Exclui o cabeçalho (5 primeiras linhas)
%     for j=1:24
%         data(i,j)=str2num(strjoin(y(i,j)));%% Converte os dados de cell-->str-->num
%     end
% end
%data=y;
%maxM=max(max(abs(data(50:numel(x)-6,19:21))));
data=[
    data(:,25)... % STAGE
    data(:,21), data(:,22), data(:,23)... % M (Magnetization)
    data(:,15), data(:,16), data(:,17)... % B (Zeeman Field)
    data(:,12), data(:, 1)... % Demag e Total Energy
    data(:, 13)... % Zeeman Energy
    data(:, 10)... % Stage Max Spin Angle
    ];
end