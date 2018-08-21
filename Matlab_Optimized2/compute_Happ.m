function [h_app]=compute_Happ(N,s)
% Computa o campo aplicado
% seguindo a logica Hrange do OOMMF
% N : é o numero de passos de simulação
%          | Inicial  |   Final   |Duracao|
%          |x_i y   z | x_f y   z | steps |
%          |1   2   3 | 4   5   6 |  7    |
% s=    [   0   0   0   1   0   0   221
%           1   0   0   1   1   0   221
%           1   1   0   0   1   0   221
%           0   1   0   0   0   0   221
%           0   0   0   0   0   0   221];
if sum(s(:,7))~=N
    error('N ~= Total shape steps');
end

global Ms;
t2am=7.943e5; % converte T para A/m
h_app=zeros(N,3);

for i=1:size(s,1)
    % Define os passos em cada direcao
    hx=(s(i,4)-s(i,1))/s(i,7);% passo em x
    hy=(s(i,5)-s(i,2))/s(i,7);% passo em y
    hz=(s(i,6)-s(i,3))/s(i,7);% passo em z

    % Testa se há variação entre inicial e final
    %X
    if hx==0
        dcx=1;
    else
        dcx=0;
    end
    %Y
    if hy==0
        dcy=1;
    else
        dcy=0;
    end
    %Z
    if hz==0
        dcz=1;
    else
        dcz=0;
    end
    
    % Determina o numero de posicoes ja utilizadas em cada direcao
    if i==1
        in=1; % indice inicial = 1
    else
        in=sum(s(1:i-1,7)); % acumula os indices
    end

    % Cria o campo aplicado baseado em s(i,:)
    %X
    if dcx
        h_app(in:in+s(i,7)-1,1)=ones(s(i,7),1)*s(i,1); %DC
    else
        h_app(in:in+s(i,7)-1,1)=s(i,1):hx:s(i,4)-hx; %rampa
    end
    %Y
    if dcy
        h_app(in:in+s(i,7)-1,2)=ones(s(i,7),1)*s(i,2); %DC
    else
        h_app(in:in+s(i,7)-1,2)=s(i,2):hy:s(i,5)-hy; %rampa
    end
    %Z
    if dcz
        h_app(in:in+s(i,7)-1,3)=ones(s(i,7),1)*s(i,3); %DC
    else
        h_app(in:in+s(i,7)-1,3)=s(i,3):hz:s(i,6)-hz; %rampa
    end
            
end
h_app(N+1,:)=h_app(N,:);% Replica a ultima posicao
h_app=h_app*t2am/Ms; %converte T --> A/m e normaliza por Ms (A/m)
end