function [i_s]=compute_Ispin(N,s)
% Computa a corrente de spin aplicada
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

i_s=zeros(N,3);

for i=1:size(s,1)
    % Define os passos em cada direcao
    isx=(s(i,4)-s(i,1))/s(i,7);% passo em x
    isy=(s(i,5)-s(i,2))/s(i,7);% passo em y
    isz=(s(i,6)-s(i,3))/s(i,7);% passo em z
    
    % Testa se há variação entre inicial e final
    %X
    if isx==0
        dcx=1;
    else
        dcx=0;
    end
    %Y
    if isy==0
        dcy=1;
    else
        dcy=0;
    end
    %Z
    if isz==0
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
    
    % Cria a corrente aplicada baseado em s(i,:)
    %X
    if dcx
        i_s(in:in+s(i,7)-1,1)=ones(s(i,7),1)*s(i,1); %DC
    else
        i_s(in:in+s(i,7)-1,1)=s(i,1):isx:s(i,4)-isx; %rampa
    end
    %Y
    if dcy
        i_s(in:in+s(i,7)-1,2)=ones(s(i,7),1)*s(i,2); %DC
    else
        i_s(in:in+s(i,7)-1,2)=s(i,2):isy:s(i,5)-isy; %rampa
    end
    %Z
    if dcz
        i_s(in:in+s(i,7)-1,3)=ones(s(i,7),1)*s(i,3); %DC
    else
        i_s(in:in+s(i,7)-1,3)=s(i,3):isz:s(i,6)-isz; %rampa
    end
    
end
i_s(N+1,:)=i_s(N,:);% Replica a ultima posicao
end