function [Nc,Nd,V]=compute_Tensores(px,py,d_or,th,part_n,platform,compute_PAR)
%% Função que computa os tensores de Acoplamento e Desmagnetização
% ENTRADAS:
% px py d_or th - dados das partículas
% platform - 'win' ou 'lin' para seleção dos executáveis
%
% SAÍDAS:
% Nc - Tensor de acoplamento (3,3,part_n,part_n)
% Nd - Tensor de desmagnetização (3,3,part_n)
% V - volume das partículas
tic;

Nc=zeros(3,3,part_n,part_n);
fprintf('\nCalculando Tensores de Acoplamento:\n')
if ~compute_PAR
    for i=1:part_n-1
        fprintf('%d -> ',i)
        for j=1:part_n
            if j~=i
                fprintf('%d ',j)
                [Nc(:,:,i,j), V2]=write_FileDipolar3D_noDAT([px(i,:);px(j,:)], ...
                    [py(i,:);py(j,:)], [th(i);th(j)],...
                    [d_or(j,:);d_or(i,:)], platform);
                %K = 4*pi*V2/1e-27;
                %Nc(:,:,j,i)=Nc(:,:,i,j)*K(1)/K(2);
            end
        end
        fprintf('\n')
    end
else
    for i=1:part_n
        fprintf('%d -> ',i)
        parfor j=1:part_n
            if j~=i
                [Nc(:,:,i,j), ~]=write_FileDipolar3D_noDAT([px(i,:);px(j,:)], ...
                    [py(i,:);py(j,:)], [th(i);th(j)],...
                    [d_or(j,:);d_or(i,:)], platform);
            end
        end
        fprintf('OK!\n')
    end
end

toc;

tic;

Nd=zeros(3,3,part_n);
V=zeros(1,part_n);
fprintf('\nCalculando Tensores de Desmagnetização:\n')
if ~compute_PAR
    for i=1:part_n
        fprintf('%d ',i)
        [Nd(:,:,i), V(i)]=write_FileDemag3D_noDAT(px(i,:), py(i,:), th(i), platform);
    end
    fprintf('\n')
else
    parfor i=1:part_n
        fprintf('%d ',i)
        [Nd(:,:,i), V(i)]=write_FileDemag3D_noDAT(px(i,:), py(i,:), th(i), platform);
    end
end
toc;
end