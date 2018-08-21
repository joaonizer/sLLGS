function [Nc,Nd,V]=compute_Tensores(px,py,d_or,th,part_n,platform,compute_PAR,Nc,Nd)
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

V=zeros(1,part_n);
fprintf('\nCalculando Tensores de Desmagnetização:\n')
if ~compute_PAR
    fprintf('Serial Mode: \n')
    for i=1:part_n
        fprintf('%d ',i)
        index=3*(i-1)+1:3*(i-1)+3;
        [temp1, V(i)]=write_FileDemag3D_noDAT(px(i,:), py(i,:), th(i), platform);
        Nd(index,index)=temp1';
    end
    fprintf('OK!\n')
else
    fprintf('Parallel Mode: ')
    parfor i=1:part_n
        [temp1(i,:,:), V(i)]=write_FileDemag3D_noDAT(px(i,:), py(i,:), th(i), platform);
    end
    for i=1:part_n
        index=3*(i-1)+1:3*(i-1)+3;
        Nd(index,index)=squeeze(temp1(i,:,:))';
    end
    fprintf(' OK!\n')
end
toc;

tic;
fprintf('\nCalculando Tensores de Acoplamento:\n')
if ~compute_PAR
    fprintf('Serial Mode: \n')
    for i=1:part_n-1
        fprintf('%d -> ',i)
        for (j=i:1:part_n)
            if j~=i
                fprintf('%d ',j)
                [Nc(:,:,i,j), V2]=write_FileDipolar3D_noDAT([px(i,:);px(j,:)], ...
                    [py(i,:);py(j,:)], [th(i);th(j)],...
                    [d_or(j,:);d_or(i,:)], platform);
                K = 4*pi*V2/1e-27;
                Nc(:,:,j,i)=Nc(:,:,i,j)*K(1)/K(2);
            end
        end
        fprintf('\n')
    end
else
    fprintf('Parallel Mode: \n')
    for i=1:part_n
        fprintf('%d -> ',i)
        parfor j=i+1:part_n
            [temp2(:,:,i,j), ~]=write_FileDipolar3D_noDAT([px(i,:);px(j,:)], ...
                [py(i,:);py(j,:)], [th(i);th(j)],...
                [d_or(j,:);d_or(i,:)], platform);
        end
        fprintf('OK!\n')
    end
    for i=1:part_n
        i_index=3*(i-1)+1:3*(i-1)+3;
        K(1) = 4*pi*V(i)/1e-27;
        for j=i+1:part_n
            K(2) = 4*pi*V(j)/1e-27;
            j_index=3*(j-1)+1:3*(j-1)+3;
            %temp2(temp2<0.00001)=0;
            Nc(i_index,j_index)=squeeze(temp2(:,:,i,j));
            Nc(j_index,i_index)=Nc(i_index,j_index)*K(1)/K(2);
        end
    end
    
end

toc;
end