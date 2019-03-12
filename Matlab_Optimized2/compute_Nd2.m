function [Nd, V]=compute_Nd2(px,py,th,part_n,compute_PAR,platform)
%% Create the Lists
%data=[px py th'];
%list=unique(data,'rows');
%[~, particle_type]=ismember(data,list,'rows');
%% Compute Demag Tensors for the List Members
tic;
V=zeros(1,part_n);
fprintf('\nCalculando Tensores de Desmagnetização:\n')
temp_Nd=zeros(3,3,part_n);
if ~compute_PAR
    fprintf('Serial Mode: \n')
    for i=1:part_n
        fprintf('%d ',i)
        [temp_Nd(:,:,i), temp_V(i)]=write_FileDemag3D_noDAT(px(i,:), py(i,:), th(i), platform);
    end
    fprintf('OK!\n')
else
    fprintf('Parallel Mode: ')
    parfor i=1:part_n
        fprintf('%d ',i)
        [temp_Nd(:,:,i), temp_V(i)]=write_FileDemag3D_noDAT(px(i,:), py(i,:), th(i), platform);
    end
    fprintf(' OK!\n')
end
toc;
%% Replicates Computed Nd's
Nd=sparse(3*part_n,3*part_n);
for i=1:part_n
    index=3*(i-1)+1:3*(i-1)+3;
    Nd(index,index)=temp_Nd(:,:,i);
    V(i)=temp_V(i);
end
end