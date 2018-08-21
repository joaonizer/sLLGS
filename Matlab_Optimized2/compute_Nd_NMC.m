function [Nd, V]=compute_Nd_NMC(px,py,th,part_n,compute_PAR,platform,NMC)
%% Create the Lists
data=[px py th'];
list=unique(data,'rows');
[~, particle_type]=ismember(data,list,'rows');
%% Compute Demag Tensors for the List Members
tic;
V=zeros(1,part_n);
fprintf('\nCalculando Tensores de Desmagnetização:\n')
temp_Nd=zeros(3,3,size(list,1));
if ~compute_PAR
    fprintf('Serial Mode: \n')
    for i=1:size(list,1)
        fprintf('%d ',i)
        
        if (min(abs(list(i,2:4))/abs(list(i,1))==1) && min(abs(list(i,6:8))/abs(list(i,5))==1))
            fprintf('rectangular')
            [temp_Nd(:,:,i), temp_V(i)]=write_FileDemag3D_rec_NMC(list(i,1:4), list(i,5:8), list(i,9), platform,NMC);
        else
            [temp_Nd(:,:,i), temp_V(i)]=write_FileDemag3D_noDAT_NMC(list(i,1:4), list(i,5:8), list(i,9), platform,NMC);
        end
    end
    fprintf('OK!\n')
else
    fprintf('Parallel Mode: ')
    parfor i=1:size(list,1)
        fprintf('%d ',i)
        if (min(abs(list(i,2:4))/abs(list(i,1))==1) && min(abs(list(i,6:8))/abs(list(i,5))==1))
            fprintf('retangular ')
            [temp_Nd(:,:,i), temp_V(i)]=write_FileDemag3D_rec_NMC(list(i,1:4), list(i,5:8), list(i,9), platform,NMC);
        else
            [temp_Nd(:,:,i), temp_V(i)]=write_FileDemag3D_noDAT_NMC(list(i,1:4), list(i,5:8), list(i,9), platform,NMC);
        end
    end
    fprintf(' OK!\n')
end
toc;
%% Replicates Computed Nd's
Nd=zeros(3*part_n,3*part_n);
for i=1:part_n
    index=3*(i-1)+1:3*(i-1)+3;
    Nd(index,index)=temp_Nd(:,:,particle_type(i));
    V(i)=temp_V(particle_type(i));
end
end