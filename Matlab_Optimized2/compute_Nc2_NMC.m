function [Nc]=compute_Nc2_NMC(px,py,th,d_or,part_n,compute_PAR,platform, NMC)
tic;
Nc=sparse(3*part_n,3*part_n);

fprintf('\nCalculando Tensores de Acoplamento:\n')
for i=1:part_n
    fprintf('%d -> ',i)
    index_i=3*(i-1)+1:3*(i-1)+3;
    for j=1:part_n
        index_j=3*(j-1)+1:3*(j-1)+3;
        if j~=i
            fprintf('%d ',j)
            [Nc(index_i,index_j), V2]=write_FileDipolar3D_noDAT_NMC([px(i,:);px(j,:)], ...
                [py(i,:);py(j,:)], [th(i);th(j)],...
                [d_or(j,:);d_or(i,:)], platform, NMC);
        else
            Nc(index_i,index_j)=zeros(3,3);
        end
    end
    fprintf('\n')
end
toc;
end