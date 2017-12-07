clc
%% Distances d_or
grid=[
    1   1   1   0   0   0   0
    0   0   2   1   1   0   0
    1   1   1   0   1   0   0
    0   0   0   0   1   0   0
    0   0   0   0   3   1   1
    0   0   0   0   1   0   0
    1   1   1   0   1   0   0
    0   0   2   1   1   0   0
    1   1   1   0   0   0   0
    ];
n=size(grid,1);
m=size(grid,2);
count=1;
for i=1:m
    for j=1:n
        if grid(j,i)
            fprintf('dx(%d)\tdy(%d)\t0\t%% P%d\n',i,n-j+1,count);
            count=count+1;
        end
    end
end


%% Plot Place
for i=1:m
    for j=1:n
        if grid(j,i)
            fprintf('%d,',i+(j-1)*m);
            count=count+1;
        end
    end
    fprintf(' ... \n')
end

%% Cortes Y
fprintf('cortes_y = [\n');
count=1;
for i=1:m
    for j=1:n
        if grid(j,i)==1 % normal
            fprintf('0\t0\t0\t0\t%% P%d\n',count);
            count=count+1;
        elseif grid(j,i)==2 %and
            fprintf('0\t0\t0\t-10\t%% P%d\n',count);
            count=count+1;
        elseif grid(j,i)==3 %or
            fprintf('10\t0\t0\t0\t%% P%d\n',count);
            count=count+1;
        end
    end
end
fprintf('];\n');
%
% Nc(:,:,15,[9,11])=zeros(3,3,1,2);
% Nc(:,:,[9,11],15)=zeros(3,3,2,1);
% Nc(:,:,16,[12,14])=zeros(3,3,1,2);
% Nc(:,:,[12,14],16)=zeros(3,3,2,1);
% Nc(:,:,24,[19,21])=zeros(3,3,1,2);
% Nc(:,:,[19,21],24)=zeros(3,3,2,1);
% Nc(:,:,16,22)=zeros(3,3,1,1);
% Nc(:,:,22,16)=zeros(3,3,1,1);
% Nc(:,:,15,18)=zeros(3,3,1,1);
% Nc(:,:,18,15)=zeros(3,3,1,1);
