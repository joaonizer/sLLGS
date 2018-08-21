function [Nc]=compute_Nc(px,py,th,d_or,radius,grid,part_n,compute_PAR,platform)
%% Cria o Grid com índice de Partícula TopBottom-LeftRight
count=1;
for j=1:size(grid,2) %% create a grid with particle ID
    for i=1:size(grid,1)
        if grid(i,j)>0
            grid_part(i,j)=count;
            count=count+1;
        end
    end
end

%% Criacao da Lista de Dados de Acoplamentos Únicos
count=1;
list2=zeros(0,24);
coupling_grid=zeros(part_n,part_n);
for j=1:size(grid,2)
    for i=1:size(grid,1)
        if grid(i,j)>0
            
            
            clear temp_list index
            index = grid_part(i,j); % indice da particula atual
            y2 = i-radius:i+radius;
            x2 = j-radius:j+radius;
            y1=y2(y2>0); % remove valores menores que 1
            x1=x2(x2>0);
            y=y1(y1<=size(grid,1)); % Remove valores superiores ao tamanho do grid
            x=x1(x1<=size(grid,2));
            temp_count=1;
            for nnn=1:length(y) % Converte os indices para número de particulas
                for mmm=1:length(x)
                    %disp((n-1)*length(y)+m)
                    if (grid_part(y(nnn),x(mmm))>0 && grid_part(y(nnn),x(mmm))~=index)
                        temp_list(temp_count) = grid_part(y(nnn),x(mmm));
                        temp_count=temp_count+1;
                    end
                end
            end
            
            % Cria Lista e o Grid de Acoplamento
            for nnn=1:length(temp_list)
                vec=[px(index,:) py(index,:) th(index) 0 0 0, ...
                    px(temp_list(nnn),:) py(temp_list(nnn),:) th(temp_list(nnn)), ...
                    d_or(temp_list(nnn),:)-d_or(index,:)];
                [Lia,Loc]=ismember(vec,list2,'rows');
                if Lia % Ja faz parte da lista
                    coupling_grid(temp_list(nnn),index)=Loc; % armazena o índice
                else % nao faz parte da lista
                    coupling_grid(temp_list(nnn),index)=count; % armazena o indice
                    list2(count,:)=vec; % inclui o vetor de caracteristicas na lista
                    count=count+1;
                end
            end
            
            
            
        end
    end
end
%% Verificacao de Dados em Disco
count=1;
index_to_compute=zeros(0,1);
clear Nc_list
% Verifica se extistem dados de acoplamento armazenados no disco
if isfile('./coupling_data.mat') % Existem
    load('coupling_data.mat','coupling_data'); % faz load dos dados
    
    for nnn=1:size(list2,1) % compara a lista com os dados do load
        [Lia,Loc]=ismember(list2(nnn,:),coupling_data(:,1:24),'rows');
        if Lia % se faz parte da lista
            Nc_list(nnn,:)=coupling_data(Loc,25:33); % armazena o resultados dos acoplamentos
        else %se nao faz parte
            index_to_compute(count)=nnn; % armazena o indice que precisa ser calculado
            count=count+1;
        end
    end
else % se nao exitem dados no disco
    index_to_compute=1:size(list2,1); % indica o calculo de todas as posicoes da lista
end


%% Cálculo dos Tensores Necessários
tic;
fprintf('\nCalculando Tensores de Acoplamento:\n')
if ~compute_PAR % Computaçao Serial
    fprintf('Serial Mode: \n')
    for i=1:size(list2,1) % Itera sobre os membros da lista
        if ismember(i,index_to_compute) % se esse membro é pra ser calculado
            [temp, V2]=write_FileDipolar3D_noDAT(...
                [list2(i,1:4);list2(i,13:16)], ...
                [list2(i,5:8);list2(i,17:20)], [list2(i,9);list2(i,21)],...
                [list2(i,10:12);list2(i,22:24)], platform);
            K = 4*pi*V2/1e-27;
            Nc_list(i,:) = [reshape(temp',1,9)];
            fprintf('%d ',i);
        end
    end
else % computacao Paralela
    fprintf('Parallel Mode: \n')
    parfor i=1:size(list2,1) % Itera sobre a lista
        if ismember(i,index_to_compute) % se o membro precisa ser calculado
            [temp, V2]=write_FileDipolar3D_noDAT(...
                [list2(i,1:4);list2(i,13:16)], ...
                [list2(i,5:8);list2(i,17:20)], [list2(i,9);list2(i,21)],...
                [list2(i,10:12);list2(i,22:24)], platform);
            K = 4*pi*V2/1e-27;
            Nc_list(i,:) = [reshape(temp',1,9)];
            fprintf('%d ',i);
        end
    end
    
end
fprintf('OK! \n')
toc;
%% Cria o Tensor Nc
Nc=sparse(3*part_n,3*part_n);
for i=1:part_n
    for j=1:part_n
        index_i=3*(i-1)+1:3*(i-1)+3;
        index_j=3*(j-1)+1:3*(j-1)+3;
        if coupling_grid(i,j)
            Nc(index_i,index_j)=reshape(Nc_list(coupling_grid(i,j),:),[3,3]);
        else
            Nc(index_i,index_j)=zeros(3,3);
        end
    end
end


%% Salva os dados Novos Computados na Variavel 'coupling_data'
if isfile('./coupling_data.mat')
    coupling_data=[coupling_data;list2(index_to_compute,:) Nc_list(index_to_compute,:)];
else
    coupling_data=[list2 Nc_list];
end
save('coupling_data','coupling_data');


