Nc_new=Nc;
% Nc_new(Nd>1.5)=0;
clc
%% Rectangle
dim='xyzxyzxyz';
NMC=[0.5 1 5 10];
fprintf('Rectangular');
fprintf('\nAmostras\tDim.\tMédia\tVar.\t2*SD\t2SD/Media\n');
for j=1:4
    for i=1:3
        rec(1,j,i)=mean(Nc_new(j,:,i,i+3));
        rec(2,j,i)=var(Nc_new(j,:,i,i+3));
        fprintf('%1.1fM\t%c%c\t%f\t%f\t%f\t%2.2f%%\n',NMC(j),dim(i),dim(i),rec(1,j,i), ...
            (rec(2,j,i)),2*sqrt(rec(2,j,i)),2*sqrt(rec(2,j,i))/rec(1,j,i)*100);
    end
end
% %% Slanted Up
% Nd_new=Nd;
% fprintf('\n\nSlanted Up 20nm');
% fprintf('\nAmostras\tDim.\tMédia\tVar.\t2*SD\t2SD/Media\n');
% for j=1:4
%     for i=4:6
%         rec(1,j,i)=mean(Nd_new(j,:,i,i));
%         rec(2,j,i)=var(Nd_new(j,:,i,i));
%         fprintf('%1.1fM\t%c%c\t%f\t%f\t%f\t%2.2f%%\n',NMC(j),dim(i),dim(i),rec(1,j,i), ...
%             (rec(2,j,i)),2*sqrt(rec(2,j,i)),2*sqrt(rec(2,j,i))/rec(1,j,i)*100);
%     end
% end
% %% Slanted Up and Down
% Nd_new=Nd;
% fprintf('\n\nSlanted Up and Down 20nm');
% fprintf('\nAmostras\tDim.\tMédia\tVar.\t2*SD\t2SD/Media\n');
% for j=1:4
%     for i=7:9
%         rec(1,j,i)=mean(Nd_new(j,:,i,i));
%         rec(2,j,i)=var(Nd_new(j,:,i,i));
%         fprintf('%1.1fM\t%c%c\t%f\t%f\t%f\t%2.2f%%\n',NMC(j),dim(i),dim(i),rec(1,j,i), ...
%             (rec(2,j,i)),2*sqrt(rec(2,j,i)),2*sqrt(rec(2,j,i))/rec(1,j,i)*100);
%     end
% end