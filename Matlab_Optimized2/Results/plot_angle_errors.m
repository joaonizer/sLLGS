% i=2;
% plot_Energias(squeeze(m(:,:,i)), squeeze(hT(:,:,i)),...
%     squeeze(i_s(:,:,i)), squeeze(hc(:,:,i)), Nd(1:3,1:3),V(i),t,50,100)
%%
for j=1:2:14
    for i=1:6
        if ~mod(i,2)
            error(j,i)=sum(angles(1,j,:,i)<-80)/250;
        else
            error(j,i)=sum(angles(1,j,:,i)>80)/250;
        end
    end
end
for j=1:2:14
    plot(1:6,error(j,:),'-o')
    hold on
end
%1 2 5 6 10 11 15 16 20 21 25 26 30 31
%legend('1nm','2nm','5nm','6nm','10nm','11nm','15nm','16nm','20nm','21nm','25nm','26nm','30nm','31nm')
legend('1nm','5nm','10nm','15nm','20nm','25nm','30nm')

sdf('P1')