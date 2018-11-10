clear Ms Mc
% n=1:12;
% Ms=zeros(12,length(1:2:2^12));
% Mc=Ms;
%for i=1:n
for j=0:48*48
    Mc(j+1)=48*48*64/8;
    Ms(j+1)=j*16+32*8;
end
%end

plot([0:512]/512,Ms./Mc);
hold on
plot(Ms)
