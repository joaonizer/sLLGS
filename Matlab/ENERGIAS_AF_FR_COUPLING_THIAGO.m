clc
Nc1=squeeze(Nc(:,:,1,[2 3 6 7 8 11 12 13]));
q = 1.60217662e-19; % carga do eletron C
mu0=4*pi*1e-7;      % H/m ou T.m/A
Ms=796000;      % Saturation Magnetization A/m
V = 50*100*10*1e-27;
m0=[0 1 0];
m=zeros(3,8);
m(:,1:2)=[
    0 -1 0
    0 -1 0]';
m(:,3:8)=[
    0 1 0
    0 1 0
    0 1 0
    0 1 0
    0 1 0
    0 1 0
    ]';
max_yy=max(max(abs(Nc(2,2,:))));
for i =1:length(m)
    hc = -squeeze(m(:,i))'*squeeze(Nc1(:,:,i));
    temp=-mu0*V*Ms^2*m0.*hc/q;
    Ec(i)= sum(temp,2);
    fprintf('& %.2f\t& %.2f\t& %.2f\n',Ec(i),Nc1(2,2,i)*1000,-Nc1(2,2,i)/max_yy);
end