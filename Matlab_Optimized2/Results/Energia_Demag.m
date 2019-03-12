%% Total
figure;
clear Ed1 Ed2 Ec1 Ec2 Et theta1 theta2
theta1=0:0.1:2*pi;
theta2=0:0.1:2*pi;
p1=1; ip1=3*(p1-1)+1:3*(p1-1)+3;
p2=2; ip2=3*(p2-1)+1:3*(p2-1)+3;
for i=1:length(theta1)
    for j=1:length(theta2)
        M1=[cos(theta1(i)),sin(theta1(i)),0];
        M2=[cos(theta2(j)),sin(theta2(j)),0];
        Ed1(i,j)=sum(0.5*V(p1)*mu0*Ms^2*M1*Nd(ip1,ip1)*M1');
        Ed2(i,j)=sum(0.5*V(p2)*mu0*Ms^2*M2*Nd(ip2,ip2)*M2');
        Ec1(i,j)=sum(0.5*V(p1)*mu0*Ms^2*M1*Nc(ip1,ip2)*M2');
        Ec2(i,j)=sum(0.5*V(p2)*mu0*Ms^2*M2*Nc(ip2,ip1)*M1');
        
    end
end
Et=Ed1+Ed2+Ec1+Ec2;

plot(theta2(1:end)*180/pi,Ed2(ceil(90*pi/180/0.1),:)'/q,'--')
hold on
plot(theta2(1:end)*180/pi,Ec2(ceil(90*pi/180/0.1),:)'/q,'-.')
plot(theta2(1:end)*180/pi,Ed2(ceil(90*pi/180/0.1),:)'/q+Ec2(ceil(90*pi/180/0.1),:)'/q)
xticks([0 90 180 270 360])
xticklabels({'0','90','180', '270', '360'})
xlabel('\theta_2 [degree]')
ylabel('U [eV]')
xlim([0 360])
legend('Demagnetizing','Coupling', 'Total')
sdf('P1')

%%
%subplot(1,2,1)
surf(theta1(1:end)*180/pi,theta2(1:end)*180/pi,Et'/q)
xticks([0 90 180 270 360])
xticklabels({'0','90','180', '270', '360'})
yticks([0 90 180 270 360])
yticklabels({'0','90','180', '270', '360'})
xlabel('\theta_1 [degree]')
ylabel('\theta_2 [degree]')
zlabel('U_t [eV]')
xlim([0 360])
ylim([0 360])
hcb=colorbar;
title(hcb,'U_t [eV]')
sdf('P1')


%% OR AND - polarizado
figure;
clear Ed Ec Et theta1 theta2 theta3 theta4
theta1=-pi/2; %AND=pi/2 OR = -pi/2;
theta2=[3*pi/2 pi/2 3*pi/2 pi/2];
theta3=0:0.1:2*pi-0.1;
theta4=[3*pi/2 3*pi/2 pi/2 pi/2];
p1=1; ip1=3*(p1-1)+1:3*(p1-1)+3;
p2=3; ip2=3*(p2-1)+1:3*(p2-1)+3;
p3=4; ip3=3*(p3-1)+1:3*(p3-1)+3;
p4=5; ip4=3*(p4-1)+1:3*(p4-1)+3;

for j=1:4
    for i=1:length(theta3)
        M1=[cos(theta1),sin(theta1),0];
        M2=[cos(theta2(j)),sin(theta2(j)),0];
        M3=[cos(theta3(i)),sin(theta3(i)),0];
        M4=[cos(theta4(j)),sin(theta4(j)),0];
        Ed(i,j)=sum(0.5*V(p3)*mu0*Ms^2*M3*Nd(ip3,ip3)*M3');
        Ec(i,j)=sum(0.5*V(p3)*mu0*Ms^2*M3*Nc(ip2,ip3)*M2')...
            +sum(0.5*V(p3)*mu0*Ms^2*M3*Nc(ip4,ip3)*M4')...
            +sum(0.5*V(p3)*mu0*Ms^2*M3*Nc(ip1,ip3)*M1');
        Et(i,j)=Ed(i,j)+Ec(i,j);
    end
end
plot(theta3(1:end)*180/pi,Et(:,1)/q,'-','MarkerIndices',1:2:length(theta3))
hold on
plot(theta3(1:end)*180/pi,Et(:,2)/q,'--o','MarkerIndices',1:2:length(theta3))
plot(theta3(1:end)*180/pi,Et(:,3)/q,'--x','MarkerIndices',1:3:length(theta3))
plot(theta3(1:end)*180/pi,Et(:,4)/q,'--','MarkerIndices',1:2:length(theta3))
xticks([0 90 180 270 360])
xticklabels({'0','90','180', '270', '360'})
xlabel('\theta_O [graus]')
ylabel('U_r [eV]')
legend('\downarrow\downarrow','\downarrow\uparrow','\uparrow\downarrow','\uparrow\uparrow')
grid on
sdf('P1')