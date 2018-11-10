%%
theta2=0:0.1:2*pi-0.1;
clear Ed
for j=1:length(theta2)
    M1=[cos(theta2(j)),sin(theta2(j)),0];
    Ed(j)=sum(0.5*V(1)*mu0*Ms^2*M1*Nd(4:6,4:6)*M1');
end
plot(theta2*180/pi,Ed/q,'o','MarkerIndices',1:1:length(theta2))
xticks([0 90 180 270 360])
xticklabels({'0','90','180', '270', '360'})
xlabel('\theta [graus]')
ylabel('E_d [eV]')
sdf('P1')

%%
figure;
theta1=pi/2;
theta2=0:0.1:2*pi-0.1;
clear Ec
p1=1;
p2=4;
for i=1:length(theta2)
    M1=[cos(theta1),sin(theta1),0];
    M2=[cos(theta2(i)),sin(theta2(i)),0];
    Ec(i)=sum(0.5*V(p2)*mu0*Ms^2*M1*Nc([1:3]*p1,[1:3]*p2)*M2');
end
plot(theta2(1:end)*180/pi,Ec/q,'o','MarkerIndices',1:1:length(theta2))
xticks([0 90 180 270 360])
xticklabels({'0','90','180', '270', '360'})
xlabel('\theta_2 [graus]')
ylabel('E_c [eV]')
sdf('P1')
%% Total
figure;
clear Ec
theta1=0:0.1:2*pi-0.1;
theta2=0:0.1:2*pi-0.1;
p1=1;
p2=4;
for i=1:length(theta1)
    for j=1:length(theta2)
        M1=[cos(theta1(i)),sin(theta1(i)),0];
        M2=[cos(theta2(j)),sin(theta2(j)),0];
        Ed1=sum(0.5*V(p1)*mu0*Ms^2*M1*Nd([1:3]*p1,[1:3]*p1)*M1');
        Ed2=sum(0.5*V(p2)*mu0*Ms^2*M2*Nd([1:3]*p2,[1:3]*p2)*M2');
        Ec1=sum(0.5*V(p1)*mu0*Ms^2*M1*Nc([1:3]*p1,[1:3]*p2)*M2');
        Ec2=sum(0.5*V(p2)*mu0*Ms^2*M2*Nc([1:3]*p2,[1:3]*p1)*M1');
        Et(i,j)=Ed1+Ed2+Ec1+Ec2;
    end
end
%subplot(1,2,1)
surf(theta1(1:end)*180/pi,theta2(1:end)*180/pi,Et'/q)
xticks([0 90 180 270 360])
xticklabels({'0','90','180', '270', '360'})
yticks([0 90 180 270 360])
yticklabels({'0','90','180', '270', '360'})
xlabel('\theta_1 [graus]')
ylabel('\theta_2 [graus]')
zlabel('E_t [eV]')
sdf('P1')

% subplot(1,2,2)
% surf(theta1(1:end)*180/pi,theta2(1:end)*180/pi,Et'/q)
% xticks([0 90 180 270 360])
% xticklabels({'0','90','180', '270', '360'})
% yticks([0 90 180 270 360])
% yticklabels({'0','90','180', '270', '360'})
% xlabel('\theta_1 [graus]')
% ylabel('\theta_2 [graus]')
% zlabel('E_t [eV]')
% sdf('P1')
%% AND - shape
figure;
clear Ed Ec Et theta1 theta2 theta3 theta4
theta1=[3*pi/2 pi/2 3*pi/2 pi/2];
theta2=0:0.01:2*pi-0.01;
theta3=[3*pi/2 3*pi/2 pi/2 pi/2];
p1=9;
p2=10;
p3=11;
for j=1:4
    for i=1:length(theta2)
        M1=[cos(theta1(j)),sin(theta1(j)),0];
        M2=[cos(theta2(i)),sin(theta2(i)),0];
        M3=[cos(theta3(j)),sin(theta3(j)),0];
        Ed(i,j)=sum(0.5*V(p2)*mu0*Ms^2*M2*Nd([1:3]*p2,[1:3]*p2)*M2');
        Ec(i,j)=sum(0.5*V(p1)*mu0*Ms^2*M2*Nc([1:3]*p2,[1:3]*p1)*M1')...
            +sum(0.5*V(p3)*mu0*Ms^2*M2*Nc([1:3]*p2,[1:3]*p3)*M3');
        Et(i,j)=Ed(i,j)-Ec(i,j);
    end
end
plot(theta2(1:end)*180/pi,Et(:,1)/q,'-','MarkerIndices',1:10:length(theta2))
hold on
plot(theta2(1:end)*180/pi,Et(:,2)/q,'--o','MarkerIndices',1:20:length(theta2))
plot(theta2(1:end)*180/pi,Et(:,3)/q,'--x','MarkerIndices',1:30:length(theta2))
plot(theta2(1:end)*180/pi,Et(:,4)/q,'--','MarkerIndices',1:40:length(theta2))
xticks([0 90 180 270 360])
xticklabels({'0','90','180', '270', '360'})
xlabel('\theta_O [graus]')
ylabel('U_r [eV]')
legend('\downarrow\downarrow','\downarrow\uparrow','\uparrow\downarrow','\uparrow\uparrow')
grid on
sdf('P1')
%% OR - shape
figure;
clear Ed Ec Et theta1 theta2 theta3 theta4
theta1=[3*pi/2 pi/2 3*pi/2 pi/2];
theta2=0:0.1:2*pi-0.1;
theta3=[3*pi/2 3*pi/2 pi/2 pi/2];

for j=1:4
    for i=1:length(theta2)
        M1=[cos(theta1(j)),sin(theta1(j)),0];
        M2=[cos(theta2(i)),sin(theta2(i)),0];
        M3=[cos(theta3(j)),sin(theta3(j)),0];
        Ed(i,j)=sum(0.5*V(2)*mu0*Ms^2*M2*Nd([1:3]*p2,[1:3]*p2)*M2');
        Ec(i,j)=sum(0.5*V(2)*mu0*Ms^2*M2*Nc([1:3]*p2,[1:3]*p1)*M1')...
            +sum(0.5*V(2)*mu0*Ms^2*M2*Nc([1:3]*p2,[1:3]*p3)*M3');
        Et(i,j)=Ed(i,j)+Ec(i,j);
    end
end
plot(theta2(1:end)*180/pi,Et(:,1)/q,'-','MarkerIndices',1:2:length(theta2))
hold on
plot(theta2(1:end)*180/pi,Et(:,2)/q,'--o','MarkerIndices',1:2:length(theta2))
plot(theta2(1:end)*180/pi,Et(:,3)/q,'--x','MarkerIndices',1:3:length(theta2))
plot(theta2(1:end)*180/pi,Et(:,4)/q,'--','MarkerIndices',1:2:length(theta2))
xticks([0 90 180 270 360])
xticklabels({'0','90','180', '270', '360'})
xlabel('\theta_O [graus]')
ylabel('U_r [eV]')
legend('\downarrow\downarrow','\downarrow\uparrow','\uparrow\downarrow','\uparrow\uparrow')
grid on
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