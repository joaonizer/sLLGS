m1=squeeze(m(:,:,1));
m2=squeeze(m(:,:,2));
Ec=0.5*Ms^2*V(1)*sum(mu0*m1.*(m2*Nc(1:3,4:6)),2)/abs(q);
plot(Ec);
%%
clear m1 m2 Ec
for i=0:2:359
    for j=0:2:359
        m1=[cos(i*pi/180) sin(i*pi/180) 0];
        m2=[cos(j*pi/180) sin(j*pi/180) 0];
        Ec(i/2+1,j/2+1)=0.5*Ms^2*V(1)*sum(mu0*m1.*(m2*Nc(1:3,4:6)),2)/abs(q);
    end
end
surf(0:2:359,0:2:359,Ec')
xlabel('\theta_1')
ylabel('\theta_2')
zlabel('U_{c_1}')