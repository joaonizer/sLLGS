[x,y,z]=sphere;
f1=mesh(x,y,z,'Marker','.',...
    'FaceColor','none','EdgeColor',[189 191 193]/255);
hold on
f2=plot3(squeeze(m(1:25:end,1,1)),...
    squeeze(m(1:25:end,2,1)),...
    squeeze(m(1:25:end,3,1)),...
    'Linewidth',2);
xlabel('m_x');
ylabel('m_y');
zlabel('m_z');
grid on;
sdf('P1')
pbaspect([1,1,1]);
set(f1,'Linewidth',0.1)
set(f2,'Linewidth',1.5)