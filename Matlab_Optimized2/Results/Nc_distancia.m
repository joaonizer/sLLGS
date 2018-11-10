%% Caso ferromagnético
n=50;
part_n=2;
compute_PAR=0;
plataform='win';
px=[-25    25    25   -25
   -25    25    25   -25];
py=[50    50   -50   -50
    50    50   -50   -50];
th=[5     5];
for i=1:n
dy=i;
d_or=[0   100+dy     0
     0     0     0];
Nc=compute_Nc2(px,py,th,d_or,part_n,compute_PAR,platform);
Nc_ferro(i)=Nc(2,5);
end
%%
part_n=2;
compute_PAR=0;
plataform='win';
px=[-25    25    25   -25
   -25    25    25   -25];
py=[50    50   -50   -50
    50    50   -50   -50];
th=[5     5];
for i=1:n
dx=i;
d_or=[0   0     0
     50+dx     0     0];
Nc=compute_Nc2(px,py,th,d_or,part_n,compute_PAR,platform);
Nc_antiferro(i)=Nc(2,5);
end
%% Plot
plot(1:n,abs(Nc_ferro));
hold on
plot(1:n,abs(Nc_antiferro));
plot([10:24],[ones(1,15)*Nc_antiferro(10)],'--')
xlabel('dx,dy (nm)')
ylabel('|Nc_{yy}|')
legend('Ferromagnético','Antiferromagnético');
sdf('P1')