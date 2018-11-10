%%
figure;

i=11;
    Nc2=Nc;
    Nc2(Nc2<10^(-i))=0;
    Nc2=sparse(Nc2);
    mem=whos('Nc2');
    un(i)=mem.bytes;
nz(i)=nnz(Nc2);
spy(Nc2)
%hold on
%spy(Nc<0.000001,'ro')
xlabel('Índice i')
ylabel('Índice j')
title([num2str(nz(i)) ' Elementos > 10^{-' num2str(i) '}'])
sdf('P1')

%%
close all
figure;
spy(Nd)
xlabel('Índice i')
ylabel('Índice j')
nz=sum(sum(Nd==0));
title([num2str(nz) ' Elementos Zero'])
sdf('P1')
figure;
for i=6:16
nz(i-5)=sum(sum(Nd<10^(-i)));
end
semilogx(10.^(-[6:16]),nz,'o-')
xlabel('Threshold')
ylabel('Quantidade de Elementos Zero')
grid on
sdf('P1')
%%
clear nz Ms Mc un
figure;
spy(Nc)
%hold on
%spy(Nc<0.000001,'ro')
xlabel('Índice i')
ylabel('Índice j')
for i=1:16
    Nc2=Nc;
    Nc2(Nc2<10^(-i))=0;
    Nc2=sparse(Nc2);
    mem=whos('Nc2');
    un(i)=mem.bytes;
nz(i)=48*48-nnz(Nc2);
end
title([num2str(48*48-nz(6)) ' Elementos < 10^{-6}'])
sdf('P1')
figure;
semilogx(10.^(-[1:16]),nz,'o-')
xlabel('Threshold')
ylabel('Qtd de Elementos Não Zero')
grid on
sdf('P1')

for j=1:length(nz)
    Mc(j)=48*48*64/8;
    Ms(j)=nz(j)*16+nz(j)*16;
end
%end
figure
semilogx(10.^(-[1:16]),Mc/1000,'.-');
hold on
semilogx(10.^(-[1:16]),un/1000,'-o');
xlabel('Threshold')
ylabel('Memória (kB)')
grid on
legend('Cheia','Esparsa')
sdf('P1')
