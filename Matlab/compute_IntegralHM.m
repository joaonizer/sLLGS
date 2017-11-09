function [A]=compute_IntegralHM(H,m,V)
% Calcula as areas das curvas HxM em Joules
% H - Campo externo aplicado em T - kg/A.s2
% m - Magnetizacao da particula A/m
% V - Volume da particula em m^3
% dt - Passo de tempo da simulacao
%

global kbT;
A=zeros(3,1);
%  figure;
%  subplot(1,2,1)
% plot(m(2:N,:)-m(1:N-1,:));
%  subplot(1,2,2)
% plot(h_app(2:N,:)-h_app(1:N-1,:));
% for i=1:length(H)/2
%     A=A+abs(m(i,:)-m(i+1,:)).*abs(H(i,:)-H(i+1,:));
% end
% for i=length(H)/2:length(H)/2-1
%     A=A-abs(m(i,:)-m(i+1,:)).*abs(H(i,:)-H(i+1,:));
% end
A(1)=abs(sum(trapz(H(:,1),m(:,1))));
A(2)=abs(sum(trapz(H(:,2),m(:,2))));
A(3)=abs(sum(trapz(H(:,3),m(:,3))));
A=A*V;
fprintf('\nLandauer''s Limit:   %e J\n', kbT*log(2));
fprintf('Area x: %e J\n',A(1));
fprintf('Area y: %e J\n',A(2));
fprintf('Area z: %e J\n',A(3));
fprintf('TOTAL: %e J \n %.2f%% of Landauer''s Limit\n',sum(A), 100*sum(A)/kbT/log(2));

end