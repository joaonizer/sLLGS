function [px,py,d_or]=write_Pontos(w,l,dx,d_or,cortes_y,j)
% Create vectors 'px', 'py' e 'd_or'from the inputs:
%     w - width of particle j
%     l - length of particle j
%     d - distance between particles j e j-1
%     d_or - vector with center coordinates
%     cortes_y - cuts on the length of the particle
%     j - index of the particle
%%
% Defines x's from w
x1=-w(j)/2;
x2=w(j)/2;
x3=x2;
x4=x1;
% Describes y's from l adjusts cust from cortes_y
y1=l(j)/2;
y2=l(j)/2;
y3=-l(j)/2;
y4=-l(j)/2;
% concatenates 'px' and 'py'
px=[x1 x2 x3 x4];
py=[y1 y2 y3 y4]-cortes_y;

if j==1
    d_or(j,:)=[0 0 0]; % primeira partícula na origem
else
    d_or(j,:)=[w(j-1)/2+w(j)/2+d_or(j-1,1)+dx(j) 0 0]; %
    % sums w/2 from particle j-1
    % and w/2 of particle j
    % plus the acumulated x-distance from the origin
    % plus the distance 'd' between both particles
    % PS: Available only for wires on the x-direction
end

