function [px,py]=write_Pontos(w,l,cortes_y,j)
% Create vectors 'px', 'py' e 'd_or'from the inputs:
%     w - width of particle j
%     l - length of particle j
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
end