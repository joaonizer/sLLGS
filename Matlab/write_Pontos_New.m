function [px,py,d_or]=write_Pontos_New(w,l,r,dx,dy,grid,xshift,yshift)
% Create vectors 'px', 'py' e 'd_or'from the inputs:
%     w - width of particle j
%     l - length of particle j
%     d - distance between particles j e j-1
%     d_or - vector with center coordinates
%     cortes_y - cuts on the length of the particle
%     j - index of the particle
%%
cortes_y=[0 0 0 0];
% Defines x's from w
x1=-w/2;
x2=w/2;
x3=x2;
x4=x1;
% Describes y's from l adjusts cust from cortes_y
y1=l/2;
y2=l/2;
y3=-l/2;
y4=-l/2;
% concatenates 'px' and 'py'
px_r=[x1 x2 x3 x4];
py_r=[y1 y2 y3 y4]-cortes_y;
d_or=zeros(size(grid,2)*size(grid,1),3);
px=zeros(size(grid,2)*size(grid,1),4);
py=px;
for i=1:size(grid,1) % iterate over row
    for j=1:size(grid,2) % iterate over columns
        index = (i-1)*size(grid,2) + j;
        switch grid(i,j)
            case 0 % No particle
                d_or(index,:)=zeros(1,3);
                px(index,:) =zeros(1,4);
                py(index,:) = zeros(1,4);
            case 1 % Rectangle
                d_or(index,:)=[(j-1)*(w+dx)+xshift(i,j) (i-1)*(l+dy)+yshift(i,j) 0];
                px(index,:) = px_r;
                py(index,:) = py_r;
            case 2 % Circle
                d_or(index,:)=[(j-1)*(w+dx)+xshift(i,j)+w/2 (i-1)*(l+dy)+yshift(i,j)+l/2 0];
                px(index,:) = [r 0 0 0];
                py(index,:) = [r 0 0 0];
            otherwise % set as rectangle
                d_or(index,:)=zeros(1,3);
                px(index,:) =zeros(1,4);
                py(index,:) = zeros(1,4);
        end
    end
end

end

