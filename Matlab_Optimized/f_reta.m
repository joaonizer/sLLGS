function [yup,ydown]=f_reta(x,px,py)
% Calcula as retas superiores e inferiores para
% gerar a imagem .gif para simulacao no OOMMF
% Utilizada por 'particle_OOMMF'
%% Computa as retas
m1=(py(2)-py(1))/(px(2)-px(1)); % slope da superior
m2=(py(3)-py(4))/(px(3)-px(4)); % slope da inferior
yup=m1*(x-px(1))+py(1);         % Valores de Y superior
ydown=m2*(x-px(4))+py(4);       % Valores de Y inferior

end