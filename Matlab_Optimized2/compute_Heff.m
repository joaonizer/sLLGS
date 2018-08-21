function [h_eff]=compute_Heff(m,h_app,hT, hc,Nd)
% Funcao para o calculo do campo efetivo H_eff:
%
global HkMs n;

h_eff = h_app ...           % Campo externo aplicado (adimensional)
    +HkMs*dot(n,m)*n ...	% Anistropia magnetocristalina (adimensional)
    -m*Nd' ...              % Campo desmagnetizacao (adimensional)
    +hT...                  % Campo termico (adimensional)
    +hc;                    % Campo de Acoplamento (adimensional)
end