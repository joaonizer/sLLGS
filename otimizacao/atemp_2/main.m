clear all
clc

options = optimoptions(@fmincon,...
    'Display','iter-detailed','Algorithm','interior-point','DiffMinChange',1,'DiffMaxChange',2);

th=10; %nm
%   w   l  dx cortes_y
x0=[0 0 0 0 0 0];
%% Boundaries:
% LB
lb=[-20 -20 -20 -20 -20 -20];
% UB
ub=[0 0 0 0 0 0];

%% 
tic;
[x,fval] = fmincon(@f_obj,x0,...
    [],[],[],[],lb,ub,@myConstraints,options)
toc;

%% Plot Initial and Final
figure();
subplot(1,2,1)
 plot_Final_Particles(x0);
 title('Inicial')
 subplot(1,2,2)
 plot_Final_Particles(x);
 title('Final')
 sdf('P1');