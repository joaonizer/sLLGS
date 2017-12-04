clear all
clc

options = optimoptions(@fmincon,...
    'Display','iter-detailed','Algorithm','sqp','DiffMinChange',1,'DiffMaxChange',2);

th=10; %nm
%   w   l  dx cortes_y
x0=[50 100 10 0 0 0 0 ...
    50 100 24 0 0 0 0 ];
%% Boundaries:
% LB
lb=[ 50  90 10 -20 -20 -20 -20 50  90 24 -20 -20 -20 -20 ];
% UB
ub=[ 60 100 10  20  20  20  20 60 100 24  20  20  20  20 ];

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