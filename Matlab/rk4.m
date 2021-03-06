function [m_new]=rk4(m,h_eff,i_s,dt)
% Implementation of the Range-Kuta method
% For solving the s-LLGS
%   m - particles magnetization at time instant i
%   h_eff - effective field applied to the particle at instant i
%   i_s - spin current applied to the particle
%   dt - normalized time step

% Step 1
k(1,:)=dm(m,h_eff,i_s);
mm(1,:)=m+k(1,:)*dt/2;

% Step 2
k(2,:)=dm(mm(1,:),h_eff,i_s);
mm(2,:)=m+k(2,:)*dt/2;

% Step 3
k(3,:)=dm(mm(2,:),h_eff,i_s);
mm(3,:)=m+k(3,:)*dt;

% Step 4
k(4,:)=dm(mm(3,:),h_eff,i_s);

% New Magnetization
m_new(:)=m+(k(1,:) + 2*(k(2,:) + k(3,:)) + k(4,:))*dt/6;

end
