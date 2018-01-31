function [dm]=dm_difusion(m,hT)
% Function to compute the drift term required for the RK4 algorithm
% RHS of the normalized s-LLGS
%   m - magnetization
%   h_eff - Effective field
%   i_s = spin current

global alpha alpha_l;

m_hT = cross(m,hT);

mm_hT = cross(m,m_hT);

dm = -alpha_l*m_hT-alpha_l*alpha*mm_hT;

end