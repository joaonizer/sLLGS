function [dm]=dm(m,h_eff,i_s)
% Function to compute the dm required for the RK4 algorithm
% RHS of the normalized s-LLGS
%   m - magnetization
%   h_eff - Effective field
%   i_s = spin current

global alpha alpha_l;

m_is = cross(m,i_s);

mm_is = cross(m,m_is);

m_heff = cross(m,h_eff);

mm_heff = cross(m,m_heff);

dm=-alpha_l*(m_heff+mm_is+alpha*(mm_heff-m_is));

end