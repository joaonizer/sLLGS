function [b]=B_term_RK_w2(m, h_eff, i_s, v, dt, dW)
global alpha alpha_l;
% Difusion Term
m_v = cross(m,v);

mm_v = cross(m,m_v);

b = -alpha_l*(m_v+alpha*mm_v);
end