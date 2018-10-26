function [a]=A_term_RK_w2(m, h_eff, i_s, v, dt, dW)
global alpha alpha_l;
m_is = cross(m,i_s);

mm_is = cross(m,m_is);

m_heff = cross(m,h_eff);

mm_heff = cross(m,m_heff);

a=-alpha_l*(m_heff+alpha*(mm_heff)-mm_is-v.^2.*m); %a
end