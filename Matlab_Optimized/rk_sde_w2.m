function [m_new] = rk_sde_w2(m, h_eff, i_s, v, dt, dW)
global alpha alpha_l;

% Drift Term
m_is = cross(m,i_s);

mm_is = cross(m,m_is);

m_heff = cross(m,h_eff);

mm_heff = cross(m,m_heff);

a=-alpha_l*(m_heff+mm_is+alpha*(mm_heff-m_is));%-alpha_l*v.^2.*m; %a

% Difusion Term
m_v = cross(m,v);

mm_v = cross(m,m_v);

b = -alpha_l*m_v-alpha_l*alpha*mm_v; %b

% Auciliary Terms
u      = m + a*dt + b.*dW;
u_plus = m + a*dt + b*sqrt(dt);
u_minus= m + a*dt - b*sqrt(dt);

% 1st Term
m_is = cross(u,i_s);

mm_is = cross(u,m_is);

m_heff = cross(u,h_eff);

mm_heff = cross(u,m_heff);

a_u=-alpha_l*(m_heff+mm_is+alpha*(mm_heff-m_is));%-alpha_l*v.^2.*m; %a

m_1 = m + .5*(a_u+a)*dt;

% 2nd Term
m_v = cross(u_plus,v);

mm_v = cross(u_plus,m_v);

b_plus = -alpha_l*m_v-alpha_l*alpha*mm_v; %b

m_v = cross(u_minus,v);

mm_v = cross(u_minus,m_v);

b_minus = -alpha_l*m_v-alpha_l*alpha*mm_v; %b

m_2 = .25*(b_plus+b_minus+2*b).*dW;

% 3rd Term

m_3 = .25*(b_plus - b_minus).*(dW.^2-dt)/sqrt(dt);

% New m -> m_(i+1)
m_new = m_1 + m_2 + m_3;

end
    

