function [m_new] = rk_sde(m, h_eff, i_s, v, dt, dW)
%% RK Weak 2
a=A_term_RK_w2(m, h_eff, i_s, v, dt, dW);

b=B_term_RK_w2(m, h_eff, i_s, v, dt, dW);

u=m+a*dt+b.*dW;
u_plus=m+a*dt+b*sqrt(dt);
u_minus=m+a*dt-b*sqrt(dt);

a_u=A_term_RK_w2(u, h_eff, i_s, v, dt, dW);

b_u_plus=B_term_RK_w2(u_plus, h_eff, i_s, v, dt, dW);
b_u_minus=B_term_RK_w2(u_minus, h_eff, i_s, v, dt, dW);

m_new = m + 0.5*(a_u+a)*dt...
        + 0.25*(b_u_plus+b_u_minus+2*b).*dW...
        + 0.25*(b_u_plus-b_u_minus).*(dW.^2-dt)/sqrt(dt);
end
    

