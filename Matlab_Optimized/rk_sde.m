function [m_new] = rk_sde(m, h_eff, i_s, v, dt, dW)
global alpha alpha_l;

% Drift Term
% m_is = cross(m,i_s);
% 
% mm_is = cross(m,m_is);
% 
% m_heff = cross(m,h_eff);
% 
% mm_heff = cross(m,m_heff);
% 
% a=-alpha_l*(m_heff+alpha*(mm_heff)-mm_is)-alpha_l*v.^2.*m; %a
% 
% % Difusion Term
% m_v = cross(m,v);
% 
% mm_v = cross(m,m_v);
% 
% b = -alpha_l*m_v-alpha_l*alpha*mm_v; %b
% 
% % Upsilon term
% Y = m + a*dt + b*sqrt(dt);
% 
% % b(Y)
% Y_v = cross(Y,v);
% 
% YY_v = cross(Y,Y_v);
% 
% bY = -alpha_l*Y_v-alpha_l*alpha*YY_v; %bY
% 
% % Update on 'm'
% m_new = (m + a*dt + b.*dW ...
%         + 0.5*(bY-b).*(dW.^2-dt)/sqrt(dt));
%% 
m_is = cross(m,i_s);

mm_is = cross(m,m_is);

m_heff = cross(m,h_eff);

mm_heff = cross(m,m_heff);

a=-alpha_l*(m_heff+alpha*(mm_heff)-mm_is)+alpha_l*v.^2.*m; %a

% Difusion Term
m_v = cross(m,v);

mm_v = cross(m,m_v);

b = -alpha_l*m_v-alpha_l*alpha*mm_v;

m_new = m + a*dt + b.*dW;

%% RK4-Heun
% % d1
% m_is = cross(m,i_s);
% 
% mm_is = cross(m,m_is);
% 
% m_heff = cross(m,h_eff);
% 
% mm_heff = cross(m,m_heff);
% 
% d1=-alpha_l*(m_heff+mm_is+alpha*(mm_heff-m_is));% -alpha_l*v.^2.*m; %a
% 
% 
% % s1
% m_v = cross(m,v);
% 
% mm_v = cross(m,m_v);
% 
% s1 = -alpha_l*m_v-alpha_l*alpha*mm_v; %b
% 
% % d2
% m_is = cross(m+(d1+s1)*dt/2,i_s);
% 
% mm_is = cross(m+(d1+s1)*dt/2,m_is);
% 
% m_heff = cross(m+(d1+s1)*dt/2,h_eff);
% 
% mm_heff = cross(m+(d1+s1)*dt/2,m_heff);
% 
% d2=-alpha_l*(m_heff+mm_is+alpha*(mm_heff-m_is));% -alpha_l*v.^2.*m; %a
% 
% 
% % s2
% m_v = cross(m+d1*dt+s1*sqrt(dt).*dW,v);
% 
% mm_v = cross(m+d1*dt+s1*sqrt(dt).*dW,m_v);
% 
% s2 = -alpha_l*m_v-alpha_l*alpha*mm_v; %b
% 
% 
% % d3
% m_is = cross(m+(d2+s1)*dt/2,i_s);
% 
% mm_is = cross(m+(d2+s1)*dt/2,m_is);
% 
% m_heff = cross(m+(d2+s1)*dt/2,h_eff);
% 
% mm_heff = cross(m+(d2+s1)*dt/2,m_heff);
% 
% d3=-alpha_l*(m_heff+mm_is+alpha*(mm_heff-m_is));% -alpha_l*v.^2.*m; %a
% 
% % d4
% m_is = cross(m+(d3+s1)*dt,i_s);
% 
% mm_is = cross(m+(d3+s1)*dt,m_is);
% 
% m_heff = cross(m+(d3+s1)*dt,h_eff);
% 
% mm_heff = cross(m+(d3+s1)*dt,m_heff);
% 
% d4=-alpha_l*(m_heff+mm_is+alpha*(mm_heff-m_is));% -alpha_l*v.^2.*m; %a
% 
% S = (s1+s2)/2;
% D = (d1+2*d2+2*d3+d4)/6;
% 
% m_new = m + D*dt + S*sqrt(dt).*dW;

end
    

