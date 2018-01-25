epsilon = 1;

while epsilon > 1e-6



J = I - dt/2*( (m+m_n)/2*Jh );


epsilon= sum(m(n+1,:)-m(n,:) -dt*f(m(n,:),m(n+1,:),h_eff(n,:),h_eff(n+1,:),dt));
end

function [ff]=f(m_n,m,h_eff_n,h_eff,dt)
    m12=(m_n+m)/2;
    
    h_eff12=(h_eff_n+h_eff)/2;
    
    ff =  cross(-m12,h_eff12)         ...
        + alpha*cross(m12,(m-m_n)/dt) ...
        - cross(-m12,cross(m12,i_s));
end