function [Ec] = compute_Ec(m,Nc,V,kt)
global Ms mu0 q;
n=size(m,2);
%k=zeros(n,1);
k=-V*Ms^2*mu0/2; % Constant to Scale energy and convert J to eV

for i=1:n
    for j=1:n
       Ec(i,j) = m(:,i)'*Nc(:,:,i,j)*m(:,j)*k(j)/kt;
        
    end
end
        


end