function maxre=max_essential_spec_SH_cubic_JR(u_ini,p,mesh_params,sigmax)
% return maximal real part of eigenvalue of the Floquet operator for the cubic SH eqn
% B(sigma;mu,k) = -(1 + k^2*(dz + i*sigma)^2)^2 + mu - 3*up^2

mu = p(1);
c = 0;
k = p(2); 

p(1) = mu;
p(2) = k;
nz = mesh_params.nz;
up = u_ini(1:nz);

%sigmax = linspace(-k/2,k/2,501);
DN = @(u) -3*u.^2;
maxre=-1e15;
for i = 1:length(sigmax)
        sigma = 1i*sigmax(i);
        L = -((k*mesh_params.Dz + sigma*speye(nz))^2+(1)*speye(nz))^2 ...
        + spdiags(ones(nz,1)*mu,0,nz,nz)  ...
        + spdiags(DN(up),0,nz,nz);
    
        d = eig(full(L)); maxrei=max(real(d));
        if (maxrei>maxre) maxre=maxrei; end
end