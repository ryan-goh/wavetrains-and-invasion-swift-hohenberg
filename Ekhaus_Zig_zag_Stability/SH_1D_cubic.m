function [F,J] = SH_1D_cubic(uin,p,mesh_params)

  % Rename parameters
  mu = p(1);
  k  = p(2);
  
  % Auxiliary variables
  n = mesh_params.nz;
  u = uin(1:n);
  c = uin(n+1);

  uz   = k*mesh_params.Dz*u;
  uzz  = k^2*mesh_params.D2z*u;
  uzzzz= k^4*mesh_params.D4z*u;
  
  F = -uzzzz - 2*uzz - u + mu*u - u.^3 + c*uz;
  
  loc = mesh_params.w0z;
  F(n+1) =  mesh_params.wz*(loc.*(u  - mesh_params.w0));
  
  % Jacobian
  if nargout > 1
     J = sparse(n,n);

     J = -k^4*mesh_params.D4z - 2*k^2*mesh_params.D2z + spdiags((-1+mu)*ones(n,1) - 3*u.^2,0,n,n) ...
         + c*k*mesh_params.Dz;
      
     J(n+1,1:n)  =mesh_params.wz*spdiags(loc,0,n,n);

     epsiF = 1e-8;
     dF = SH_1D_cubic([uin(1:n); uin(n+1) + epsiF],p,mesh_params,1);
     J(:,n+1) = (dF - F)/epsiF;

  end
      

end
