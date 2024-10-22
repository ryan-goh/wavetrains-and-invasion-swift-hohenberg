function [F,J] = SH_1D_cubic_zig_zag(uin,p,mesh_params)

  % Rename parameters
  mu = p(1);
  
  % Auxiliary variables
  n = mesh_params.nz;
  u = uin(1:n);
  c = uin(n+1);
  k = uin(n+2);

  uz   = k*mesh_params.Dz*u;
  uzz  = k^2*mesh_params.D2z*u;
  uzzzz= k^4*mesh_params.D4z*u;
  int_weights = mesh_params.wz;
  
  F = -uzzzz - 2*uzz - u + mu*u - u.^3 + c*uz;
  F(n+1) = int_weights*(mesh_params.w0z.*(u  - mesh_params.w0)); % phase condition
  F(n+2) = int_weights*( ((k^2*mesh_params.D2z + speye(n))*uz).*uz )./(mesh_params.wz*(uz.^2)); % zig-zag eigenvalue = 0
     
end
