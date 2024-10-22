% solve the eigenvalue value problem 
function [F,J] = SH_1D_eig_x(uin,p,mesh_params)

  % Rename parameters
  sigma = p(1);
  k  = p(2);
  mu = p(3);

  
  % Auxiliary variables
  n  = mesh_params.nz;
  u  = uin(1:n);
  w  = uin(1+n:2*n);
  
  c    = uin(2*n+1);
  lam  = uin(2*n+2);

  u_ref = mesh_params.w0;
  u_refz= mesh_params.w0z;
  int_weights = mesh_params.wz;
  
  Dz = mesh_params.Dz; D2z = mesh_params.D2z; D4z = mesh_params.D4z;
  u4z  = k^4*D4z*u;
  u2z  = k^2*D2z*u;
  uz   = k*Dz*u;
  
  B = -(k^2*(Dz+ 1i*sigma*speye(n))^2+speye(n))^2 + ... 
      spdiags(mu*ones(n,1) - 3*u.^2,0,n,n);

  F1 = -u4z - 2*u2z - u + mu*u - u.^3 + c*uz;
  F2 = B*w   - lam*w; 
  F3 = int_weights*(u_refz.*(u  - u_ref));
  F4 = norm(w,2).^2-1;
  
  F = [F1;F2;F3;F4];
  
  

