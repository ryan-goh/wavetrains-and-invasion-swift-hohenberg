% solve the Eckhaus eigenvalue value problem for a given mu and k
function [F,J] = SH_1D_eig_JR(uin,p,mesh_params)

  % Rename parameters
  mu = p(1);
  k  = 0.3320;%p(2);
  
  % Auxiliary variables
  n  = mesh_params.nz;
  u  = uin(1:n);
  w  = uin(1+n:2*n);
  wg = uin(2*n+1:3*n);
  wgg= uin(3*n+1:4*n);
  
  c    = uin(4*n+1);
  lam  = uin(4*n+2);
  lamg = uin(4*n+3);
  lamgg= uin(4*n+4);

  u_ref = mesh_params.w0;
  u_refz= mesh_params.w0z;
  int_weights = mesh_params.wz;
  
  Dz = mesh_params.Dz; D2z = mesh_params.D2z; D4z = mesh_params.D4z;
  u4z  = k^4*D4z*u;
  u2z  = k^2*D2z*u;
  uz   = k*Dz*u;
  
  L = -k^4*D4z - 2*k^2*D2z + spdiags((-1+mu)*ones(n,1) - 3*u.^2,0,n,n);

  F1 = -u4z - 2*u2z - u + mu*u - u.^3 + c*uz;
  F2 = L*w   - lam*w; 
  F3 = L*wg  - lamg*w - 4*(speye(n) + k^2*D2z)*k^2*Dz*w;  
  F4 = L*wgg - lamgg*w + 2*lamg*wg + 8*(speye(n) + k^2*D2z)*k^2*Dz*wg + 4*k^2*w + 12*k^2*D2z*w;
  F5 = int_weights*(u_refz.*(u  - u_ref));
  F6 = norm(w,2).^2-1;
  F7 = int_weights*(w.*wg);
  F8 = int_weights*(w.*wgg);

  F = [F1;F2;F3;F4;F5;F6;F7;F8];
  
  

