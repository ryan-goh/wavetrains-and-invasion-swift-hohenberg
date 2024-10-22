function F = SolutionMeasures_cubic(step,u,p,mesh_params)
  
mu = p(1);
nu  = p(2);

  % Auxiliary variables
  n = mesh_params.nz;
  u_out = u(1:n);
  


  F = [mesh_params.wz*(u_out.^2)/mesh_params.Lz, u(end)];

end
