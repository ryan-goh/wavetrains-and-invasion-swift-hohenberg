function plotHandle = PlotSolution_cubic(u,p,parentHandle,mesh_params)

   if isempty(parentHandle)
     scrsz = get(0,'ScreenSize');
     plotHandle = figure('Position',[2/4*scrsz(3) scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     parentHandle = plotHandle;
   else
     plotHandle = parentHandle;
   end 

   figure(parentHandle);
   
     % Rename parameters
mu = p(1);
k  = p(2);

  
  % Auxiliary variables
  n = mesh_params.nz;
  u_out     = u(1:n);

  plot(mesh_params.z,u_out,'b');xlabel('x');ylabel('u');
drawnow;


   print -dtiff state.tiff

end
