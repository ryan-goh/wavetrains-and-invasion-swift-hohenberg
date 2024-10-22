%%%%profile-plotter
clear all
close all
%plot with cx fixed, or  ky fixed,
fix_cx = 1;

%%pick a moduli slice, and plot data points along
slice_name = 'mod_space_mu4.00e-01_par_direc_1.mat';
slice_name_down = 'mod_space_mu4.00e-01_par_direc_-1.mat';
data_names = {'data/data_mu_1.11e-01.mat','data/data_mu_5.05e-01.mat','data/data_mu_8.57e-01.mat'};
num_datapoint =length(data_names);
title_pl = 1;


  load(slice_name);
  ELX = gather(par_grids(:,4));
  ELX = ELX(1:steps-1);
  CX = gather(par_grids(:,5));
  CX = CX(1:steps-1);
  MU = gather(par_grids(:,6));
  MU = MU(1:steps-1);

  figure(3)
  plot(MU,CX)

  clear par_grids
    load(slice_name_down);
  ELX = [gather(par_grids(steps-1:-1:1,4));ELX];
  CX =  [gather(par_grids(steps-1:-1:1,5));CX];
  MU =  [gather(par_grids(steps-1:-1:1,6));MU];
    
  figure(3)
   plot(MU,CX)


  figure(1)%plot the slice of the moduli space
  plot(MU,CX,'LineWidth',2.5)
  xlabel('$\mu$','Interpreter','Latex')
  ylabel('$c_\mathrm{lin}$','Interpreter','Latex')
  ax = gca;
  ax.FontSize = 18;
  ax.TickLabelInterpreter = 'latex';



  figure(2)%plot the slice of the moduli space
  plot(MU,ELX,'LineWidth',2.5)
  xlabel('$\mu$','Interpreter','Latex')
  ylabel('$k_\mathrm{lin}$','Interpreter','Latex')
  ax = gca;
  ax.FontSize = 18;
  ax.TickLabelInterpreter = 'latex';


  %%add theoretical linear predictions into figures
  %%%Known pinched double root parameters (spreading speed cxlin, decay rate nulin, wavenumber ellxlin) for given initial mu :
  mu = [0:0.005:1]';
mu6 = sqrt(1+6*mu);
mu16 = sqrt(-1+mu6);
cxlin = 4*(2+mu6).*mu16/(3*sqrt(3));

%nulin = 
%-( 2*2^(1/3)*(1-sqrt(3)*i)/(-384*sqrt(3)*mu16-192*sqrt(3)*mu6*mu16 + sqrt(442368+(-384*sqrt(3)*mu16-192*sqrt(3)*mu6*mu16)^2  ) )^(1/3)   )...
%    +(1/(24*2^(1/3)))*(1+sqrt(3)*i)*(-384*sqrt(3)*mu16-192*sqrt(3)*mu6*mu16 + sqrt(442368+(-384*sqrt(3)*mu16-5192*sqrt(3)*mu6*mu16)^2  ) )^(1/3)  

% (-6i-2*sqrt(3)-3*2^(1/3)*i*(-(2+mu16)*mu16+sqrt(mu6*(3+mu6)))^(2/3) + ...
%         2^(1/3)*sqrt(3)*(-(2+mu6)*mu16 + sqrt(mu6*(3+mu6))))/...
%         (6*2^(2/3)*(-(2+mu6)*mu16+sqrt(mu6*(3+mu6)))^(1/3))
s3 = sqrt(3);
nulin = (6i-2*s3 + 3*2^(1/3)*i*(-2*mu16 - mu6.*mu16+sqrt((1+6*mu).*(3+mu6)) ).^(2/3) + 2^(1/3)*s3.*(-2*mu16-mu6.*mu16+ sqrt((1+6*mu).*(3+mu6)) ).^(2/3))./...
    (6*2^(2/3)*(-2*mu16-mu6.*mu16+ sqrt((1+6*mu).*(3+mu6)) ).^(1/3));
omlin = imag(-(1+nulin.^2).^2+mu+cxlin.*nulin);

ellxlin =  3*(3+mu6).^(3/2)./8./(2+mu6);

figure(1)
hold on 
plot(mu,cxlin,'k--','LineWidth',2)
hold off


figure(2)
hold on
plot(mu,ellxlin,'k--','LineWidth',2)
hold off

clear mu nulin cxlin omlin ellxlin;

him = figure(3)


  for ii = 1:num_datapoint
    if num_datapoint ==1
      load(data_names{ii})
      %form meshgrid
      X0 = linspace(-Lx,Lx,Nx+1);X0=X0(1:end-1)';
      X0=X0-Xi;
      Y0 = linspace(-Ly,Ly,Ny+1);Y0=Y0(1:end-1)';
      [X, Y] = meshgrid(X0,Y0);

      
     

      figure(11)
      hold on
      surf(X,Y,uall)
      shading interp
      view([0 90])
      daspect([3 1 1])
      xlim([min(X0) max(X0)])
      ylim([min(Y0) max(Y0)])
      xlabel('$\zeta$','Interpreter','Latex')
      ylabel('$\tau$','Interpreter','Latex')
      yticks([-pi 0 pi-3*Ly/Ny ])
      yticklabels({'$0$','$\pi$','$2\pi$'})
      %xticks([ceil(-Lx-Xi) 0 floor(Lx-Xi)])
      ax = gca;
      ax.FontSize = 18;
      ax.TickLabelInterpreter = 'latex';
      if title_pl
        title(['$(\mu,c_x,k_x) = ( $' num2str(mu) ',' num2str(cx) ',' num2str(ellx) ')'],'Interpreter','latex')
      end
      hold off

      figure(2)
      hold on
      plot(mu,ellx,'o','LineWidth',2)
      %title(['$c_x = ' num2str(cx)],'Interpreter','latex')
      hold off

    else
      load(data_names{ii})
      %form meshgrid
      X0 = linspace(-Lx,Lx,Nx+1);X0=X0(1:end-1)';
      X0=X0-Xi;
      Y0 = linspace(-Ly,Ly,Ny+1);Y0=Y0(1:end-1)';
      [X, Y] = meshgrid(X0,Y0);

      fig11 = figure(11)
      hold on
      subplot(num_datapoint,1,ii)
      surf(X,Y,uall)
      shading interp
      view([0 90])
      %daspect([3 1 1])
      xlim([min(X0) max(X0)])
      ylim([min(Y0) max(Y0)])
      xlabel('$\zeta$','Interpreter','Latex')
      ylabel('$\tau$','Interpreter','Latex')
      yticks([-pi 0 pi-3*Ly/Ny ])
      yticklabels({'$0$','$\pi$','$2\pi$'})
      %xticks([ceil(-Lx-Xi) 0 floor(Lx-Xi)])
      ax = gca;
      ax.FontSize = 18;
      colorbar
      if ii == 1
        caxis([-0.5 0.5])
      else
           caxis([-1.1 1.1])
      end
      ax.TickLabelInterpreter = 'latex';
      if title_pl
          title(['$( $' num2str(ii) '$ ) $'],'Interpreter','latex')
       %title(['$(\mu,c_x,k_x) = ( $' num2str(mu) ',' num2str(cx) ',' num2str(ellx) ')'],'Interpreter','latex')
      end
      hold off
        


      fig1 = figure(1)
      hold on
      plot(mu,cx,'o','LineWidth',2)
      text(mu+0.02,cx,['(' num2str(ii) ')'],'Interpreter','latex','FontSize',20)
         xlim([0.00 0.98])
      % xlim([0.12 0.9])
     

      fig2 = figure(2)
      hold on
      plot(mu,ellx,'o','LineWidth',2)
      text(mu+0.02,ellx,['(' num2str(ii) ')'],'Interpreter','latex','FontSize',20)
      %title(['$c_x = $' num2str(cx)],'Interpreter','latex')
       xlim([0.00 0.98])
      hold off
        
    % if ii ==2
    % 
    % % cutoff supported towards +L and commutators
    %     eps  = -1.5;-0.35;  X0=-20;
    %     EX   =  exp(min(-(X-X0)/eps,500));
    %     %  EX   =  exp(-(X-X0)/eps);
    %     CHI  =  1./(1+EX);
    %     CHIm=  1-CHI; %CHIp = CHIm;
    %     EXp   =  exp(min((X-1)/eps,500));
    %     CHIp  =  1./(1+EXp);
    %     CHIpm = 1-CHIp;
    % 
    %     w0 = uall-(CHI.*ur+CHIp.*real((alpha*X+beta).*exp(nu*X/ellx + 1i*Y)));
    %     set(0,'CurrentFigure',him);
    %     subplot(3,1,1)
    %     imagesc(X(1:Ny:N),Y(1:Nx:N),reshape(w0(1:N),Ny,Nx)+CHI.*ur+CHIp.*real((alpha*X+beta).*exp(nu*X/ellx+1i*Y)));
    %     colorbar
    %      subplot(3,1,2)
    %      imagesc(X(1:Ny:N),Y(1:Nx:N),reshape(w0(1:N),Ny,Nx));
    %      colorbar
    %       subplot(3,1,3)
    %       imagesc(X(1:Ny:N),Y(1:Nx:N),CHIp.*real((alpha*X/ellx+beta).*exp(nu*X+1i*Y)));
    % end

    end

  end

  figure(1)
  legend({'AUTO Numerics','Linear Prediction'},'Location','southeast','Interpreter','latex')

  figure(2)
   legend({'AUTO Numerics','Linear Prediction'},'Location','southeast','Interpreter','latex')


         saveas(fig11,'profile_mu.eps','epsc')

          saveas(fig1,'speed_mu.eps','epsc')

           saveas(fig2,'wavenumber_mu.eps','epsc')

%%%%
