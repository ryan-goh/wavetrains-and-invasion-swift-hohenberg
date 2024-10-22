function plotting(plotsol,plotpar,h,h1,him,w0,Nx,Ny,CHI,CHIp,ur,Lx,Ly,X,Y,par_grids,steps,tail_high);
  %%   par_grids = [alphagrid,betagrid,nugrid,ellxgrid,cxgrid,mugrid];

  cxgrid = par_grids(:,5);
  mugrid = par_grids(:,6);
  ellxgrid = par_grids(:,4);
  alphagrid = par_grids(:,1);
  betagrid = par_grids(:,2);
  N = Nx*Ny;

    if plotsol

        alpha = par_grids(steps,1);
        beta = par_grids(steps,2);
        nu = par_grids(steps,3);
        ellx=par_grids(steps,4);
        cx =par_grids(steps,5);
        mu=par_grids(steps,6);
        set(0,'CurrentFigure',him);
        subplot(3,1,1)
       % imagesc(X(1:Ny:N),Y(1:Nx:N),reshape(w0(1:N),Ny,Nx)+CHI.*ur+CHIp.*real((alpha*X+beta).*exp(nu*X/ellx+1i*Y)));
       surf(X,Y,reshape(w0(1:N),Ny,Nx)+CHI.*ur+CHIp.*real((alpha*X+beta).*exp(nu*X/ellx+1i*Y)));
       shading interp
        view([0 90])
         xlim([X(1) X(end)])
         ylim([Y(1) Y(end)])
            %   xlabel('$\zeta$','Interpreter','Latex','FontSize',18)
         ylabel('$\tau$','Interpreter','Latex','FontSize',18)
        yticks([-pi 0 pi-3*Ly/Ny ])
        yticklabels({'$0$','$\pi$','$2\pi$'})
              ax = gca;
      ax.FontSize = 16;
      ax.TickLabelInterpreter = 'latex';
        %xticks([ceil(-Lx-Xi) 0 floor(Lx-Xi)])
          % xlabel('$x$','Interpreter','Latex')
          % ylabel('$y$','Interpreter','Latex')
          % yticks([-pi 0 pi-3*Ly/Ny ])
          % yticklabels({'$0$','$\pi$','$2\pi$'})
          % xticks([ceil(-Lx-Xi) 0 floor(Lx-Xi)])
          title('$u_f$','Interpreter','latex','FontSize',20)
        colorbar
         subplot(3,1,2)
         %imagesc(X(1:Ny:N),Y(1:Nx:N),reshape(w0(1:N),Ny,Nx));
         surf(X,Y,reshape(w0(1:N),Ny,Nx));
                shading interp
        view([0 90])
         xlim([X(1) X(end)])
         ylim([Y(1) Y(end)])
        %  xlabel('$\zeta$','Interpreter','Latex','FontSize',18)
         ylabel('$\tau$','Interpreter','Latex','FontSize',18)
        yticks([-pi 0 pi-3*Ly/Ny ])
        yticklabels({'$0$','$\pi$','$2\pi$'})
                      ax = gca;
      ax.FontSize = 16;
      ax.TickLabelInterpreter = 'latex';
        %xticks([ceil(-Lx-Xi) 0 floor(Lx-Xi)])
          %      xlim([min(X0) max(X0)])
          % ylim([min(Y0) max(Y0)])
          % xlabel('$x$','Interpreter','Latex')
          % ylabel('$y$','Interpreter','Latex')
          % yticks([-pi 0 pi-3*Ly/Ny ])
          % yticklabels({'$0$','$\pi$','$2\pi$'})
          % xticks([ceil(-Lx-Xi) 0 floor(Lx-Xi)])
          title('$w$','Interpreter','latex','FontSize',20)
         colorbar
          subplot(3,1,3)
          %imagesc(X(1:Ny:N),Y(1:Nx:N),CHIp.*real((alpha*X+beta).*exp(nu*X/ellx+1i*Y)));
                   surf(X,Y,CHIp.*real((alpha*X+beta).*exp(nu*X/ellx+1i*Y)));
                shading interp
        view([0 90])
         xlim([X(1) X(end)])
         ylim([Y(1) Y(end)])
         xlabel('$\zeta$','Interpreter','Latex','FontSize',18)
         ylabel('$\tau$','Interpreter','Latex','FontSize',18)
        yticks([-pi 0 pi-3*Ly/Ny ])
        yticklabels({'$0$','$\pi$','$2\pi$'})
                      ax = gca;
      ax.FontSize = 16;
      ax.TickLabelInterpreter = 'latex';
        %xticks([ceil(-Lx-Xi) 0 floor(Lx-Xi)])
          % xlim([min(X0) max(X0)])
          % ylim([min(Y0) max(Y0)])
          % xlabel('$x$','Interpreter','Latex')
          % ylabel('$y$','Interpreter','Latex')
          % yticks([-pi 0 pi-3*Ly/Ny ])
          % yticklabels({'$0$','$\pi$','$2\pi$'})
          % xticks([ceil(-Lx-Xi) 0 floor(Lx-Xi)])
          title('$\chi_+u_{\nu,\alpha,\beta}$','Interpreter','latex','FontSize',20)
         colorbar
          % y-length vs number of data points
      %  daspect([gather(Ly/elly/Ny),gather(Lx/ellx/Nx),1]);
        drawnow;%shg
        set(0,'CurrentFigure',h);
        plot(X(1:Ny:N),w0(1:Ny:N));
        %axis([gather(X(1,1)) gather(X(1,end)) -tail_high tail_high])
        drawnow
     end
     if plotpar
         set(0,'CurrentFigure',h1);
         subplot(3,1,1)
         plot(mugrid(1:steps),ellxgrid(1:steps),'b-o');
         xlabel(['mu'])
         ylabel(['kx'])
         subplot(3,1,2)
         plot(mugrid(1:steps),cxgrid(1:steps),'b-o');
         xlabel(['mu'])
         ylabel(['cx'])
         subplot(3,1,3)
         plot(mugrid(1:steps),[real(alphagrid(1:steps)),imag(alphagrid(1:steps)),...
             real(betagrid(1:steps)),imag(betagrid(1:steps))],'-o');
         xlabel(['mu'])
         ylabel(['alpha,beta'])
         %title(['c_x = ', num2str(cx), ', beta = ', num2str(beta)])
     end

 end
