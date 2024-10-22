% main routine for secant continuation in ffcore for directional quenching in SH
function [] = main_fcn(mu0,gpuon,pardirec)
%mu initial onset parameter value
%%gpuon=1;     % use gpu (one that is nvidia cuda enabled)
% pardirec %which way to continue in mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  global problem parameter
del = 0.5;     % steepness of parameter ramp (currently should not be used)


mu = mu0;
% output control and cpu/gpu
plotsol=1;  % if plotting during continuation
plotpar=1;  % plot wavenumber continuation curve while computing


verbose=2;  % 1: output after each completed Newton step; 2 output of each Newton step
% load initial data
load_data=0; % load initial condition from file to be specified below
save_on = 1; %save data while continuing
dstep=3; % save every dstep number of continuation steps

size(mu)
%%%Known pinched double root parameters (spreading speed cxlin, decay rate nulin, wavenumber ellxlin) for given initial mu :
mu6 = sqrt(1+6*mu);
mu16 = sqrt(-1+mu6);
cxlin = 4*(2+mu6)*mu16/(3*sqrt(3))

%nulin =
%-( 2*2^(1/3)*(1-sqrt(3)*i)/(-384*sqrt(3)*mu16-192*sqrt(3)*mu6*mu16 + sqrt(442368+(-384*sqrt(3)*mu16-192*sqrt(3)*mu6*mu16)^2  ) )^(1/3)   )...
%    +(1/(24*2^(1/3)))*(1+sqrt(3)*i)*(-384*sqrt(3)*mu16-192*sqrt(3)*mu6*mu16 + sqrt(442368+(-384*sqrt(3)*mu16-5192*sqrt(3)*mu6*mu16)^2  ) )^(1/3)

% (-6i-2*sqrt(3)-3*2^(1/3)*i*(-(2+mu16)*mu16+sqrt(mu6*(3+mu6)))^(2/3) + ...
%         2^(1/3)*sqrt(3)*(-(2+mu6)*mu16 + sqrt(mu6*(3+mu6))))/...
%         (6*2^(2/3)*(-(2+mu6)*mu16+sqrt(mu6*(3+mu6)))^(1/3))
s3 = sqrt(3);
nulin = (6i-2*s3 + 3*2^(1/3)*i*(-2*mu16 - mu6*mu16+sqrt((1+6*mu)*(3+mu6)) )^(2/3) + 2^(1/3)*s3*(-2*mu16-mu6*mu16+ sqrt((1+6*mu)*(3+mu6)) )^(2/3))/...
    (6*2^(2/3)*(-2*mu16-mu6*mu16+ sqrt((1+6*mu)*(3+mu6)) )^(1/3))
omlin = imag(-(1+nulin^2)^2+mu+cxlin*nulin)

ellxlin =  3*(3+mu6)^(3/2)/8/(2+mu6)

%just a random guess, hope it converges
alpha0 = 0.0*4*sqrt(mu)/3/2;
beta0 = 0.0*4*sqrt(mu)/3/2;%1/2;

ar0 = alpha0;
ai0 = 0;
br0 = beta0;
bi0 = 0;

nulin_r = real(nulin); nulin_i = imag(nulin);

for par_direc = pardirec

mod_pars_file=['mod_space_mu' num2str(mu0,'%.2e') '_par_direc_' num2str(par_direc) '.mat']



  %%%%%%%%%%%%  INITIALIZE data and grid, operators  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if load_data % overwrites the choice of elly from above but loads other parameters
      load(data_file); % procures  'uall','ellx','elly','n','Ny','Nx','Ly','Lx','mu','del','cx'; 'uall' has CHI ur added back in, full solution, in matrix form
  else
      %Set parameters and define large grid
      [n,Nx,Ny,Lx,Ly,Xi]=set_domain_size;
  end

  [x,k,dsh,pc,ur0,flag] = init_small_grid(n,mu0);  % defines grid and operators and computes a first solution %%%!!!%% 6/24-left off here:
  if flag ~= 0 fprintf('computation of unit roll solution failed initially'); disp(flag); return; end

  % define large grid and ffcore operators in SH
  [X,Y,Kx,Ky,DSH,PC,EP,CHI,CHIp,DcommSH,DcommSHp,phase] = init_SH_discr(Nx,Ny,Lx,Ly,mu0,del,cxlin,Xi);

  N = Nx*Ny;
  % set Newton parameters (many not used, due to omission of adaptive stepping)
[stepend,steps,ds,ntol,nnitermax,nminstep,ngmrestol,newtonflag,maxstepsize,gminner,maxit,arclength_max,tail_low,tail_high,xrefine,yrefine]=init_Newton;

% figure(13)
% imagesc(real(CHI+CHIp))


% set up initial condition for given mu value
%%if load_data %%not currently working
%  w0=[reshape(uall-CHI.*roll_ext(ur0,n,Nx,Ny,Lx,Ly),Ny*Nx,1)-CHIp.*real(((ar+1i*ai).*X+(br+1i*bi)).*exp((nur+1i*nui)X+Y));ar;ai;br;bi;nur;nui;ellx;cx;mu]; %%fix this
%else
%  w0=w_init(X,Y,del,ur0,n,Nx,Ny,Lx,Ly,CHI,CHIp,alpha,beta,nulin,ellxlin;cxlin;mu);
%end

w0=w_init(X,Y,del,ur0,n,Nx,Ny,Lx,Ly,CHI,CHIp,alpha0,beta0,nulin,ellxlin,cxlin,mu0);


up = reshape(w0(1:N),Ny,Nx)+CHIp.*real((beta0*X+alpha0).*exp(nulin*X/ellxlin+1i*Y));
% figure(23)
% subplot(2,1,1)
% imagesc(up)
% colorbar
% subplot(2,1,2)
% imagesc( reshape(w0(1:N),Ny,Nx))
% colorbar
% break

% initialize secant as stepping in the par_direc of the main parameter elly
sec=[zeros(Nx*Ny+8,1);par_direc]; % negative direction in parameter

%initial continuation parameter grids (possibly complex)
alphagrid = zeros(stepend,1);
betagrid = zeros(stepend,1);
nugrid = zeros(stepend,1);
ellxgrid=zeros(stepend,1);  % ellx results along solutions curves
cxgrid = zeros(stepend,1);
mugrid=zeros(stepend,1); % mu results along solutions curve
if plotsol    h=figure(1); him=figure(3); end
if plotpar    h1=figure(2); end


if gpuon %%%%!!!need to EDITl, not functioning right now!!!
    [CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxitf] = movedatatogpu(CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit);
    gminner=min(gminner,floor(2^31/(Nx*Ny)/2.)); % make sure theres enough memory on device
end


cont_start=tic;
newton_control='none';

%% march one step in the direction of the previous secant
wold=w0;
winit=w0; % only march after step 1

%%%%%%%%%%%%%%%%%% MAIN CONTINUATION LOOP %%%
while (steps<stepend) && w0(end)>1e-2 && w0(end)<0.98 && ds > 1e-12%&&(max(w0(20/32*Nx*Ny:24/32*Nx*Ny),[],'all')>1e-4)


tic
  % define weights and preconditioner for this continuation step
  ellx=w0(end-2);
  EP0=EP(ellx);
    N = Nx*Ny;
  % weight for secant
    sece=sec.*EP0.^(-2);
    % compute residual after marching with added row corresponding to secant condition
    [nrhs,ur0,ur,ur1]=f0(w0,CHI,CHIp,X,Y,Lx,Ly,Nx,Ny,DcommSH,DcommSHp,DSH,phase,ur0,n,x,k,dsh,pc,sece,winit);
  % for linear solver of Newton iteration apply exponential weights
  nrhse = nrhs./EP0;

  %figure(30)
  %imagesc(reshape(nrhse(1:N),Ny,Nx))
  %colorbar
  %drawnow

  % symmetrize: make odd function
%    nrhse = [1/2*(nrhse(1:N) - reshape(circshift(reshape(nrhse(1:N),Ny,Nx),Ny/2,1),Ny*Nx,1));nrhse(end-5:end)];
    nresidual = norm(nrhse);                % estimate residual
    % initialize counter for number of Newton steps
    nniter = 1;
    fprintf(['starting Newton with mu=' num2str(winit(end)) ', '])
    %%%%%%%%%%%%%%%%%%%%%%%% Newton loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cx = w0(end-1);
ellx = w0(end-2);
mu = w0(end);
while (nresidual>ntol)
% prepare the Jacobian by computing parameter derivative only once

dp= 1e-6;



dfar = (f0(w0+dp*double(1:Nx*Ny+9==(N+1))',CHI,CHIp,X,Y,Lx,Ly,Nx,Ny,DcommSH,DcommSHp,DSH,phase,ur0,n,x,k,dsh,pc,sece,winit)-nrhs)/dp;

dfai = (f0(w0+dp*double(1:Nx*Ny+9==(N+2))',CHI,CHIp,X,Y,Lx,Ly,Nx,Ny,DcommSH,DcommSHp,DSH,phase,ur0,n,x,k,dsh,pc,sece,winit)-nrhs)/dp;

dfbr = (f0(w0+dp*double(1:Nx*Ny+9==(N+3))',CHI,CHIp,X,Y,Lx,Ly,Nx,Ny,DcommSH,DcommSHp,DSH,phase,ur0,n,x,k,dsh,pc,sece,winit)-nrhs)/dp;

dfbi = (f0(w0+dp*double(1:Nx*Ny+9==(N+4))',CHI,CHIp,X,Y,Lx,Ly,Nx,Ny,DcommSH,DcommSHp,DSH,phase,ur0,n,x,k,dsh,pc,sece,winit)-nrhs)/dp;

dfnur = (f0(w0+dp*double(1:Nx*Ny+9==(N+5))',CHI,CHIp,X,Y,Lx,Ly,Nx,Ny,DcommSH,DcommSHp,DSH,phase,ur0,n,x,k,dsh,pc,sece,winit)-nrhs)/dp;

dfnui = (f0(w0+dp*double(1:Nx*Ny+9==(N+6))',CHI,CHIp,X,Y,Lx,Ly,Nx,Ny,DcommSH,DcommSHp,DSH,phase,ur0,n,x,k,dsh,pc,sece,winit)-nrhs)/dp;

dfellx = (f0(w0+dp*double(1:Nx*Ny+9==(N+7))',CHI,CHIp,X,Y,Lx,Ly,Nx,Ny,DcommSH,DcommSHp,DSH,phase,ur0,n,x,k,dsh,pc,sece,winit)-nrhs)/dp;

dfcx = (f0(w0+dp*double(1:Nx*Ny+9==(N+8))',CHI,CHIp,X,Y,Lx,Ly,Nx,Ny,DcommSH,DcommSHp,DSH,phase,ur0,n,x,k,dsh,pc,sece,winit)-nrhs)/dp;
;

dfmu = (f0(w0+dp*double(1:Nx*Ny+9==(N+9))',CHI,CHIp,X,Y,Lx,Ly,Nx,Ny,DcommSH,DcommSHp,DSH,phase,ur0,n,x,k,dsh,pc,sece,winit)-nrhs)/dp;
;
% now define function handle to combine ell derivative with on the fly evaluation of jac applied to w; include exponential conjugation
        dfroot = @(dw)  f1(w0,EP0.*dw,CHI,CHIp,X,Y,Nx,Ny,DSH,phase,ur,ur1,dfar,dfai,dfbr,dfbi,dfnur,dfnui,dfellx,dfcx,dfmu,sece)./EP0;

        %%%%%%%%%%%% heres computing the Newton increment; all the work is done here
        if verbose==2 fprintf(['now solving linear system with res ' num2str(nresidual) '...']); end

        [nnincre,flag,relres,iter,resvec]=gmres(dfroot,nrhse,gminner,ngmrestol,maxit,@(w) PC(w,ellx,cx,mu));
    % figure(99)
    % plot(resvec)
    % figure(100)
    % imagesc(reshape(abs(nnincre(1:N)).*EP0(1:N),Ny,Nx))
    % colorbar
    %nnincre(end-5:end)
                % if linear solver did not converge
                if flag>0   disp(['gmres did not converge, residual is ' num2str(nresidual) ' after ' num2str(iter) ' iterations']); newtonflag=2;
                    break
                end
                % actual Newton step
                w0=w0-nnincre.*EP0;



        % symmetrize ro remove round-off errors
       % w0 = [1/2*(w0(1:N) - reshape(circshift(reshape(w0(1:N),Ny,Nx),Ny/2,1),Ny*Nx,1));w0(N+1:end)]; % make odd function
        nniter=nniter+1;
        % recompute residual and estimate norm
        [nrhs,ur0,ur,ur1]=f0(w0,CHI,CHIp,X,Y,Lx,Ly,Nx,Ny,DcommSH,DcommSHp,DSH,phase,ur0,n,x,k,dsh,pc,sece,winit);

        nrhse=nrhs./EP0;
       % nrhse = 1/2*(nrhse - [reshape(circshift(reshape(nrhse(1:N),Ny,Nx),Ny/2,1),Ny*Nx,1);nrhse(N+1:end)]); % make odd function

        nresidual=norm(nrhse);% estimate residual

        %figure(30)
        % imagesc(reshape(abs(nrhse(1:N)),Ny,Nx))
        % colorbar
        % drawnow
        % figure(31)
        % plot(nrhs(N+1:N+9))
        % drawnow
        %
        % figure(31)
        % subplot(2,1,1)
        % plot(log(abs(nrhse)))
        % subplot(2,1,2)
        % plot(log(abs(w0)))


        if verbose>=2
            disp(['nres=' num2str(nresidual)  ', last gmresiters outer=' num2str(iter(1)) ', inner=' num2str(iter(2))])
        end
        %
        if ((nniter>nnitermax)&&steps>1 )|| (nniter>30)
            disp(['Maximal number of iterations reached in large domain problem, giving up; residual is ' num2str(nresidual)])
            newtonflag=1;
            break
        end
        %
        if  norm(nnincre)<nminstep || (abs(ds)<1e-14)
            disp(['Newton step is ineffective, giving up; residual is ' num2str(nresidual)])
            newtonflag=1;
            break
        end
        %
        if  (nresidual>100)&&(steps>1)
            disp(['Residual too large, is ' num2str(nresidual)])
            newtonflag=1;
            break
        end
        %


end
toc

%%%%%%%%%%%%%%% end Newton loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%f
%%%%%%%%%%%%%
cx = w0(end-1); mu = w0(end);

%%%adaptive re-meshing/stepping:
if verbose>=1
         disp(['#Newton iters=' num2str(nniter-1)  ', mu=' num2str(w0(end)) ', ellx=' num2str(w0(end-2)) ', cx=' num2str(w0(end-1)) ', Nx=' num2str(Nx)  ', Ny=' num2str(Ny) ', Lx=' num2str(Lx) ', time=' num2str(toc)]) % output summary from this step
    end
    %
    % adaptive stepsize and gridsize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newton_control='none'; % perform tests on necessary grid changes
    wh=fft2(reshape(w0(1:N),Ny,Nx));
    if (newtonflag==1) || (newtonflag ==2)
        newton_control='reduce_step';
    elseif  (max(abs(wh(7*Ny/16:9*Ny/16,:)),[],'all')>yrefine*sqrt(Nx*Ny)) &&(steps>3)
        newton_control='refine_y';gridflag=1;
    elseif  (max(abs(wh(:,7*Nx/16:9*Nx/16)),[],'all')>xrefine*sqrt(Nx*Ny)) &&(steps>3)
        newton_control='refine_x';gridflag=1;
    elseif  (max(abs(wh(Ny/4:3*Ny/4,:)),[],'all')<yrefine*sqrt(Nx*Ny)/10000) &&(steps>3) && Ny >=128
        newton_control='coarsen_y';gridflag=1;
    elseif  (max(abs(wh(:,Nx/4:3*Nx/4)),[],'all')<xrefine*sqrt(Nx*Ny)/1e5) &&(steps>3) && Nx >512*2
        newton_control='coarsen_x';gridflag=1;
    elseif (max(w0(1/16*Nx*Ny:6/16*Nx*Ny),[],'all')<tail_low) &&(steps>3)&&(Lx/w0(N)>1000)  % ((max(w0(N-3*Nx/8*Ny:N-Nx/8*Ny),[],'all')<tail_low) && if also check upstream
        newton_control='cut_x';
    elseif (max(w0(1/16*Nx*Ny:3/16*Nx*Ny),[],'all')>tail_high) &&(steps>3)  %   ((max(w0(N-3*Nx/8*Ny:N-Nx/8*Ny),[],'all')>tail_high) || if also check upstream
        newton_control='ext_x';
    end
    % execute necessary changes
    switch newton_control
        case 'reduce_step'
            ds=ds/2;
            disp(['Reduce step size to  ' num2str(ds)])
            newtonflag=0;
            winit=wold+(steps>2)*ds*sec+(steps==2)*1e-4*ds*sec; % only march after step 2
            w0=winit;
            %figure(40)
            %plot(abs(sec))
        case 'refine_y'
            [w0,winit,wold,sec,Nx,Ny]=grid_refine(w0,winit,wold,sec,Nx,Ny,'y');
            if verbose>= 1 disp(['refined y-grid, new Ny=' num2str(Ny) ', new dy=' num2str(Ly/Ny) ', dx='  num2str(Lx/w0(end-1)/Nx) ]); end
            [X,Y,Kx,Ky,DSH,PC,EP,CHI,CHIp,DcommSH,DcommSHp,phase]=init_SH_discr(Nx,Ny,Lx,Ly,mu,del,cx,Xi);
                        size(EP(1))
            if gpuon [CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit] = movedatatogpu(CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit); end
            % adjust tolerances
            ntol=ntol*sqrt(2); % tolerance fort the Newton solver, max of residual
            nminstep=nminstep*sqrt(2); % minimal Newton step, if making smaller steps without bringing the residual down give up
            ngmrestol=ngmrestol*sqrt(2); % tolerance for the gmres solver; shjould be smaller than the tolerance for the Newton solver
            if gpuon
                gminner=min(gminner,floor(2^31/(Nx*Ny)/2.)); % make sure there's enough memory on devise
            end
        case 'refine_x'
            [w0,winit,wold,sec,Nx,Ny]=grid_refine(w0,winit,wold,sec,Nx,Ny,'x');
            if verbose>= 1 disp(['refined x-grid, new Nx=' num2str(Nx) ', new dy=' num2str(Ly/Ny) ', dx='  num2str(Lx/w0(end-3)/Nx) ]); end
            [X,Y,Kx,Ky,DSH,PC,EP,CHI,CHIp,DcommSH,DcommSHp,phase]=init_SH_discr(Nx,Ny,Lx,Ly,mu,del,cx,Xi);

            if gpuon  [CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit] = movedatatogpu(CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit); end
            % adjust tolerances
            ntol=ntol*sqrt(2);              % tolerance fort the Newton solver, max of residual
            nminstep=nminstep*sqrt(2);      % minimal Newton step, if making smaller steps without bringing the residual down give up
            ngmrestol=ngmrestol*sqrt(2);    % tolerance for the gmres solver; shjould be smaller than the tolerance for the Newton solver
            if gpuon
                gminner=min(gminner,floor(2^31/(Nx*Ny)/2.)); % make sure there's enough memory on devise
            end
             ds=1e-8*ds;
        case 'coarsen_x'
            [w0,winit,wold,sec,Nx,Ny]=grid_coarsen(w0,winit,wold,sec,Nx,Ny,'x');
            if verbose>= 1 disp(['coarsened x-grid, new Nx=' num2str(Nx) ', new dy=' num2str(Ly/Ny) ', dx='  num2str(Lx/w0(end-3)/Nx) ]); end
            [X,Y,Kx,Ky,DSH,PC,EP,CHI,CHIp,DcommSH,DcommSHp,phase]=init_SH_discr(Nx,Ny,Lx,Ly,mu,del,cx,Xi);
            if gpuon  [CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit] = movedatatogpu(CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit); end
            % adjust tolerances
            ntol=ntol/sqrt(2);              % tolerance fort the Newton solver, max of residual
            nminstep=nminstep/sqrt(2);      % minimal Newton step, if making smaller steps without bringing the residual down give up
            ngmrestol=ngmrestol/sqrt(2);    % tolerance for the gmres solver; shjould be smaller than the tolerance for the Newton solver
            if gpuon
                gminner=min(gminner,floor(2^31/(Nx*Ny)/2.)); % make sure there's enough memory on devise
            end
        case 'coarsen_y'
            [w0,winit,wold,sec,Nx,Ny]=grid_coarsen(w0,winit,wold,sec,Nx,Ny,'y');
            if verbose>= 1 disp(['coarsened y-grid, new Nx=' num2str(Nx) ', new dy=' num2str(Ly/Ny) ', dx='  num2str(Lx/w0(end-3)/Nx) ]); end
              [X,Y,Kx,Ky,DSH,PC,EP,CHI,CHIp,DcommSH,DcommSHp,phase]=init_SH_discr(Nx,Ny,Lx,Ly,mu,del,cx,Xi);
            if gpuon  [CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit] = movedatatogpu(CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit);  end
            % adjust tolerances
            ntol=ntol/sqrt(2);              % tolerance fort the Newton solver, max of residual
            nminstep=nminstep/sqrt(2);      % minimal Newton step, if making smaller steps without bringing the residual down give up
            ngmrestol=ngmrestol/sqrt(2);    % tolerance for the gmres solver; shjould be smaller than the tolerance for the Newton solver
            if gpuon
                gminner=min(gminner,floor(2^31/(Nx*Ny)/2.)); % make sure there's enough memory on devise
            end
        case 'cut_x'
            [w0,winit,wold,sec,Nx,Lx,Xi]=grid_cut(w0,winit,wold,sec,Nx,Ny,Lx,Xi);
            if verbose>= 1 disp(['cut x-domain, new Lx=' num2str(Lx)]); end
              [X,Y,Kx,Ky,DSH,PC,EP,CHI,CHIp,DcommSH,DcommSHp,phase]=init_SH_discr(Nx,Ny,Lx,Ly,mu,del,cx,Xi);
            if gpuon  [CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit] = movedatatogpu(CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit);  end
            % adjust tolerances
            ntol=ntol/sqrt(2); % tolerance fort the Newton solver, max of residual
            nminstep=nminstep/sqrt(2); % minimal Newton step, if making smaller steps without bringing the residual down give up
            ngmrestol=ngmrestol/sqrt(2); % tolerance for the gmres solver; shjould be smaller than the tolerance for the Newton solver
            if gpuon
                gminner=min(gminner,floor(2^31/(Nx*Ny)/2.)); % make sure there's enough memory on devise
            end
        case 'ext_x'
            [w0,winit,wold,sec,Nx,Lx,Xi]=grid_extend(w0,winit,wold,sec,Nx,Ny,Lx,Xi);
            if verbose>= 1 disp(['extended x-grid, new Lx=' num2str(Lx)]); end
              [X,Y,Kx,Ky,DSH,PC,EP,CHI,CHIp,DcommSH,DcommSHp,phase]=init_SH_discr(Nx,Ny,Lx,Ly,mu,del,cx,Xi);
            if gpuon  [CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit] = movedatatogpu(CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit);  end
            % adjust tolerances
            ntol=ntol*sqrt(2) % tolerance fort the Newton solver, max of residual
            nminstep=nminstep*sqrt(2); % minimal Newton step, if making smaller steps without bringing the residual down give up
            ngmrestol=ngmrestol*sqrt(2); % tolerance for the gmres solver; shjould be smaller than the tolerance for the Newton solver
            if gpuon
                gminner=min(gminner,floor(2^31/(Nx*Ny)/2.)); % make sure there's enough memory on devise
            end
            ds=1e-1*ds; % reduce step size since secant will be slightly wrong as wold will not be recomputed;
        otherwise % this is where everything went fine, so collect info and move on to next step
        alphagrid(steps) = w0(N+1)+1i*w0(N+2);
        betagrid(steps) = w0(N+3)+1i*w0(N+4);
        nugrid(steps) = w0(N+5)+1i*w0(N+6);
        ellxgrid(steps) = w0(N+7);
        cxgrid(steps) = w0(N+8);
        mugrid(steps) = w0(N+9);
        par_grids = [alphagrid,betagrid,nugrid,ellxgrid,cxgrid,mugrid];
        if plotsol ||plotpar plotting(plotsol,plotpar,h,h1,him,w0,Nx,Ny,CHI,CHIp,ur,Lx,Ly,X,Y,par_grids,steps,tail_high); end
        % compute new secant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if steps>0 sec=w0-wold;sec=sec/norm(sec); end  % update secant (only after initial step which is Newton without marching)
        if (nniter <4)&&(steps>2)&&(ds*norm(sec(end-2:end))/min(abs(w0(end-2:end)))<arclength_max)&&(ds<maxstepsize*sqrt(Nx*Ny))% Newton too short, increase ds for next step
            ds=min(ds*1.4); if verbose >=1 disp(['Increase step size to  ' num2str(ds)]); end
        end
        steps=steps+1;
        % save continbuation progress in wavenumber space'
        % *********************
        while (ds*norm(sec(end-1:end))/min(abs(w0(end-1:end)))>2*arclength_max)
            ds=ds/2.1;
        end
        save(mod_pars_file,'par_grids','steps'); % display with plot(ellygrid(1:steps),ellxgrid(1:steps),'o-')
        if save_on && (mod(steps,dstep)==0||steps==stepend)
            muold = mu;
            mu=w0(N+9);
            cx=w0(N+8);
            ellx = w0(N+7);
            nu = w0(N+5)+1i*w0(N+6);
            beta = w0(N+3)+1i*w0(N+4);
            alpha = w0(N+1) + 1i*w0(N+2);

            file_name=['data/data_mu_' num2str(mu,'%.2e') '.mat'];
            uall=reshape(w0(1:N),Ny,Nx)+CHI.*ur+CHIp.*real((alpha*X+beta).*exp(nu*X/ellx + 1i*Y));
            save(file_name,'uall','alpha','beta','nu','mu','ellx','cx','n','Ny','Nx','Ly','Lx','mu','del','Xi','ur');
            %clear uall; % don't need this in the code other than for record keeping
          %  cx = cxold;
        end
        % march one step in the direction of the previous secant
        wold=w0;
        %(sign(par_direc*sec(end)))
        %steps
        winit=w0+(steps>3)*ds*sec+(steps==2)*1e-1*ds*sec+(steps==3)*1e-1*ds*sec*sign(sec(end)*par_direc); % only march after step 1; in step 2 reduce marching since only in direction of parameter
        w0=winit;
        %sec(end)
        end % END of switch controls

end %%%%END OF CONTINUATION LOOP



























end
