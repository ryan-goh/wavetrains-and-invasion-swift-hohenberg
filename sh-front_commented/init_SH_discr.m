function  [X,Y,Kx,Ky,DSH,PC,EP,CHI,CHIp,DcommSH,DcommSHp,phase] =init_SH_discr(Nx,Ny,Lx,Ly,mu,del,cx,Xi);

    % large grid mesh
    N = Nx*Ny;
    X0 = linspace(-Lx,Lx,Nx+1);X0=X0(1:end-1)';
    X0=X0-Xi;
    Y0 = linspace(-Ly,Ly,Ny+1);Y0=Y0(1:end-1)';
    [X, Y] = meshgrid(X0,Y0);
    % derivative vectors (in Fourier space)
    Kx   = ([[0:Nx/2] [-Nx/2+1:-1]])'*(pi/Lx);
    Ky   = ([[0:Ny/2] [-Ny/2+1:-1]])'*(pi/Ly);
    Kx=ones(Ny,1)*Kx';
    Ky=Ky*ones(1,Nx);
    % main differential operator , here without the mu part since x dependent
    DSH =@(ellx,cx,mu) -(1-ellx^2*Kx.^2).^2+cx*ellx*i*(Kx-Ky)+mu*ones(Ny,Nx);
    %preconditioner (check on number of additional parameters)
    PC0  = @(du,ellx,cx,mu) [reshape(ifft2(fft2(reshape(du(1:N),Ny,Nx))./(DSH(ellx,cx,mu)-(1.5+mu)*ones(Ny,Nx)),'symmetric'),Ny*Nx,1);du(N+1:end)];

    %%% exponential preconditioner; also need grid preconditioner to weight ellx; not smooth at L
    %  EP = @(ellx) [1/sqrt(Lx*Ly/Nx/Ny)*ones(N,1);1;1].^(-1);        % only grid weight
%      EP=@(ellx) ones(Ny*Nx+2,1);                            % no weight
    %  EP = [reshape(1/sqrt(Lx*Ly/Nx/Ny)*cosh(sqrt(mu)/4*X).^(-1),Ny*Nx,1);.01;.01].^(-1);   % exponential as predicted from analysis in the middle of the spatial spectral gap
   %  EP = @(ellx) [reshape(1/sqrt(Lx*Ly/Nx/Ny)*cosh(min(3/ellx,14)/Lx*X).^(-1),Ny*Nx,1);.01;.01].^(-1); %
%      EP=@(ellx) [reshape(sqrt(Lx*Ly/Nx/Ny)*cosh(1/Lx*X),Ny*Nx,1);100;100]; % prefactor 3 can be ramped up to 14 but no more due to overflow
    %EP=@(ellx) [reshape(sqrt(Lx*Ly/Nx/Ny)*(exp(3*(X-X(end))/Lx)+exp(-3*(X-X(1))/Lx)),Ny*Nx,1);10;10]; % prefactor in exponent can be ramped up to 14 but no more due to overflow; this is asymmetric choice of weights for asymmetric X, Xi>0

    eta = 0.;
    EP=@(ellx) [reshape(sqrt(Lx*Ly/Nx/Ny)*(exp(eta*(X-X(end))/Lx)+exp(-eta*(X-X(1))/Lx)),Ny*Nx,1);1;1;1;1;1;1;1e0;1e0;1e-1];
  %EP=@(ellx) [reshape(sqrt(Nx*Ny/Lx/Ly)*(exp(eta*(X-X(end))/Lx)+exp(-eta*(X-X(1))/Lx)),Ny*Nx,1);1;1;1;1;1;1];

    % prefactor in exponent can be ramped up to 14 but no more due to overflow; this is asymmetric choice of weights for asymmetric X, Xi>0
    % preconditioner
    PC = @(w,ellx,cx,mu) PC0(w.*EP(ellx),ellx,cx,mu)./EP(ellx);



    % cutoff supported towards +L and commutators
    eps  = -1.5;-0.35;  X0=-40;
    EX   =  exp(min(-(X-X0)/eps,500));
    %  EX   =  exp(-(X-X0)/eps);
    CHI  =  1./(1+EX);
    CHIm=  1-CHI; %CHIp = CHIm;
    EXp   =  exp(min((X-1)/eps,500));
    CHIp  =  1./(1+EXp);
    CHIpm = 1-CHIp;
    
    figure(111)
    surf(X,Y,CHI+CHIp)
    shading interp
    view([0 90])

    CHI1 =  CHI.*CHIm/eps;
    CHI2 =  (2*CHIm.^2-CHIm).*CHI/eps^2;
    CHI3 =  (6*CHIm.^3-6*CHIm.^2+CHIm).*CHI/eps^3;
    CHI4 =  (24*CHIm.^4-36*CHIm.^3+14*CHIm.^2-CHIm).*CHI/eps^4;

    CHI1p =  -CHIp.*CHIpm/eps;
    CHI2p =  (2*CHIpm.^2-CHIpm).*CHIp/eps^2;
    CHI3p =  -(6*CHIpm.^3-6*CHIpm.^2+CHIpm).*CHIp/eps^3;
    CHI4p =  (24*CHIpm.^4-36*CHIpm.^3+14*CHIpm.^2-CHIpm).*CHIp/eps^4;

%    % DcommSH = @(u,ellx,cx) -((ellx^4*CHI4+2*ellx^2*CHI2).*u+...
%                4*(ellx^3*CHI3 +ellx*CHI1).*ifft2(i*ellx*Kx.*fft2(u),'symmetric')+...
%                6*ellx^2*CHI2.*ifft2(-ellx^2*Kx.^2.*fft2(u),'symmetric')+...
%                4*ellx*CHI1.*ifft2(-i*ellx^3*Kx.^3.*fft2(u),'symmetric')+...
%                2*ellx^2*CHI2.*ifft2(-elly^2*Ky.^2.*fft2(u),'symmetric')+...
%                4*ellx*CHI1.*ifft2(-i*ellx*elly^2*Ky.^2.*Kx.*fft2(u),'symmetric'))+...
%                cx*ellx*CHI1.*u;

    DcommSH = @(u,ellx,cx)  -((ellx^4*CHI4+2*ellx^2*CHI2).*u+...
            4*(ellx^3*CHI3 +ellx*CHI1).*ifft2(i*ellx*Kx.*fft2(u),'symmetric')+...
            6*ellx^2*CHI2.*ifft2(-ellx^2*Kx.^2.*fft2(u),'symmetric')+...
            4*ellx*CHI1.*ifft2(-1i*ellx^3*Kx.^3.*fft2(u),'symmetric'))+...
            cx*ellx*CHI1.*u;

    DcommSHp = @(u,ellx,cx)  -((ellx^4*CHI4p+2*ellx^2*CHI2p).*u+...
            4*(ellx^3*CHI3p +ellx*CHI1p).*ifft2(1i*ellx*Kx.*fft2(u),'symmetric')+...
            6*ellx^2*CHI2p.*ifft2(-ellx^2*Kx.^2.*fft2(u),'symmetric')+...
            4*ellx*CHI1p.*ifft2(-1i*ellx^3*Kx.^3.*fft2(u),'symmetric'))+...
            cx*ellx*CHI1p.*u;

    phase = (X>-(15+2*pi)).*(X<-15);

end









