function w= w_init(X,Y,del,ur0,n,Nx,Ny,Lx,Ly,CHI,CHIp,alpha,beta,nu,ellx,cx,mu)
% generates initial guess by extending roll solution, modulating with the parameter ramp, subtracting the cutoff from ff-core
%    alpha = ar+1i*ai; beta = br+1i*bi; nu = nur+1i*nui;
%    w = ((1+tanh(-(X)/del))/2-CHI).*roll_ext(ur0,n,Nx,Ny,Lx,Ly)-0*CHIp.*real((alpha*X+beta).*exp(nu*X/ellx+1i*Y) );%...
        %-CHIp.*real((alpha*X+beta).*exp(nu*X/ellx+Y) );
    w = (1-CHI).*(1-CHIp).*roll_ext(ur0,n,Nx,Ny,Lx,Ly)-0*CHIp.*real((alpha*X/ellx+beta).*exp(nu*X/ellx+1i*Y) );%...
        %-CHIp.*real((alpha*X+beta).*exp(nu*X/ellx+Y) );
ar = real(alpha);
ai = imag(alpha);
br = real(beta);
bi = imag(beta);
nur = real(nu);
nui = imag(nu)
% figure(29)
% imagesc(w)
% colorbar
% drawnow

%      w=zeros(N,1);%if one wanted to do it with purely real variables
%    ar = real(alpha);
%    ai = imag(alpha);
%    br = real(beta);
%    bi = imag(beta);
%    nur = real(nu);
%    nui = imag(nu);
    w= reshape(w,Nx*Ny,1); % flatten into clumn vector
    w=[w;ar;0;br;0;nur;nui;ellx;cx;mu]; % append the guess for the wavenumber
end
