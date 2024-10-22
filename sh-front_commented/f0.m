%evaluate SH with phase condition, as well as additional equation

function [f,u,ur,ur1] = f0(w,CHI,CHIp,X,Y,Lx,Ly,Nx,Ny,DcommSH,DcommSHp,DSH,phase,ur0,n,x,k,dsh,pc,sece,winit)
N=Nx*Ny;
ar = w(N+1); ai = w(N+2); br=w(N+3); bi = w(N+4); nur = w(N+5); nui = w(N+6);
 ellx = w(N+7); cx = w(N+8);mu = w(N+9);

 alpha = ar+1i*ai;
 beta = br+1i*bi;
 nu = nur+1i*nui;
w=reshape(w(1:N),Ny,Nx); % rearrange as matrix
[u,flag]=unit_roll(ellx,mu,ur0,x,k,dsh,pc);
ur=roll_ext(u,n,Nx,Ny,Lx,Ly);
ueps=ifft(exp(i*k*1e-3).*fft(u),'symmetric');
ur1=(roll_ext(ueps,n,Nx,Ny,Lx,Ly)-ur)/1e-3;
%unu = (1./(1+exp(-5*(X+1)))).*real((alpha*X+beta).*exp(nu*X/ellx+1i*Y));
unu = (ones(Ny,Nx)-CHI).*real((alpha*X+beta).*exp(nu*X/ellx+1i*Y));%heaviside(X+10).*real((alpha*X+beta).*exp(nu*X/ellx+1i*Y));

% figure(11)
% %imagesc(real(CHIp.*unu))
% plot(real(CHIp(1,:).*unu(1,:)))

% evaluate
f = ifft2(DSH(ellx,cx,mu).*fft2(w),'symmetric')...
  +DcommSH(ur,ellx,cx)+DcommSHp(unu,ellx,cx)...
   -(w.^3+3*(CHI.*ur+CHIp.*unu).*w.^2+3*(CHI.^2.*ur.^2+2*CHI.*CHIp.*ur.*unu+CHIp.^2.*unu.^2).*w+CHIp.^3.*unu.^3)...
   -3*CHI.*CHIp.*unu.*ur.*(CHIp.*unu+CHI.*ur)+(CHI-CHI.^3).*ur.^3;
%   -(w+CHI.*ur + CHIp.*unu).^3 +CHI.*ur.^3;



%    -(w.^3+3*(CHI.*ur+CHIp.*unu).*w.^2+3*(CHI.^2.*ur.^2+CHIp.^2.*unu.^2).*w+CHI.^3.*ur.^3+CHIp.^3.*unu.^3)...
%   -CHIp.*ur.^3;
     %-(w+CHI.*ur + CHIp.*unu).^3 +CHI.*ur.^3;
  %or (w+CHI.*ur+CHIp.*unu).^3 - CHI.*ur.^3;
d0r = -1+mu+cx*nur - 2*nur^2 - nur^4 + 2*nui^2 + 6*nur^2*nui^2 - nui^4;
d0i = -cx*ellx + nui*(cx-4*nur*(1+nur^2-nui^2));
d1r = cx - 4*nur*(1+nur^2 - 3*nui^2);
d1i = -4*nui*(1+3*nur^2 - nui^2);
%b0 =sum(w(:,end-7:end),'all');
%b0=sum(w.*CHIp,'all');
b0 = w(3*end/4,end)+w(3*end/4,end-1);
%sum(w(:,end)+w(:,end-1),'all');%
%b1 =  sum(w(:,end).*sin(Y(:,end)),'all') ;
%b1 = sum(w.*sin(Y).*((X>(0-2*pi)).*(X<0)),"all");
b1 = sum((0*w.*sin(Y)+ w.*cos(Y)).*((X>((5/2)*pi-2*pi)).*(X<(5/2)*pi)),"all");
%b1 =sum(w.*exp(nur*X/ellx).*cos(nui*X/ellx+Y).*((X>(20-2*pi)).*(X<20)),'all');
%b1 =sum(w.*unu.*CHIp,'all');
%b1 =sum(w.*unu.*((X>(10-2*pi)).*(X<10)),'all');
%sum(w(:,end).*sin(Y(:,end))+w(:,end-1).*sin(Y(:,end-1))w(:,end-2).*sin(Y(:,end-2)),'all');%
%sum((w(:,end)+w(:,end-1)).*cos(Y(:,end)));%sum(w(end,:)+w(end-1,:));
p0 = sum(exp(-(X+1.5).^2).*(w+1*CHI.*ur+1*CHIp.*unu),'all');
p1 = sum(phase.*ur1.*w,'all');

% figure(12)
% %imagesc(real(CHIp.*unu))
% %plot(real(CHIp(1,:).*unu(1,:)))
% imagesc(unu)%DcommSHp(unu,ellx,cx))
% colorbar

% reshape as column vector, append phase condition and secant condition
%f=[reshape(f,Nx*Ny,1);sum(phase.*ur1.*w,'all');sum(([reshape(w,Nx*Ny,1);ellx;cx]-winit).*sece)];

f = [reshape(f,N,1);d0r;d0i;d1r;d1i;b0;b1;p0;p1;sum(([reshape(w,N,1);ar;ai;br;bi;nur;nui;ellx;cx;mu]-winit).*sece)];

end
