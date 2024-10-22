%derivative eval
function df=f1(w,dw,CHI,CHIp,X,Y,Nx,Ny,DSH,phase,ur,ur1,dfar,dfai,dfbr,dfbi,dfnur,dfnui,dfellx,dfcx,dfmu,sece)

N = Nx*Ny;
ar = w(N+1); ai = w(N+2); br=w(N+3); bi = w(N+4); nur = w(N+5); nui = w(N+6); 
 ellx = w(N+7); cx = w(N+8);mu = w(N+9);
 alpha = ar+1i*ai;
 beta = br+1i*bi;
 nu = nur+1i*nui;

dar = dw(N+1); dai = dw(N+2); dbr=dw(N+3); dbi = dw(N+4); dnur = dw(N+5); dnui = dw(N+6); 
 dellx = dw(N+7); dcx = dw(N+8); dmu = dw(N+9);


%unu = (1./(1+exp(-5*(X+1)))).*real((alpha*X+beta).*exp(nu*X/ellx+1i*Y));
unu = (ones(Ny,Nx)-CHI).*real((alpha*X+beta).*exp(nu*X/ellx+1i*Y));%heaviside(X+10).*real((alpha*X+beta).*exp(nu*X/ellx+1i*Y));




w = reshape(w(1:N),Ny,Nx);
dw=reshape(dw(1:N),Ny,Nx); % rearrange as matrix



dfw = ifft2(DSH(ellx,cx,mu).*fft2(dw),'symmetric')-(3*(CHI.^2.*ur.^2+2*CHI.*CHIp.*ur.*unu+CHIp.^2.*unu.^2)+6*(CHI.*ur+CHIp.*unu).*w+3*w.^2).*dw;
%dfw = ifft2(DSH(ellx,cx,mu).*fft2(dw),'symmetric') - 3*(w.^3+CHI.*ur + CHIp.*unu).^2.*dw;
% dampen even modes; .3 is (artificial) damping factor
%dfw=dfw-0.1*(dw+circshift(dw,Ny/2,1));

%dd1 =-12*nu.^2 - 4;
%dd0 =  -4*nu.*(1+nu.^2).*dnu; ;

%dumb0 = zeros(Ny,Nx); db0(end-1:end,:) = 1;
 %reshape(dumb0,Nx*Ny,1);
db0 =  dw(3*end/4,end)+dw(3*end/4,end-1);%
%db0=sum(dw.*CHIp,'all');
%db1 =  sum(dw(:,end-2:end).*cos(Y(:,end-2:end)),'all');%
%db1 = sum(dw.*sin(Y).*((X>(0-2*pi)).*(X<0)),'all');
db1 = sum((0*dw.*sin(Y)+ dw.*cos(Y)).*((X>((5/2)*pi-2*pi)).*(X<(5/2)*pi)),'all');

%db1 = sum(dw.*exp(nur*X/ellx).*cos(nui*X/ellx+Y).*((X>(20-2*pi)).*(X<20)),'all');
%db1 = sum(dw.*unu.*((X>(10-2*pi)).*(X<10)),'all');
%db1 =sum(dw.*unu.*((X>(20-2*pi)).*(X<20)),'all');
%db1 = sum(dw.*unu.*CHIp,'all');


%sum(dw(:,end-2:end).*cos(Y(:,end-2:end)),'all');
%sum(dw(:,end).*sin(Y(:,end))+dw(:,end-1).*sin(Y(:,end-1))+dw(:,end-2).*sin(Y(:,end-2)),'all');
%sum((dw(:,end)+dw(:,end-1)).*cos(Y(:,end)),"all");
%sum((dw(:,end)+dw(:,end-1)).*cos(Y(:,end)));%sum(w(end,:)+w(end-1,:));
dp0 = sum(exp(-(X+1.5).^2).*dw,'all');
dp1 = sum(phase.*ur1.*dw,'all');




% flatten into column vector, append linearized phase condition and secant condition
dfw=[reshape(dfw,N,1);0;0;0;0;db0;db1;dp0;dp1; sum([reshape(dw,N,1);dar;dai;dbr;dbi;dnur;dnui;dellx;dcx;dmu].*sece,'all')];

% add in the corrections from the parameter derivatives
% size(dfw)
% size(dfb)
% size(dfnu)
% size(dfellx)
% size(dfcx)
% size(dfmu)
    df=dfw+dar*dfar + dai*dfai + dbr*dfbr + dbi*dfbi + dnur*dfnur + dnui*dfnui...
        + dellx*dfellx + dcx*dfcx + dmu*dfmu;







end
