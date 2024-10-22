function [w,w1,wold,sec,Nx,Ny]=grid_refine(w,w1,wold,sec,Nx,Ny,flag);


% double grid in x ad y
N = Nx*Ny;
wp=w(N+1:end);
wp1=w1(N+1:end);
secp=sec(N+1:end);
woldp=wold(N+1:end);

w=reshape(w(1:N),Ny,Nx);
w1=reshape(w1(1:N),Ny,Nx);
sec=reshape(sec(1:N),Ny,Nx);
wold=reshape(wold(1:N),Ny,Nx);

if flag=='xy'
    wh=fft2(w);
    wh=[wh(:,1:Nx/2),zeros(Ny,Nx),wh(:,Nx/2+1:Nx)]; % refine in x
    wh=[wh(1:Ny/2,:);zeros(Ny,2*Nx);wh(Ny/2+1:Ny,:)]; % refine in y
    w=ifft2(wh,'symmetric');
    wh=fft2(w1);
    wh=[wh(:,1:Nx/2),zeros(Ny,Nx),wh(:,Nx/2+1:Nx)]; % refine in x
    wh=[wh(1:Ny/2,:);zeros(Ny,2*Nx);wh(Ny/2+1:Ny,:)]; % refine in y
    w1=ifft2(wh,'symmetric');
    wh=fft2(sec);
    wh=[wh(:,1:Nx/2),zeros(Ny,Nx),wh(:,Nx/2+1:Nx)]; % refine in x
    wh=[wh(1:Ny/2,:);zeros(Ny,2*Nx);wh(Ny/2+1:Ny,:)]; % refine in y
    sec=ifft2(wh,'symmetric');
    wh=fft2(wold);
    wh=[wh(:,1:Nx/2),zeros(Ny,Nx),wh(:,Nx/2+1:Nx)]; % refine in x
    wh=[wh(1:Ny/2,:);zeros(Ny,2*Nx);wh(Ny/2+1:Ny,:)]; % refine in y
    wold=ifft2(wh,'symmetric');

%      w=Nx/n*Ny/n/(Lx/pi*Ly/pi)*w;
    w=4*w;
    w1=4*w1;
    sec=4*sec;
    wold=4*wold;
    Nx=Nx*2;
    Ny=Ny*2;
end
if flag=='x'
    wh=fft2(w);
    wh=[wh(:,1:Nx/2),zeros(Ny,Nx),wh(:,Nx/2+1:Nx)]; % refine in x
    w=ifft2(wh,'symmetric');
    wh=fft2(w1);
    wh=[wh(:,1:Nx/2),zeros(Ny,Nx),wh(:,Nx/2+1:Nx)]; % refine in x
    w1=ifft2(wh,'symmetric');
    wh=fft2(sec);
    wh=[wh(:,1:Nx/2),zeros(Ny,Nx),wh(:,Nx/2+1:Nx)]; % refine in x
    sec=ifft2(wh,'symmetric');
    wh=fft2(wold);
    wh=[wh(:,1:Nx/2),zeros(Ny,Nx),wh(:,Nx/2+1:Nx)]; % refine in x
    wold=ifft2(wh,'symmetric');

%      w=Nx/n*Ny/n/(Lx/pi*Ly/pi)*w;
    w=2*w;
    w1=2*w1;
    sec=2*sec;
    wold=2*wold;
    Nx=Nx*2;
end
if flag=='y'
    wh=fft2(w);
    wh=[wh(1:Ny/2,:);zeros(Ny,Nx);wh(Ny/2+1:Ny,:)]; % refine in y
    w=ifft2(wh,'symmetric');
    wh=fft2(w1);
    wh=[wh(1:Ny/2,:);zeros(Ny,Nx);wh(Ny/2+1:Ny,:)]; % refine in y
    w1=ifft2(wh,'symmetric');
    wh=fft2(sec);
    wh=[wh(1:Ny/2,:);zeros(Ny,Nx);wh(Ny/2+1:Ny,:)]; % refine in y
    sec=ifft2(wh,'symmetric');
    wh=fft2(wold);
    wh=[wh(1:Ny/2,:);zeros(Ny,Nx);wh(Ny/2+1:Ny,:)]; % refine in y
    wold=ifft2(wh,'symmetric');

%      w=Nx/n*Ny/n/(Lx/pi*Ly/pi)*w;
    w=2*w;
    w1=2*w1;
    sec=2*sec;
    wold=2*wold;
    Ny=Ny*2;
end

w=[reshape(w,Nx*Ny,1);wp];
w1=[reshape(w1,Nx*Ny,1);wp1];
sec=[reshape(sec,Nx*Ny,1);secp];
wold=[reshape(wold,Nx*Ny,1);woldp];


end
