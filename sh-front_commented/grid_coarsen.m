function [w,w1,wold,sec,Nx,Ny]=grid_coarsen(w,w1,wold,sec,Nx,Ny,flag);
    N = Nx*Ny;
    % halven grid in x and/or y
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
        wh=[wh(:,1:Nx/4),wh(:,3*Nx/4+1:Nx)]; % coarsen in x
        wh=[wh(1:Ny/4,:);wh(3*Ny/4+1:Ny,:)]; % coarsen in y
        w=ifft2(wh,'symmetric');
        w=1/4*w;
        wh=fft2(w1);
        wh=[wh(:,1:Nx/4),wh(:,3*Nx/4+1:Nx)]; % coarsen in x
        wh=[wh(1:Ny/4,:);wh(3*Ny/4+1:Ny,:)]; % coarsen in y
        w1=ifft2(wh,'symmetric');
        w1=1/4*w1;
        wh=fft2(sec);
        wh=[wh(:,1:Nx/4),wh(:,3*Nx/4+1:Nx)]; % coarsen in x
        wh=[wh(1:Ny/4,:);wh(3*Ny/4+1:Ny,:)]; % coarsen in y
        sec=ifft2(wh,'symmetric');
        sec=1/4*sec;
        wh=fft2(wold);
        wh=[wh(:,1:Nx/4),wh(:,3*Nx/4+1:Nx)]; % coarsen in x
        wh=[wh(1:Ny/4,:);wh(3*Ny/4+1:Ny,:)]; % coarsen in y
        wold=ifft2(wh,'symmetric');
        wold=1/4*wold;
        Nx=Nx/2;
        Ny=Ny/2;
    end
    if flag=='x'
        wh=fft2(w);
        wh=[wh(:,1:Nx/4),wh(:,3*Nx/4+1:Nx)]; % coarsen in x
        w=ifft2(wh,'symmetric');
        w=1/2*w;
        wh=fft2(w1);
        wh=[wh(:,1:Nx/4),wh(:,3*Nx/4+1:Nx)]; % coarsen in x
        w1=ifft2(wh,'symmetric');
        w1=1/2*w1;
        wh=fft2(sec);
        wh=[wh(:,1:Nx/4),wh(:,3*Nx/4+1:Nx)]; % coarsen in x
        sec=ifft2(wh,'symmetric');
        sec=1/2*sec;
        wh=fft2(wold);
        wh=[wh(:,1:Nx/4),wh(:,3*Nx/4+1:Nx)]; % coarsen in x
        wold=ifft2(wh,'symmetric');
        wold=1/2*wold;
        Nx=Nx/2;
    end
    if flag=='y'
        wh=fft2(w);
        wh=[wh(1:Ny/4,:);wh(3*Ny/4+1:Ny,:)]; % coarsen in y
        w=ifft2(wh,'symmetric');
        w=1/2*w;
        wh=fft2(w1);
        wh=[wh(1:Ny/4,:);wh(3*Ny/4+1:Ny,:)]; % coarsen in y
        w1=ifft2(wh,'symmetric');
        w1=1/2*w1;
        wh=fft2(sec);
        wh=[wh(1:Ny/4,:);wh(3*Ny/4+1:Ny,:)]; % coarsen in y
        sec=ifft2(wh,'symmetric');
        sec=1/2*sec;
        wh=fft2(wold);
        wh=[wh(1:Ny/4,:);wh(3*Ny/4+1:Ny,:)]; % coarsen in y
        wold=ifft2(wh,'symmetric');
        wold=1/2*wold;
        Ny=Ny/2;
    end
    w=[reshape(w,Nx*Ny,1);wp];
    w1=[reshape(w1,Nx*Ny,1);wp1];
    sec=[reshape(sec,Nx*Ny,1);secp];
    wold=[reshape(wold,Nx*Ny,1);woldp];

end
