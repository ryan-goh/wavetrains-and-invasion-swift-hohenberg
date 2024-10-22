function [CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit]= movedatatogpu(CHI,Kx,Ky,sec,ur0,w0,X,Y,phase,x,k,Lx,Ly,del,cx,beta,mu,ellxgrid,ellygrid,gminner,ngmrestol,maxit)
    beta = gpuArray(beta);
    CHI=gpuArray(CHI);
    Kx=gpuArray(Kx);
    Ky=gpuArray(Ky);
    sec=gpuArray(sec);
    ur0=gpuArray(ur0);
    w0=gpuArray(w0);
    X=gpuArray(X);
    Y=gpuArray(Y);
    phase=gpuArray(phase);
    x=gpuArray(x);
    k=gpuArray(k);
    Lx=gpuArray(Lx);
    Ly =gpuArray(Ly);
    del=gpuArray(del);
    cx=gpuArray(cx);
    mu=gpuArray(mu);
    ellxgrid=gpuArray(ellxgrid);
    ellygrid=gpuArray(ellygrid);
    gminner=gpuArray(gminner);
    ngmrestol=gpuArray(ngmrestol);
    maxit=gpuArray(maxit);
end
