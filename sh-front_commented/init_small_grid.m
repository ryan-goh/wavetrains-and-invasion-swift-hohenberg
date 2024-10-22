function [x,k,dsh,pc,ur0,flag]=init_small_grid(n,mu0);
    % construct small grid and differential operators
    x   = linspace(-pi, pi,n+1)'; x=x(1:end-1);
    %derivative vectors (in Fourier space)
    k   = ([[0:n/2] [-n/2+1: -1]])';
    dsh = @(ell,mu) -(ell^4*k.^4 - 2*ell^2*k.^2 + 1)+mu*ones(n,1);% main differential operator
    pc  = @(du,mu) [ifft(fft(du(1:end-1))./(dsh(1,mu)-(1+mu)*ones(n,1)),'symmetric');du(end)]; %preconditioner
    %initial shape
    ell=1;
    u = sqrt(4*(mu0-(ell-1)^2)/3)*cos(x);
    % now find roll solution
    [ur0,flag] = unit_roll(ell,mu0,u,x,k,dsh,pc);
    if norm(u)/sqrt(abs(mu0))<1e-6  flag='zero_solution'; end

end
