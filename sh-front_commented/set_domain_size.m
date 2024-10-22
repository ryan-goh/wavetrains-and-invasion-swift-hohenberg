function [n,Nx,Ny,Lx,Ly,Xi]=set_domain_size;

    n =16*1;  % for small grid, discretization, n=16 gives error 1e-8; for larger mu may need n=32, also helpful for small kx...
    m=16*2;   % number of SH periods in large grid, should see about half of those behind the quenching line
    Nx=m*n*1; % number of Fourier modes for large grid; has the same spatial resolution as small grid initially;
    Lx=m*pi;% half the large domain size; note that the domain size is actually Lx/ellx due to rescaling
    Ny=n;   % number of Fourier modes across strip thus far need to use initially same discretization in y as in x;
    Ly=pi;  % half the domain size across; formailty since we always rescale across; actual domain size is Ly/elly
    Xi=Lx/3;% position where front interface is locked.  (should be multiple of 2*pi)

end
