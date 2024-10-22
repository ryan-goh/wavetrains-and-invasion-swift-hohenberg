function [w,w1,wold,sec,Nx,Lx,Xi]=grid_cut(w,w1,wold,sec,Nx,Ny,Lx,Xi);
    % halven domain length Lx, adding zeros; double grid in x
    % ends up with grid Nx/2,Ny
    Nl=round(Nx/2*(1+Xi/Lx));
    Nr=Nx-Nl;
    N = Nx*Ny;
    % filter in FOurier space to smooth out discontinuities that may arise from padding wqith zeros
    w =[w(Ny*Nx/4+1:Ny*3*Nx/4);w(N+1:end)];
    w1=[w1(Ny*Nx/4+1:Ny*3*Nx/4);w1(N+1:end)];
    wold=[wold(Ny*Nx/4+1:Ny*3*Nx/4);wold(N+1:end)];
    sec=[sec(Ny*Nx/4+1:Ny*3*Nx/4);sec(N+1:end)];
    Nx=1/2*Nx;
    Lx=1/2*Lx;
    Xi=Xi/2;
end
