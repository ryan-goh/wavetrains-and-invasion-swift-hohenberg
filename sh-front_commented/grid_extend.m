function [w,w1,wold,sec,Nx,Lx,Xi]=grid_extend(w,w1,wold,sec,Nx,Ny,Lx,Xi);

% double domain length Lx, adding zeros; double grid in x
% ends up with grid 2*Nx,Ny
    Nl=round(Nx/2*(1+Xi/Lx));
    Nr=Nx-Nl;
    N = Nx*Ny;
    % filter in FOurier space to smooth out discontinuities that may arise from padding wqith zeros
    smooth=ones(Ny,2*Nx); smooth(:,7/16*(2*Nx):9/16*(2*Nx))=0;
    % now pad with zeros and smooth
    w=     [reshape(ifft2(smooth.*fft2([zeros(Ny,Nl),reshape(   w(1:N),Ny,Nx),zeros(Ny,Nr)]),'symmetric') ,Ny*2*Nx,1);    w(N+1:end)];
    w1=    [reshape(ifft2(smooth.*fft2([zeros(Ny,Nl),reshape(  w1(1:N),Ny,Nx),zeros(Ny,Nr)]),'symmetric') ,Ny*2*Nx,1);   w1(N+1:end)];
    wold=  [reshape(ifft2(smooth.*fft2([zeros(Ny,Nl),reshape(wold(1:N),Ny,Nx),zeros(Ny,Nr)]),'symmetric') ,Ny*2*Nx,1); wold(N+1:end)];
    sec=   [reshape(ifft2(smooth.*fft2([zeros(Ny,Nl),reshape( sec(1:N),Ny,Nx),zeros(Ny,Nr)]),'symmetric') ,Ny*2*Nx,1);  sec(N+1:end)];
    Nx=2*Nx;
    Lx=2*Lx;
    Xi=2*Xi;
end
