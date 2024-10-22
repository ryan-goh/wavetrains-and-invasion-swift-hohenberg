function U2=roll_ext(u,n,Nx,Ny,Lx,Ly);

    U2=n*ifft2(diag(fft(u)),'symmetric'); % shear/rotate to diagonal. forms roll on 2d grid u(ellx x +  y) with size n in both directions
       % now extend to large domain
    U2=U2*kron(ones(round(Lx/pi),1),speye(n,n))'; % repeat L/pi times to get solution on large domain
    N1=round(Lx/pi)*n;
    U2h=fft2(U2);
    U2h=[U2h(:,1:N1/2),zeros(n,Nx-N1),U2h(:,N1/2+1:N1)]; % refine in x
    U2h=[U2h(1:n/2,:);zeros(Ny-n,Nx);U2h(n/2+1:n,:)]; % refine in y
    U2=ifft2(U2h,'symmetric');
    U2=Nx/n*Ny/n/(Lx/pi*Ly/pi)*U2;
    % now refine grid

end
