
function [u,ff]=unit_roll(ell,mu,u,x,k,dsh,pc)
    % n=16 is usually enough to get rolls to 1e-8
    % input
    % wavenumber ell, initial guess u, x, derivative k, dsh is sh linear part, pc is fourth order preconditioner for gmres, tol is tolerance for Newton residual, maxiter is limit on number of Newton iterations
    tol     = 1e-8;
    maxiter = 20;

    shroll=@(u) ([ifft((dsh(ell,mu)-i*u(end)*k).*fft(u(1:end-1)),'symmetric')-u(1:end-1).^3 ;sin(x)'*u(1:end-1)]); %nonlinear function with phase condition
    %                  u(1:end-1).^3;sin(x)'*u(1:end-1)];
    dshrolla= @(u,du) ...
        ([ifft((dsh(ell,mu)-i*k*u(end)).*fft(du(1:end-1)),'symmetric')-3*u(1:end-1).^2.*du(1:end-1)+ ...
        + ifft(-i*k.*fft(u(1:end-1)),'symmetric')*du(end);...
         sin(x)'*du(1:end-1)]);
    %
    pc0 = @(du) pc(du,mu);
    u=[u;0];

    % Set linear solver tolerance to be about 10 times smaller than the Newton
    gmrestol=tol*1e-1;

    % initialize first Newton step
    rhs=shroll(u);                      % compute initial residual
    residual=norm(rhs);                 % estimate norm of residual

    ff=0; % default flag output if successful

    % Newton loop
    niter = 1;
    while (residual>tol)&&(niter<maxiter)
    %
        dshroll=@(du) dshrolla(u,du);       % form Jacobian
        [nincr,flag]=gmres(dshroll,rhs,10,gmrestol,15,pc0);
                                            % gmres solve for increment
        if flag==1
            sprintf(['gmres did not converge in small domain problem, residual is ' num2str(residual) ' after ' num2str(niter) ' iterations'])
            ff='gmres';
        break
        end
        u=u-nincr;                         % Newton step
        niter=niter+1;                % keep track of number of iterations
        %
        % recompute residual
        rhs=shroll(u);                     % compute residual
        residual=norm(rhs);           % estimate residual
        %
        if niter>maxiter
            sprintf(['Maximal number of iterations reached in small domain problem, giving up; residual is ' num2str(npar.residual)])
            ff='newton';
            break
        end
    end
    u=u(1:end-1);
end
