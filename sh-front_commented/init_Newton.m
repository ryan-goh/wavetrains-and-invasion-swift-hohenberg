function [stepend,steps,ds,ntol,nnitermax,nminstep,ngmrestol,newtonflag,maxstepsize,gminner,maxit,arclength_max,tail_low,tail_high,xrefine,yrefine]=init_Newton;

    % tolerances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ntol=5e-4; % tolerance fort the Newton solver, max of residual
    nminstep=1e-6; % minimal Newton step, if making smaller steps without bringing the residual down give up
    ngmrestol=1e-8; % tolerance for the gmres solver; shjould be smaller than the tolerance for the Newton solver

    % max iteration counters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nnitermax= 6; % max number of Newton iterations after first initial iteration
    gminner=800; % max numer of inner iteration of gmres; limited by memory;
                    % 1500 inner iterations, 128*8192 OK, 256*8192 not OK; limit of variable size is  2048*128,8192 (but would run out of memory), 2^31
    maxit=8; % max number of gmres restarts; does not seem to help much but may need this if computing on GPU with less memory

    % initialize flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newtonflag=0; % initialize the error flag from Newton as "no issue"

    % continuation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ds=1e-1;% arclength stepsize
    maxstepsize=1e1; % max ds in secant continuation (gets weighted with the grid size to account for L^2 norm as one of the components
    arclength_max=1e-1; % maximal arclength distance in wavenumber space
    steps=1; % initialize step counter
    stepend=1500; % max number of continuation steps

    % domain and grid refinement %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tail_low=1e-10; % reduce domain if tail gets too small, below this threshold
    tail_high=1e-2; % extend domain if tail gets to large, above this threshold
%      dxmax=.3;   % grid refinement when grid coarser --- needs to be small when del is small;
%      dymax=.3;   % grid refinement when grid coarser
    yrefine=1e-1;   % max abs value in last quarter of Fourier modes in y allowed before grid refinement sets in
    xrefine=1e2;   % max abs value in last quarter of Fourier modes in x allowed before grid refinement sets in
    % coarsening starts when amplitude in last half of Fourier modes is less than xreinfe/100 (or yrefine/100)
end
