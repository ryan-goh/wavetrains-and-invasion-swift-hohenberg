% compute zig-zag instability boundary of the cubic SH equation

clear all, close all, clc;

% Spatial coordinates: z direction
nz = 40; Lz = pi; hz = 2*pi/nz;  z = hz*(1:nz); z = Lz*(z-pi)/pi;
% Fourier psuedo-spectral differentiation matrices
column = [0 .5*(-1).^(1:nz-1).*cot((1:nz-1)*hz/2)]';
Dz  = toeplitz(column,column([1 nz:-1:2]));
D2z = toeplitz([-pi^2/(3*hz^2)-1/6 .5*(-1).^(2:nz)./sin(hz*(1:nz-1)/2).^2]);
D4z = D2z^2;
wz  = 2*pi*ones(1,nz)/nz; % integration weights for trapzoid rule - mean
z   = z';
Iz = speye(nz);

mesh_params.nz  = nz;  mesh_params.Lz  = Lz;  mesh_params.z   = z;    
mesh_params.Iz  = Iz;  mesh_params.wz  = wz;  
mesh_params.Dz  = Dz;  mesh_params.D2z = D2z; mesh_params.D4z = D4z;

 % Swift-Hohenberg equation parameters
mu = 0.1;
c = 0;
k = 0.95; 

p(1) = mu;
p(2) = k;

% initial guess
w0 = cos(z(:));
mesh_params.w0 = w0; % reference profiles
mesh_params.w0z= Dz*w0;

u0 = [w0; c];

figure;plot(z,w0);title('Initial guess');xlabel('z');ylabel('u');drawnow;

%% converge to wavetrain up
disp('Converging periodic orbit up');
my_rhs = @(u) SH_1D_cubic(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',500,'Algorithm','levenberg-marquardt');
[u_out,fval,exitflag,output,jacobian] = fsolve(my_rhs,u0,options);

figure; hold on;
up = u_out(1:nz);
plot(z,up);title('converged stationary SH solution \mu=0.1');xlabel('z');ylabel('u');

%% converge zig-zag boundary problem
disp('Converging zig-zag boundary')
u0 = [u_out;k];
my_rhs = @(u) SH_1D_cubic_zig_zag(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',500,'Algorithm','levenberg-marquardt');
[u_out2,fval,exitflag,output,jacobian] = fsolve(my_rhs,u0,options);

%% Continue in the bifurcation parameter mu using Daniele's continuation code
% Various function handles
% Note computation for small k has significant numerical error
problemHandle            = @(u,p)  SH_1D_cubic_zig_zag(u,p,mesh_params);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution_cubic(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures_cubic(step,u,p,mesh_params);
computeEigenvaluesHandle = [];
plotSpetcrumHandle       = [];
stepperPars.iContPar      = 1;
stepperPars.s0            = -0.01; % minus sign here steps backwards in parameter
stepperPars.sMin          = 1e-8;
stepperPars.sMax          = .1;
stepperPars.pMin          = -1.0;
stepperPars.pMax          = 10;
stepperPars.maxSteps      = 20000;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 100;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'Jacobian','off',...
                                     'MaxIter',15);
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = 'Zig_zag_run_1';
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;      
stepperPars.PlotBranchVariableId = 2; % plot the selected k value
stepperPars.uzstop= @(v1,v0) uzstop_val(v1,v0,10);
branch = SecantContinuation(problemHandle,u_out2,p,stepperPars);
