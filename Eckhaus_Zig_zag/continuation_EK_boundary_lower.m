%% See continuation_EK_boundary.m -- here small k
clear all, close all, clc;

% Spatial coordinates: z direction
nz = 80; Lz = pi; hz = 2*pi/nz;  z = hz*(1:nz); z = Lz*(z-pi)/pi;
% Fourier pseudo-spectral differentiation matrices
column = [0 .5*(-1).^(1:nz-1).*cot((1:nz-1)*hz/2)]';
Dz  = toeplitz(column,column([1 nz:-1:2]));
D2z = Dz^2; %JR =toeplitz([-pi^2/(3*hz^2)-1/6 .5*(-1).^(2:nz)./sin(hz*(1:nz-1)/2).^2]);
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
k = 0.95; % initial choice for lower EK-boundary

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
%% plot 1D essential spectrum (here stable with quadratic shape)
plot_essential_spec_SH_cubic(up,p,mesh_params,1,linspace(-1,1,101));
%% converge EK eigenvalue problem
disp('Converging EK eigenvalue problem')
[W,D] = eigs(jacobian(1:nz,1:nz),10,0.1);
u1 = [u_out(1:nz);  W(:,1); 0*ones(nz,1); 0*ones(nz,1); 0; D(1,1); 0; 0];
my_eig = @(u) SH_1D_eig(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',500,'Algorithm','levenberg-marquardt');
[u_out2,fval,exitflag,output,jacobian] = fsolve(my_eig,u1,options);
%% converge EK boundary problem
disp('Converging EK boundary')
u2 = u_out2; % u_out2(end) is lambda_sigsig
u2(end) = p(2); % according to SH_1D_eig_2 last entry should be k
my_eig = @(u) SH_1D_eig_2(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',500,'Algorithm','levenberg-marquardt');
[u_out3,fval,exitflag,output,jacobian] = fsolve(my_eig,u2,options);
k_ek = u_out3(end)
PlotSolution_cubic(u_out3,p,[],mesh_params);
%% plot 1D essential spectrum (now flat quartic shape)
p(2)=u_out3(end);
plot_essential_spec_SH_cubic(u_out3,p,mesh_params,1,linspace(-0.1,0.1,101));
%% Continue in (mu,k) using Avitabile's continuation code
% Various function handles
problemHandle            = @(u,p)  SH_1D_eig_2(u,p,mesh_params);
plotSolutionHandle       = []; %@(u,p,parentHandle) PlotSolution_cubic(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures_cubic(step,u,p,mesh_params);
computeEigenvaluesHandle = [];
plotSpetcrumHandle       = [];
stepperPars.iContPar      = 1;
stepperPars.s0            = 0.01; % minus sign here steps backwards in parameter
stepperPars.sMin          = 1e-8;
stepperPars.sMax          = .1;
stepperPars.pMin          = -1.0;
stepperPars.pMax          = 5;
stepperPars.maxSteps      = 100;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 100;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'Jacobian','off',...
                                     'MaxIter',15);
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = 'Eckhaus_run_1_lower';
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;      
stepperPars.PlotBranchVariableId = 2; % plot the selected k value
stepperPars.uzstop= @(v1,v0) uzstop_val(v1,v0,10);
branch = SecantContinuation(problemHandle,u_out3,p,stepperPars);

figure;
plot(branch(:,3),branch(:,6));xlabel('\mu');ylabel('k_{ek}');axis tight;
%% continue further
s=load('Eckhaus_run_1_lower/solution_0000100.mat');
u_ini=s.u; p_ini=s.p;
stepperPars.dataFolder    = 'Eckhaus_run_2_lower';
branch = SecantContinuation(problemHandle,u_ini,p_ini,stepperPars);
%% continue further
s=load('Eckhaus_run_2_lower/solution_0000077.mat');
u_ini=s.u; p_ini=s.p;
stepperPars.dataFolder    = 'Eckhaus_run_3_lower';
stepperPars.pMax          = 10;
branch = SecantContinuation(problemHandle,u_ini,p_ini,stepperPars);
%% continue negative direction
s=load('Eckhaus_run_1_lower/solution_0000000.mat');
stepperPars.dataFolder    = 'Eckhaus_run_1-_lower';
stepperPars.s0            = -0.01; % minus sign here steps backwards in parameter
branch = SecantContinuation(problemHandle,s.u,s.p,stepperPars);
%% continue further
s=load('Eckhaus_run_1-_lower/solution_0000100.mat');
stepperPars.dataFolder    = 'Eckhaus_run_1-2_lower';
stepperPars.s0            = -0.01; % minus sign here steps backwards in parameter
stepperPars.maxSteps      = 500;
stepperPars.nPrint        = 10;
stepperPars.nSaveSol      = 500;
branch = SecantContinuation(problemHandle,s.u,s.p,stepperPars);
%% plot 1D essential spectrum at largest mu: here one loop/interval near origin (due to symmetry seemingly two loops)
s=load('Eckhaus_run_3_lower/solution_0000091.mat');
s.p(2)=s.u(end);
[fh,~]=plot_essential_spec_SH_cubic(s.u,s.p,mesh_params,2,linspace(-1,1,101));
%% plot 2D essential spectrum: stable
s=load('Eckhaus_run_3_lower/solution_0000091.mat');
s.p(2)=s.u(end);
plot_2D_essential_spec_SH_cubic(s.u,s.p,mesh_params,1,linspace(-0.5,0.5,49),linspace(0,0.1,5));
%% plot 2D essential spectrum: unstable
s=load('Eckhaus_run_3_lower/solution_0000000.mat');
s.p(2)=s.u(end);
plot_2D_essential_spec_SH_cubic(s.u,s.p,mesh_params,1,linspace(-0.5,0.5,49),linspace(0,0.1,5));
%% essential spectrum near origin by continuation with given up
% Various function handles
s=load('Eckhaus_run_3_lower/solution_0000091.mat');
u_x_ini=[mesh_params.Dz*s.u(1:nz);0]; %u_x_ini=[s.u(1:2*nz);0;0]; % 0 for c and lam
p_x_ini=[0,s.u(4*nz+4),s.p(1)]; % sigma, k, mu
fprintf('k=%f, mu=%f\n',p_x_ini(2),p_x_ini(3));

problemHandle            = @(u,p)SH_1D_eig_x(u,p,mesh_params,s.u(1:nz));
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution_cubic(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures_cubic(step,u,p,mesh_params);
computeEigenvaluesHandle = [];
plotSpetcrumHandle       = [];
stepperPars.iContPar      = 1;
stepperPars.s0            = 0.01; % minus sign here steps backwards in parameter
stepperPars.sMin          = 1e-8;
stepperPars.sMax          = 0.02;
stepperPars.pMin          = -1.0;
stepperPars.pMax          = 1;
stepperPars.maxSteps      = 20000;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 100;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'Jacobian','off',...
                                     'MaxIter',15);
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = 'Spec_x';
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;      
stepperPars.PlotBranchVariableId = 2; % plot the lam value
stepperPars.uzstop=@(v1,v0)uzstop_val(v1,v0,1);

branch = SecantContinuation(problemHandle,u_x_ini,p_x_ini,stepperPars);
%% plot into figure of direct computation
figure(fh); hold on; plot(real(branch(:,3)),real(branch(:,6)),'r', 'LineWidth',3);
