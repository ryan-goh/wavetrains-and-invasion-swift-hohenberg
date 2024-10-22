%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code locates the Eckhaus stability boundary of the cubic SH equation
% Code does this in 3 stages:
% 1. Computes just a spatially periodic orbit - via a Fourier
% psuedo-spectral code 
% 2. Computes the eigenvalues/eigenvectors associated with the EK boundary
% 3. Computes the EK boundary via solving for (up, w, w_σ,w_σσ,λ_σ,λ_σσ):
% −(1 + k^2*d_zz)^2[up] + μ*up − up^3 + cd_z(up) = 0, (1.15a)
% L(up)[w] = 0, (1.15b)
% L(up)[w_σ] − λ_σ*w − 4(1 + k^2*d_zz)*k^2*d_z[w] = 0, (1.15c)
% L(up)[w_σσ] - λ_σσ*w + 2*λ_σ*w_σ + 8*(1 + k^2*d_zz)*k^2*d_z[w_σ] + 4*k^2*w + 12*k^4*d_zz[w] =0, (1.15d) 
% where L(up) = −(1 + k^2*d_zz)^2 + μ − 3*up^2
% with phase conditions
% int_0^{2pi} (u_z^{ref}*(u^{ref}-up))dz = 0 (1.16a)
% int_0^{2pi} w^2 dz = 1
% int_0^{2pi} w*w_σ dz = 0
% int_0^{2pi} w*w_σσ dz = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc;

% Spatial coordinates: z direction
nz = 40; Lz = pi; hz = 2*pi/nz;  z = hz*(1:nz); z = Lz*(z-pi)/pi;
% Fourier psuedo-spectral differentiation matrices
column = [0 .5*(-1).^(1:nz-1).*cot((1:nz-1)*hz/2)]';
Dz  = toeplitz(column,column([1 nz:-1:2]));
D2z = Dz^2; %JR = toeplitz([-pi^2/(3*hz^2)-1/6 .5*(-1).^(2:nz)./sin(hz*(1:nz-1)/2).^2]);
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
k = 1.1;

p(1) = mu;
p(2) = k;

% initial guess
w0 = cos(z(:));
mesh_params.w0 = w0; % reference profiles
mesh_params.w0z= Dz*w0;

u0 = [w0; c];

figure;plot(z,w0);title('Initial guess');xlabel('z');ylabel('u');drawnow;

%% converge to up
disp('Converging periodic orbit up');
my_rhs = @(u) SH_1D_cubic(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',500,'Algorithm','levenberg-marquardt');
[u_out,fval,exitflag,output,jacobian] = fsolve(my_rhs,u0,options);

figure; hold on;
up = u_out(1:nz);
plot(z,up);title('converged stationary SH solution \mu=0.1');xlabel('z');ylabel('u');

%% converge to EK eigenvalue problem
disp('Converging EK eigenvalue problem')
[W,D] = eigs(jacobian(1:nz,1:nz),10,0.1);
u1 = [u_out(1:nz);  W(:,1); 0*ones(nz,1); 0*ones(nz,1); 0; D(1,1); 0; 0];
my_eig = @(u) SH_1D_eig(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',500,'Algorithm','levenberg-marquardt');
[u_out2,fval,exitflag,output,jacobian] = fsolve(my_eig,u1,options);

%% converge to EK boundary problem
disp('Converging EK boundary')
u2 = u_out2;
u2(end) = p(2);
my_eig = @(u) SH_1D_eig_2(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',500,'Algorithm','levenberg-marquardt');
[u_out3,fval,exitflag,output,jacobian] = fsolve(my_eig,u2,options);
p(2) = u_out3(end)

%% Continue in the bifurcation parameter mu using Daniele's continuation code
% Various function handles
problemHandle            = @(u,p)  SH_1D_eig_2(u,p,mesh_params);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution_cubic(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures_cubic(step,u,p,mesh_params);
computeEigenvaluesHandle = [];%@(u,p) ComputeEigenvalues(u,p,t);
plotSpetcrumHandle       = []; %@(d,p,parentHandle) PlotSpectrum(d,p,parentHandle); 
stepperPars.iContPar      = 1;
stepperPars.s0            = 0.01; % minus sign here steps backwards in parameter
stepperPars.sMin          = 1e-8;
stepperPars.sMax          = .1;
stepperPars.pMin          = -1.0;
stepperPars.pMax          = 5;
stepperPars.maxSteps      = 20000;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 100;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'Jacobian','off',...
                                     'MaxIter',15);
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = 'Eckhaus_run_1';
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
s=load('Eckhaus_run_1/solution_0000117.mat');
u_ini=s.u; p_ini=s.p;
stepperPars.dataFolder    = 'Eckhaus_run_2';
stepperPars.uzstop= @(v1,v0) uzstop_val(v1,v0,10);
stepperPars.pMax          = 10;
branch2 = SecantContinuation(problemHandle,u_ini,p_ini,stepperPars);
%% go into other direction
stepperPars.s0            = -0.01; % minus sign here steps backwards in parameter
stepperPars.maxSteps      = 400;
stepperPars.nPrint        = 20;
stepperPars.dataFolder    = 'Eckhaus_run_1-';
branch = SecantContinuation(problemHandle,u_out3,p,stepperPars);
