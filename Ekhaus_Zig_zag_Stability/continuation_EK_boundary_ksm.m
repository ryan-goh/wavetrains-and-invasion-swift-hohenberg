%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code locates the Ekhaus stability boundary of the cubic SH equation
% Code does this in 3 stages:
% 1. Computes just a spatially periodic orbit - via a Fourier
% psuedo-spectral code 
% 2. Computes the eigenvalues/eigenvectors associated with the EK boundary
% 3. Computes the EK boundary via solving for (up, w, w_σ,w_σσ,λ_σ,λ_σσ):
% −(1 + k^2*d_zz)^2[up] + μ*up − up^3 + cd_z(up) = 0, (1.15a)
% L(up)[w] = 0, (1.15b)
% L(up)[w_σ] − λ_σ*w − 4(1 + k^2*d_zz)*k^2*d_z[w] = 0, (1.15c)
% L(up)[w_σσ] - λ_σσ*w + 2*λ_σ*w_σ + 8*(1 + k^2*d_zz)*k^2*d_z[w_σ] + 4*k^2*w + 12*k^2*d_zz[w] =0, (1.15d) 
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
k = 0.9;%1.0222;

p(1) = mu;
p(2) = k;

% initial guess
w0 = cos(z(:));
mesh_params.w0 = w0; % reference profiles
mesh_params.w0z= Dz*w0;

u0 = [w0; c];

figure;plot(z,w0);title('Initial guess');xlabel('z');ylabel('u');drawnow;pause;

%% converge up
disp('Converging periodic orbit up');
my_rhs = @(u) SH_1D_cubic(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',500,'Algorithm','levenberg-marquardt');
[u_out,fval,exitflag,output,jacobian] = fsolve(my_rhs,u0,options);

figure; hold on;
up = u_out(1:nz);
plot(z,up);title('converged stationary SH solution \mu=0.1');xlabel('z');ylabel('u');

%% converge EK eigenvalue problem
disp('Converging EK eigenvalue problem')
[W,D] = eigs(jacobian(1:nz,1:nz),10,0.1);
u1 = [u_out(1:nz);  W(:,1); 0*ones(nz,1); 0*ones(nz,1); 0; D(1,1); 0; 0];
my_eig = @(u) SH_1D_eig(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',500,'Algorithm','levenberg-marquardt');
[u_out2,fval,exitflag,output,jacobian] = fsolve(my_eig,u1,options);

%% converge EK boundary problem
disp('Converging EK boundary')
u2 = u_out2;
u2(end) = p(2);
my_eig = @(u) SH_1D_eig_2(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',500,'Algorithm','levenberg-marquardt');
[u_out3,fval,exitflag,output,jacobian] = fsolve(my_eig,u2,options);
k_ek = u_out3(end)

%% Continue in (mu,k) using Daniele's continuation code
% Various function handles
problemHandle            = @(u,p)  SH_1D_eig_2(u,p,mesh_params);
plotSolutionHandle       = []; %@(u,p,parentHandle) PlotSolution_cubic(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures_cubic(step,u,p,mesh_params);
computeEigenvaluesHandle = [];%@(u,p) ComputeEigenvalues(u,p,t);
plotSpetcrumHandle       = []; %@(d,p,parentHandle) PlotSpectrum(d,p,parentHandle); 
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
stepperPars.dataFolder    = 'Eckhaus_run_2_ksm';
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;      
stepperPars.PlotBranchVariableId = 2; % plot the selected k value
branch = SecantContinuation(problemHandle,u_out3,p,stepperPars);

figure;
plot(branch(:,3),branch(:,6));xlabel('\mu');ylabel('k_{ek}');axis tight;
%%
s=load('Eckhaus_run_2_ksm/solution_0000100.mat');
u_ini=s.u; p_ini=s.p;
stepperPars.dataFolder    = 'Eckhaus_run_2_ksm2';
branch = SecantContinuation(problemHandle,u_ini,p_ini,stepperPars);
%%
s=load('Eckhaus_run_2_ksm2/solution_0000100.mat');
u_ini=s.u; p_ini=s.p;
stepperPars.dataFolder    = 'Eckhaus_run_2_ksm3';
branch = SecantContinuation(problemHandle,u_ini,p_ini,stepperPars);
%%
s=load('Eckhaus_run_2_ksm3/solution_0000100.mat');
u_ini=s.u; p_ini=s.p;
stepperPars.dataFolder    = 'Eckhaus_run_2_ksm4';
stepperPars.maxSteps      = 500;
stepperPars.nPrint        = 10;
stepperPars.nSaveSol      = 500;
branch = SecantContinuation(problemHandle,u_ini,p_ini,stepperPars);
%%
s=load('Eckhaus_run_2_ksm4/solution_0000500.mat');
u_ini=s.u; p_ini=s.p;
stepperPars.dataFolder    = 'Eckhaus_run_2_ksm5';
branch = SecantContinuation(problemHandle,u_ini,p_ini,stepperPars);
%%
figure(1);
EKm= load('Eckhaus_run_2_ksm-/branch.mat');
EK= load('Eckhaus_run_2_ksm/branch.mat');
EK2= load('Eckhaus_run_2_ksm2/branch.mat');
EK3= load('Eckhaus_run_2_ksm3/branch.mat');
EK4= load('Eckhaus_run_2_ksm4/branch.mat');
EK5= load('Eckhaus_run_2_ksm5/branch.mat');
plot(EKm.branch(:,3),EKm.branch(:,end),'b','LineWidth',3);
hold on; 
plot(EK.branch(:,3),EK.branch(:,end),'r','LineWidth',3);
plot(EK2.branch(:,3),EK2.branch(:,end),'k','LineWidth',3);
plot(EK3.branch(:,3),EK3.branch(:,end),'k','LineWidth',3);
plot(EK4.branch(:,3),EK4.branch(:,end),'k','LineWidth',3);
plot(EK5.branch(:,3),EK5.branch(:,end),'k','LineWidth',3);
hold off;
%%
s=load('Eckhaus_run_2_ksm2/solution_0000100.mat');
u_ini=s.u; p_ini=[s.p(1), u_ini(end)]; u_ini(end)=0;
maxre=max_essential_spec_SH_cubic_JR(u_ini,p_ini,mesh_params,linspace(-0.01,0.01,10))
%%
s=load('Eckhaus_run_2_ksm5/solution_0000500.mat');
u_ini=s.u; p_ini=[s.p(1), u_ini(end)]; u_ini(end)=0;
plot_essential_spec_SH_cubic_JR(u_ini,p_ini,mesh_params,2,linspace(-0.01,0.01,101));
%% Continue in with free curvature in mu
% Various function handles
problemHandle            = @(u,p)  SH_1D_eig_JR(u,p,mesh_params);
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
stepperPars.maxSteps      = 10;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 1;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'Jacobian','off',...
                                     'MaxIter',15);
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = 'Eckhaus_run_3';
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;    
stepperPars.PlotBranchVariableId = 2; % not how this is set
branch = SecantContinuation(problemHandle,u_ini,p_ini,stepperPars);








%% converge EK eigenvalue problem
disp('Converging EK eigenvalue problem')
[W,D] = eigs(jacobian(1:nz,1:nz),10,0.1);
u1 = [u_out(1:nz);  W(:,1); 0*ones(nz,1); 0*ones(nz,1); 0; D(1,1); 0; 0];
my_eig = @(u) SH_1D_eig(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',500,'Algorithm','levenberg-marquardt');
[u_out2,fval,exitflag,output,jacobian] = fsolve(my_eig,u1,options);
