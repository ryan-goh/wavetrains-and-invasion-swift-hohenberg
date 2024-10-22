%% plot both the upper Eckhaus and Zig-zag boundaries
% load zig-zag and Eckhaus boundary data
zz = load('Zig_zag_run_1/branch.mat');
EK1= load('Eckhaus_run_1/branch.mat');
EK2= load('Eckhaus_run_2/branch.mat');
EK1n= load('Eckhaus_run_1-/branch.mat');

ff=figure(1); hold on;
plot(zz.branch(:,3),zz.branch(:,end),'-r','LineWidth',3);
plot(EK1.branch(:,3),EK1.branch(:,end),'b','LineWidth',3);
plot(EK2.branch(:,3),EK2.branch(:,end),'b','LineWidth',3);
plot(EK1n.branch(:,3),EK1n.branch(:,end),'b','LineWidth',3);

xlabel('$\mu$','Interpreter','Latex');
ylabel('$k$','Interpreter','Latex');
ax = gca;
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';
axis tight;
title('Eckhaus and Zig-zag instability boundaries')
legend('Zig-zag boundary','Eckhaus bounday');
%% add lower k boundary
EK1sn= load('Eckhaus_run_1-_lower/branch.mat');
EK2sn= load('Eckhaus_run_1-2_lower/branch.mat');
EK1s= load('Eckhaus_run_1_lower/branch.mat');
EK2s= load('Eckhaus_run_2_lower/branch.mat');
EK3s= load('Eckhaus_run_3_lower/branch.mat');

EKsnb=[EK1sn.branch; EK2sn.branch]; 
EKsnb=sortrows([EKsnb(:,3), EKsnb(:,end)]); 
EKsb=[EK1s.branch;EK2s.branch;EK3s.branch];
EKsb=[EKsb(:,3),EKsb(:,end)];
EKsb=[EKsnb;EKsb];

plot(EKsb(:,1),EKsb(:,2),'b','LineWidth',3); 
% shading zigzag
zx=zz.branch(:,3); zy=zz.branch(:,end);
zxmi=min(zx); zxma=max(zx); zymi=min(zz.branch(:,end));
zx=[zx;zxma;zxmi;zxmi]; zy=[zy;zymi;zymi;max(zy)];
patch(zx,zy, 'r','FaceAlpha',0.25) % zigzag unstable

% shading EK
EK1b=sortrows([EK1n.branch(:,3), EK1n.branch(:,end)]); 
EK2b=[EK1.branch;EK2.branch]; EK2b=[EK2b(:,3),EK2b(:,end)];
EK1b=[EK1b;EK2b];

ex=EK1b(:,1); ey=EK1b(:,2); ex=[ex;min(ex);min(ex)]; 
eyma=max(ey); ey=[ey;eyma;min(ey)];
patch(ex,ey, 'b', 'FaceAlpha',0.25)
%
% small k
esx=EKsb(:,1); esy=EKsb(:,2);
esx=[esx;max(esx);min(esx)]; 
esy=[esy;zymi;zymi];
patch(esx,esy, 'b', 'FaceAlpha',0.25,'EdgeAlpha',0.25)

% plot markers
plot(0.1, 1.05, 'k.', 'LineWidth', 2, 'MarkerSize', 25);
plot(0.1, 1.12, 'k.', 'LineWidth', 2, 'MarkerSize', 25);

legend('Zig-zag boundary','Eckhaus bounday');
%
axis([0,9,0.6,inf]);
saveas(ff,'Busse.eps','epsc');
%% More plots:

%% plot 1D essential spectrum for largest mu on lower EK boundary: 
% here one loop near origin (due to symmetry seemingly two loops)
s=load('Eckhaus_run_3_lower/solution_0000091.mat');
s.p(2)=s.u(end); 
% set meshparams:
nz = 80; Lz = pi; hz = 2*pi/nz;  z = hz*(1:nz); z = Lz*(z-pi)/pi;
mesh_params.nz=nz; 
column = [0 .5*(-1).^(1:nz-1).*cot((1:nz-1)*hz/2)]';
mesh_params.Dz  = toeplitz(column,column([1 nz:-1:2]));

[fh,~]=plot_essential_spec_SH_cubic(s.u,s.p,mesh_params,2,linspace(-1,1,101));
saveas(fh,'spec1_up_crit_mu10_k0p9.eps','epsc');
%% add computation by continuation
b=load('Spec_x/branch.mat');
figure(fh); hold on; plot(real(b.branch(:,3)),real(b.branch(:,6)),'r','LineWidth',3);
figure(fh); axis([-1,1,-inf,inf]);
%% illustrate disconnected parts
[fh,~]=plot_essential_spec_SH_cubic(s.u,s.p,mesh_params,4,linspace(-1,1,101));
saveas(fh,'spec2_up_crit_mu10_k0p9.eps','epsc');

%% plot 2D essential spectrum: illustrate zigzag stable
s=load('Eckhaus_run_3_lower/solution_0000091.mat');
s.p(2)=s.u(end);
fh=plot_2D_essential_spec_SH_cubic(s.u,s.p,mesh_params,1,linspace(-0.5,0.5,49),linspace(0,0.1,5));
saveas(fh,'spec_zz_up_crit_mu10_k0p9.eps','epsc');
%% plot 2D essential spectrum: illustrate zigzag unstable
s=load('Eckhaus_run_3_lower/solution_0000000.mat');
s.p(2)=s.u(end);
fh=plot_2D_essential_spec_SH_cubic(s.u,s.p,mesh_params,1,linspace(-0.5,0.5,49),linspace(0,0.1,5));
saveas(fh,'spec_zz_up_crit_mu5_k0p8.eps','epsc');
