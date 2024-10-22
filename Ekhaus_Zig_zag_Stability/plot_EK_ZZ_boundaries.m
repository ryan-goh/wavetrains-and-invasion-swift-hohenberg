% plot both the Eckhaus and Zig-zag boundaries
% load zig-zag and Ekhaus boundary data
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
%% with lower k boundary
EK1sn= load('Eckhaus_run_1-_ksm/branch.mat');
EK2sn= load('Eckhaus_run_1-2_ksm/branch.mat');
EK1s= load('Eckhaus_run_1_ksm/branch.mat');
EK2s= load('Eckhaus_run_2_ksm/branch.mat');
EK3s= load('Eckhaus_run_3_ksm/branch.mat');

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
%%
saveas(ff,'Busse.eps','epsc');