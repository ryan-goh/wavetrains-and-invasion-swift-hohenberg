%% auto data plots

%% bif-dias
bd0p9=load("wt_bif_mu-0p9.txt",'-ascii'); 
[mv,mi]=min(abs(bd0p9(:,1)-7.2)); % sol12 near par(11)=7.2
bd0p1=load("wt_bif_mu-0p1.txt",'-ascii');
[mv,mi2]=min(abs(bd0p1(:,1)-3.3173833828)); % sol20
f1=figure(1); 
plot(pi*ones(800,1)./bd0p9(1:800,1),bd0p9(1:800,3),'-ok','LineWidth',3, ...
    'MarkerIndices',mi,'MarkerEdgeColor','red','MarkerSize',10'); 
n=numel(bd0p9(801:end,1)); 
n2=numel(bd0p1(:,1)); 
hold on; 
plot(pi*ones(n,1)./bd0p9(801:end,1),bd0p9(801:end,3),'black','LineWidth',3); 
plot(pi*ones(n2,1)./bd0p1(:,1),bd0p1(:,3),'-ob','LineWidth',3, ...
    'MarkerIndices',mi2,'MarkerEdgeColor','red','MarkerSize',10'); 
hold off;

xlabel('$k$','Interpreter','Latex')
ylabel('$\max(u)$','Interpreter','Latex')
ax = gca;
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';
%saveas(f1,'wt_bifs','pdf'); % does not crop white margins
saveas(f1,'wt_bifs.eps','epsc');
%% sols
sol12=load("wt_bif_mu-0p9_sol12.txt",'-ascii');
sol20=load("wt_bif_mu-0p1_sol20.txt",'-ascii');
f2=figure(2); 
plot(pi*(1+sol12(:,1)),sol12(:,2),'black','LineWidth',3); hold on;
plot(pi*(1-sol12(:,1)),sol12(:,2),'black','LineWidth',3);
plot(pi*(1+sol20(:,1)),sol20(:,2),'blue','LineWidth',3); 
plot(pi*(1-sol20(:,1)),sol20(:,2),'blue','LineWidth',3); hold off;
xlabel('$x$','Interpreter','Latex')
ylabel('$u(x)$','Interpreter','Latex')
%xticks([0 pi 2*pi])
ax = gca;
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';
saveas(f2,'wt_sols.eps','epsc');
