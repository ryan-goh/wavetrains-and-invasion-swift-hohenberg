function plot_2D_essential_spec_SH_cubic_JR(u_in,p,mesh_params,nev,sigmax,tauy)
% Plot leading eigenvalue of the Floquet operator for the cubic SH eqn
% B(sigma;mu,k) = -(1 + k^2*(dz + i*sigma)^2)^2 + mu - 3*up^2

mu = p(1);
c = 0;
k = p(2); 

nz = mesh_params.nz;
up = u_in(1:nz);

%sigmax = linspace(-k/2,k/2,501);
DN = @(u) -3*u.^2;
taun=numel(tauy);
dd2 = zeros(length(sigmax),nev);

figure; hold on;
for j = 1:taun
    for i = 1:length(sigmax)
        sigma = 1i*sigmax(i); tau=tauy(j);
        L = -(k^2*(mesh_params.Dz + sigma*speye(nz))^2+(1-tau^2)*speye(nz))^2 ...
            + mu*speye(nz) ... %    + spdiags(ones(nz,1)*mu,0,nz,nz)  ...
            + spdiags(DN(up),0,nz,nz);
    
        d = eig(full(L)); [~,ri]=sort(-real(d)); d=d(ri);
        %plot(sigmax(i),max(real(d)),'.b');xlabel('\sigma_1');ylabel('max eigenvalue');drawnow;
        %dd2 = [dd2; max(real(d))];
        dd2(i,:)=real(d(1:nev));
    end
    plot(sigmax,dd2,'Color',[0 0 1-(j/taun)],'LineWidth',3);
end

%figure;plot(sigmax,dd2(:,1),'b',sigmax,dd2(:,2),'m','LineWidth',3);
xlabel('$\sigma$','Interpreter','Latex');
ylabel('eigenvalues','Interpreter','Latex');
ax = gca;
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';
%legend('k=1.12','k=1.05')
% axis([-0.4 0.4 -0.4 0.1])