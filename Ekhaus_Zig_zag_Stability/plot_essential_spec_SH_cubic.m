% Plot leading eigenvalue of the Floquet operator for the cubic SH eqn
% B(sigma;mu,k) = -(1 + k^2*(dz + i*sigma)^2)^2 + mu - 3*up^2
% Spatial coordinates: x direction
nz = 20; Lz = pi; hz = 2*pi/nz;  z = hz*(1:nz); z = Lz*(z-pi)/pi;
% Fourier differentiation matrix first order for t between -pi and pi
column = [0 .5*(-1).^(1:nz-1).*cot((1:nz-1)*hz/2)]';
Dz  = toeplitz(column,column([1 nz:-1:2]));
D2z = toeplitz([-pi^2/(3*hz^2)-1/6 .5*(-1).^(2:nz)./sin(hz*(1:nz-1)/2).^2]);
D4z = D2z^2;

wz = 2*pi*ones(1,nz)/nz; % integration weights for trapzoid rule - mean
z = z';
% Dx = -Dx;
Iz = speye(nz);

mesh_params.nz  = nz;  mesh_params.Lz  = Lz;  mesh_params.z   = z;    
mesh_params.Iz  = Iz;  mesh_params.wz  = wz;  
mesh_params.Dz  = Dz; mesh_params.D2z = D2z; mesh_params.D4z = D4z;

mu = 0.1;
c = 0;
k = 1.12; % unstable

p(1) = mu;
p(2) = k;

% initial guess
w0 = 0.9*cos(z(:));

mesh_params.w0 = w0; % reference profiles
mesh_params.w0z= Dz*w0;

u0 = [w0; c];

figure;plot(z,w0);title('Initial guess');xlabel('z');ylabel('u');drawnow;

% converge periodic orbit
my_rhs = @(u) SH_1D_cubic(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',500,'Algorithm','levenberg-marquardt');
[u_out,fval,exitflag,output,jacobian] = fsolve(my_rhs,u0,options);

figure; hold on;
up = u_out(1:nz);
plot(z,up);title(['converged stationary SH solution \mu= ',num2str(mu)]);xlabel('z');ylabel('up');axis tight;
% pause;

sigmax = linspace(-k/2,k/2,501);
DN = @(u) -3*u.^2;
figure;hold on;
dd = [];

for i = 1:length(sigmax)
        sigma = 1i*sigmax(i);
        L = -((k*mesh_params.Dz + sigma*speye(nz))^2+(1)*speye(nz))^2 ...
        + spdiags(ones(nz,1)*mu,0,nz,nz)  ...
        + spdiags(DN(up),0,nz,nz);
    
        d = eig(full(L));
        plot(sigmax(i),max(real(d)),'.b');xlabel('\sigma_1');ylabel('max eigenvalue');drawnow;
        dd = [dd; max(real(d))];
end

k = 1.05; % stable

p(1) = mu;
p(2) = k;

% initial guess
w0 = 0.9*cos(z(:));

mesh_params.w0 = w0; % reference profiles
mesh_params.w0z= Dz*w0;

u0 = [w0; c];

figure;plot(z,w0);title('Initial guess');xlabel('z');ylabel('u');drawnow;

% converge periodic orbit
my_rhs = @(u) SH_1D_cubic(u,p,mesh_params);
options = optimset('Jacobian','off','Display','iter','MaxIter',500,'Algorithm','levenberg-marquardt');
[u_out,fval,exitflag,output,jacobian] = fsolve(my_rhs,u0,options);

figure; hold on;
up = u_out(1:nz);
plot(z,up);title(['converged stationary SH solution \mu= ',num2str(mu)]);xlabel('z');ylabel('up');axis tight;

dd2 = [];
figure;hold on;

for i = 1:length(sigmax)
        sigma = 1i*sigmax(i);
        L = -((k*mesh_params.Dz + sigma*speye(nz))^2+(1)*speye(nz))^2 ...
        + spdiags(ones(nz,1)*mu,0,nz,nz)  ...
        + spdiags(DN(up),0,nz,nz);
    
        d = eig(full(L));
        plot(sigmax(i),max(real(d)),'.b');xlabel('\sigma_1');ylabel('max eigenvalue');drawnow;
        dd2 = [dd2; max(real(d))];
end

figure;plot(sigmax,dd,'r',sigmax,dd2,'b','LineWidth',3);xlabel('$\sigma$','Interpreter','Latex');ylabel('max eigenvalue','Interpreter','Latex');
ax = gca;
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';
legend('k=1.12','k=1.05')
axis([-0.4 0.4 -0.4 0.1])