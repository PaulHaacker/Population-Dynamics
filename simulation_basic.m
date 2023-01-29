close all
clear

%% ------ parameters

A = 2; % max age
mu = .1; % constant mortality rate
k = @(a) 2*a.*(A-a); % birth kernel
p = 1; % output kernel
u_star = 1; % steady-state dilution rate
y0 = 1; % initial output

% parameters for IC - paper [Schmidt17]
c1 = -.066;
c2 = -.9;
x0 = @(a) c1*a + exp(c2*a);

% % parameters for IC - Thesis [Schmidt16]
% mu_par = 1.3;
% b_par = .066;
% x0_par = 1.256;
% x0 = @(a) (-b_par*a + x0_par * exp(-mu_par*a));

% % parameters for IC - Thesis [Schmidt16] --> does not seem to work
% nu_par = 1.3;
% b_par = 12*k(0)/A^4/nu_par^3*(A*nu_par*exp(-A*nu_par)+2*exp(-A*nu_par)-2+A*nu_par);
% x0_par = 6/A^2/nu_par^3/(1-exp(-nu_par*A))*(A*nu_par*exp(-A*nu_par)+2*exp(-A*nu_par)-2+A*nu_par+2*y0*nu_par);
% x0 = @(a) (-b_par*a + x0_par * exp(-nu_par*a));

%% eigenvalues
% of the form lambda = -sigma/A+-j*omega/(2*pi*A)

% % values suggested by paper in [Schmidt17] eq. (49)-(50)
% sigma(1) = 1.010;
% omega(1) = .351;
% sigma(2) = 1.250;
% omega(2) = .601;

% % values suggested by paper in [Schmidt17] figure 2
% sigma(1) = 4.1;
% omega(1) = 37.69;
% sigma(2) = 5;
% omega(2) = 62.83;

% % values I suspect by paper in [Schmidt17] eq. (49)-(50)
% % -> they match the plots - note that eigenfuctions are still
% eigenfunctions with flipped signs.
% sigma(1) = 4.04;
% omega(1) = 55.43;
% sigma(2) = 5;
% omega(2) = 94.92;

% values I found, BUT with flipped signs of real parts
% -> 
sigma(1) = -4.0335;
omega(1) = 55.4606;
sigma(2) = -4.9866;
omega(2) = 95.7048;

N_EV = 2; % number of nonzero eigenvalues considered

sign_ImaginaryPart = 1; % only works for +1
EV = -sigma/A + 1i*omega/(2*pi*A)*sign_ImaginaryPart;

%% get discretization

parameter.A = A; % max age - double
parameter.mu = mu; % constant mortality rate - double
parameter.k = k; % birth kernel - function handle
parameter.p = p; % output kernel - double
parameter.u_star = u_star; % steady-state dilution rate - double

% parameters for IC - paper [Schmidt17]
parameter.x0 = x0; % function handle

% eigenvalues of the form lambda = -sigma/A+-j*omega/(2*pi*A)

parameter.sigma(1) = sigma(1);
parameter.omega(1) = omega(1);
parameter.sigma(2) = sigma(2);
parameter.omega(2) = omega(2);

[A_mat, C_mat, phi] = getDiscretization(parameter);

%% simulate linear system - Steady State Input
% here, with steady state input u(t) == u_star
% denote the simulation state by lambda

dynamics = @(t,lambda) (A_mat-eye(size(A_mat))*u_star)*lambda;

lambda_0 = zeros(size(A_mat,1),1);
lambda_0(end) = 1;
tspan = [0 20];

[t_sample,lambda_sample] = ode45(dynamics,tspan,lambda_0);

y_sample = C_mat*lambda_sample';

%% plot results - Steady State Input

figure
tiledlayout(2,2)
nexttile
plot(t_sample,lambda_sample)
title('discretized states $\lambda(t)$ - steady state input $u(t) = u^\ast$')
xlabel('time $t$')
grid on

nexttile
plot(t_sample,y_sample)
title('output $y(t)$ - steady state input $u(t) = u^\ast$')
xlabel('time $t$')
grid on

% plot the PDE state where x(t,a) = lambda(t)'*phi(a);

% time sample from above
% define domain sample
a_sample = 0:0.1:A;

[a_mesh,t_mesh] = meshgrid(a_sample,t_sample);
x_mesh = zeros(size(a_mesh));

for ii = 1:length(t_sample)
    for jj = 1:length(a_sample)
        x_mesh(ii,jj) = lambda_sample(ii,:)*eval_phi(phi,a_sample(jj));
    end
end

axes_handle = nexttile;
surf_plot = surf(a_mesh,t_mesh,x_mesh);
LessEdgeSurf(surf_plot);

xlabel('age $a$')
ylabel('time $t$')
title('population density $x(t,a)$ - steady state input $u(t) = u^\ast$')

axes_handle.CameraPosition = [15.7853   91.8902    2.8718];

%% simulate linear system - P-controller stabilizing setpoint
% here, with controller u(t) == u_star + ln(y(t)/y_des)
% notice that y(t) = C*lambda(t)
% denote the simulation state by lambda

% choose desired setpoint for output - equivalent to choosing a desired
% equilibrium profile x^\ast(a), or better its family parameter.

y_des = .5;

% u_ctrl = @(lambda) u_star + log(C_mat*lambda/y_des);
u_ctrl = @(lambda) u_star + (C_mat*lambda-y_des)/y_des;

dynamics = @(t,lambda) (A_mat-eye(size(A_mat))*u_ctrl(lambda))*lambda;

lambda_0 = zeros(size(A_mat,1),1);
lambda_0(end) = 1;
tspan = [0 20];

[t_sample,lambda_sample] = ode45(dynamics,tspan,lambda_0);

y_sample = C_mat*lambda_sample';

u_sample = zeros(size(t_sample));
for kk = 1:size(lambda_sample,1)
    u_sample(kk) = u_ctrl(lambda_sample(kk,:)');
end

%% plot results - P-controller stabilizing setpoint
figure
tiledlayout(2,2)
nexttile
plot(t_sample,lambda_sample)
title('discretized states $\lambda(t)$ - logarithmic P-controller stabilizing setpoint')
xlabel('time t')
grid on

nexttile
hold on
plot(t_sample,ones(size(y_sample))*y_des,'--k','Linewidth',1.5)
plot(t_sample,y_sample)
title('output $y(t)$ - logarithmic P-controller stabilizing setpoint')
legend('desired output $y_{des}$','output $y(t)$')
xlabel('time $t$')
grid on

nexttile
hold on
plot(t_sample,ones(size(u_sample))*u_star,'--k','Linewidth',1.5)
plot(t_sample,u_sample)
title('control input $u(t)$ - logarithmic P-controller stabilizing setpoint')
legend('steady state input $u^\ast$','input $u(t)$')
xlabel('time $t$')
grid on

% plot the PDE state where x(t,a) = lambda(t)'*phi(a);

% time sample from above
% define domain sample
a_sample = 0:0.1:A;

[a_mesh,t_mesh] = meshgrid(a_sample,t_sample);
x_mesh = zeros(size(a_mesh));

for ii = 1:length(t_sample)
    for jj = 1:length(a_sample)
        x_mesh(ii,jj) = lambda_sample(ii,:)*eval_phi(phi,a_sample(jj));
    end
end

axes_handle = nexttile;
surf_plot = surf(a_mesh,t_mesh,x_mesh);
LessEdgeSurf(surf_plot);
axes_handle.CameraPosition = [15.7853   91.8902    2.8718];

xlabel('age $a$')
ylabel('time $t$')
title('population density $x(t,a)$ - logarithmic P-controller stabilizing setpoint')


%% functions

function val = eval_phi(phi,a)
% takes N-by-1 cell of functions phi and evaluates them at a, returns
% values as N-by-1 matrix
N = length(phi);
val = zeros(size(phi));
for kk = 1:N
    val(kk) = phi{kk}(a);
end
end
