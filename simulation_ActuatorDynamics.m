close all
clear
%% Population Dynamics with Actuator Dynamics
% this script extends the Population Dynamics system with integrator
% dynamics for the dilution rate u(t)
%       \dot u(t) = w(t),
% where w(t) is the new input. Simulate the new dynamics with the
% backstepping-based controller.

%% ------ parameters

% [Schmidt17]
A = 2; % max age
mu = @(a) .1; % mortality rate fcn
k = @(a) 2*a.*(A-a); % birth kernel
p = 1; % output kernel
manuallyProvideMuINT = false; % boolean, that switches integral of mu on or off.
u_star = 1; % steady-state dilution rate
y0 = 1; % initial output
c1 = -.066;
c2 = -.9;
x0 = @(a) c1*a + exp(c2*a); % IC
sigma(1) = -4.0335; % eigenvalues of the form lambda = -sigma/A+-j*omega/(2*pi*A)
omega(1) = 55.4606;
sigma(2) = -4.9866;
omega(2) = 95.7048;

% % [KurthSawodny21]
% A = 2; % max age
% mu = @(a) 1/(20-5*a); % mortality rate function - problem: matlab cannot find the correct integral...
% k = @(a) a; % birth kernel
% p = 1; % output kernel
% manuallyProvideMuINT = true; % boolean, that switches integral of mu on or off.
% mu_int = @(a) -log((4-a)/4)/5; % = int_0^a mu(s) ds for a \in [0,2]
% u_star = 0.4837;
% y0 = 1; % initial output
% x0 = @(a) 8-3*a; % IC
% sigma(1) = -1.8224; % eigenvalues of the form lambda = -sigma/A+-j*omega/(2*pi*A)
% omega(1) = 48.0574;
% sigma(2) = -2.3838;
% omega(2) = 87.8539;

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

N_EV = length(sigma); % number of nonzero eigenvalues considered

sign_ImaginaryPart = 1; % only works for +1
EV = -sigma/A + 1i*omega/(2*pi*A)*sign_ImaginaryPart;

%% ------ find integral of mortality rate symbolically

if ~manuallyProvideMuINT
    syms a
    mu_sym = mu(a);
    mu_int_sym = int(mu_sym,0,a);
    mu_int = matlabFunction(mu_int_sym);
end

%% get discretization

parameter.A = A; % max age - double
parameter.mu = mu; % mortality rate - function
parameter.mu_int = mu_int; % mortality rate integral - function
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

%% simulate system - P-controller stabilizing setpoint
% here, with controller w(t) = w_cancel(t) + w_stabilize(t)
% notice that y(t) = C*lambda(t)
% denote the simulation state by rho(t) = [lambda'(t),u(t)]

% choose desired setpoint for output - equivalent to choosing a desired
% equilibrium profile x^\ast(a), or better its family parameter.
% y_des = 1.5;
y_des = .5;

% --- define controller - backstepping type
c = 1; % control gain c > 0

w_cancel = @(rho) -rho(end)-p/(C_mat*rho(1:end-1))...
            *(rho(1:end-1)'*(eval_phi(phi,A)-eval_phi(phi,0)))-mu(0);
w_stabilize = @(rho) -c*(rho(end)-u_star-log(C_mat*rho(1:end-1)/y_des));
w_ctrl = @(rho) w_cancel(rho) + w_stabilize(rho);

dynamics = @(t,rho) [(A_mat-eye(size(A_mat))*rho(end))*rho(1:end-1);
                      w_ctrl(rho)];

lambda_0 = zeros(size(A_mat,1),1);
lambda_0(end) = 1;
rho_0 = [lambda_0;u_star];
tspan = [0 20];

[t_sample,rho_sample] = ode45(dynamics,tspan,rho_0);

lambda_sample = rho_sample(:,1:end-1);
u_sample = rho_sample(:,end);

y_sample = C_mat*lambda_sample';

w_ctrl_sample = zeros(size(t_sample));
w_stabilize_sample = zeros(size(t_sample));
w_cancel_sample = zeros(size(t_sample));
for kk = 1:size(lambda_sample,1)
    w_ctrl_sample(kk) = w_ctrl(rho_sample(kk,:)');
    w_stabilize_sample(kk) = w_stabilize(rho_sample(kk,:)');
    w_cancel_sample(kk) = w_cancel(rho_sample(kk,:)');
end

%% plot results - P-controller stabilizing setpoint

figure
tiledlayout(2,2)
nexttile
plot(t_sample,lambda_sample)
title('discretized states $\lambda(t)$ - backstepping controller')
xlabel('time t')
grid on

output_ax_handle = nexttile;
hold on
plot(t_sample,ones(size(y_sample))*y_des,'--k','Linewidth',1.5)
plot(t_sample,y_sample)
title('output $y(t)$ - backstepping controller')
legend('desired output $y_\mathrm{des}$','output $y(t)$')
xlabel('time $t$')
grid on
str = {'steady-state error $\Delta y = y_\mathrm{ss} - y_\mathrm{des}$ = '...
    ,num2str(y_sample(end)-y_des)};
text(max(xlim), min(ylim),str, 'Horiz','right', 'Vert','top')

% REMARK: notation in 
% - documentation is Dilution rate D(t), voltage input u(t)
% - simulations is Dilution rate u(t), voltage input w(t)

nexttile
hold on
plot(t_sample,ones(size(u_sample))*u_star,'--k','Linewidth',1.5)
plot(t_sample,u_sample)
plot(t_sample,w_ctrl_sample)
plot(t_sample,w_cancel_sample)
plot(t_sample,w_stabilize_sample)
title('dilution rate $u(t)$ and input $w(t)$ - backstepping controller')
legend('steady state dilution $u^\ast$','dilution rate $u(t)$','input $w(t) = w_\mathrm c(t) +w_\mathrm s(t)$',...
    'cancelling terms $w_\mathrm c(t)$','stabilizing terms $w_\mathrm s(t)$')
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
title('population density $x(t,a)$ - backstepping controller')


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
