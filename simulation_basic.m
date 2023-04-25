%% Population Dynamics
% this script runs the simulation of the age-dependent population dynamics
% as considered in [Schmidt18]
close all
clear

%% ------ parameters

% % [Schmidt18]
% A = 2; % max age
% mu = @(a) .1; % mortality rate fcn
% k = @(a) 2*a.*(A-a); % birth kernel
% p = 1; % output kernel
% manuallyProvideMuINT = false; % boolean, that switches integral of mu on or off.
% D_star = 1; % steady-state dilution rate
% y0 = 1; % initial output
% c1 = -.066;
% c2 = -.9;
% x0 = @(a) c1*a + exp(c2*a); % IC
% sigma(1) = -4.0335; % eigenvalues of the form lambda = -sigma/A+-j*omega/(2*pi*A)
% omega(1) = 55.4606;
% sigma(2) = -4.9866;
% omega(2) = 95.7048;

% [KurthSawodny21]
A = 2; % max age
mu = @(a) 1./(20-5*a); % mortality rate function - problem: matlab cannot find the correct integral...
k = @(a) a; % birth kernel
p = @(a) 1+.1*a.^2; % output kernel
manuallyProvideMuINT = true; % boolean, that switches integral of mu on or off.
mu_int = @(a) -log((4-a)/4)/5; % = int_0^a mu(s) ds for a \in [0,2]
D_star = 0.4837;
y0 = 1; % initial output
x0 = @(a) (8 - 3*a); % IC
sigma(1) = -1.8224; % eigenvalues of the form lambda = -sigma/A+-j*omega/(2*pi*A)
omega(1) = 48.0574;
sigma(2) = -2.3838;
omega(2) = 87.8539;

%% stash of unordered parameter sets

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

%% debug

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
parameter.D_star = D_star; % steady-state dilution rate - double

% parameters for IC
parameter.x0 = x0; % function handle

% eigenvalues of the form lambda = -sigma/A+-j*omega/(2*pi*A)

parameter.sigma(1) = sigma(1);
parameter.omega(1) = omega(1);
parameter.sigma(2) = sigma(2);
parameter.omega(2) = omega(2);

[A_mat, C_mat, phi] = getDiscretization(parameter);

%% simulate linear system ODE45 - P-controller stabilizing setpoint
% here, with controller u(t) == D_star + ln(y(t)/y_des)
% notice that y(t) = C*lambda(t)
% denote the simulation state by lambda

% choose desired setpoint for output - equivalent to choosing a desired
% equilibrium profile x^\ast(a), or better its family parameter.
% edit: also works for trajectories!

% static setpoint:
y_des = @(t) 12*ones(size(t));
y_des_d = @(t) 0*ones(size(t));

% % time signal:
% y_des = @(t) 12+sin(t);
% y_des_d = @(t) cos(t);

% D_ctrl = @(t,lambda) D_star + log(C_mat*lambda/y_des); % logarithmic P-gain
% D_ctrl = @(t,lambda) D_star + (C_mat*lambda-y_des)/y_des; % linear P-gain
D_ctrl = @(t,lambda) D_star - y_des_d(t)./y_des(t)...
        + log(C_mat*lambda./y_des(t)); % linear P-gain - dynamic FF


dynamics = @(t,lambda) (A_mat-eye(size(A_mat))*D_ctrl(t,lambda))*lambda;

lambda_0 = zeros(size(A_mat,1),1);
lambda_0(end) = 1;
tspan = [0 15];

[t_sample,lambda_sample] = ode45(dynamics,tspan,lambda_0);

y_sample = C_mat*lambda_sample';

D_sample = zeros(size(t_sample));
for kk = 1:size(lambda_sample,1)
    D_sample(kk) = D_ctrl(t_sample(kk),lambda_sample(kk,:)');
end

%% plot results ODE54 - P-controller stabilizing setpoint - DIAGNE plot
fig_handle = figure;
sgtitle('Population Dynamics with P-controller -- ODE45','Interpreter','Latex')

% nexttile
% plot(t_sample,lambda_sample)
% title('discretized states $\lambda(t)$')
% xlabel('time t')
% grid on

ax1 = subplot(2,6,1:3);
hold on
plot(t_sample,y_des(t_sample),'--k','Linewidth',1.5)
plot(t_sample,y_sample)
title('output $y(t)$')
legend('desired output $y_\mathrm{des}(t)$','output $y(t)$')
xlabel('time $t$')
grid on

ax2 = subplot(2,6,4:6);
hold on
plot(t_sample,ones(size(D_sample))*D_star,'--k','Linewidth',1.5)
plot(t_sample,D_sample)
title('control input $D(t)$')
legend('steady state input $D^\ast$','input $D(t)$')
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

axes_handle = subplot(2,6,8:11);
% surf_plot = surf(a_mesh,t_mesh,x_mesh,'FaceColor',[0 0.4470 0.7410]);
surf_plot = surf(a_mesh,t_mesh,x_mesh,'FaceColor','none');
LessEdgeSurf(surf_plot,20,20);
axes_handle.CameraPosition = [15.8393  -65.0215   37.0943];

xlabel('age $a$')
ylabel('time $t$')
title('population density $f(t,a)$')

%% simulate linear system RK4 manually - P-controller stabilizing setpoint
% here, with controller u(t) == D_star + ln(y(t)/y_des)
% notice that y(t) = C*lambda(t)
% denote the simulation state by lambda

% choose desired setpoint for output - equivalent to choosing a desired
% equilibrium profile x^\ast(a), or better its family parameter.
% edit: also works for trajectories!

% static setpoint:
y_des = @(t) 12*ones(size(t));
y_des_d = @(t) 0*ones(size(t));

% % time signal:
% y_des = @(t) 12+sin(t);
% y_des_d = @(t) cos(t);

% D_ctrl = @(t,lambda) D_star + log(C_mat*lambda/y_des); % logarithmic P-gain
% D_ctrl = @(t,lambda) D_star + (C_mat*lambda-y_des)/y_des; % linear P-gain
D_ctrl = @(t,lambda) D_star - y_des_d(t)./y_des(t)...
        + log(C_mat*lambda./y_des(t)); % linear P-gain - dynamic FF

dynamics = @(t,lambda) (A_mat-eye(size(A_mat))*D_ctrl(t,lambda))*lambda;

lambda_0 = zeros(size(A_mat,1),1);
lambda_0(end) = 1;
tspan = [0 15];
dt = .1;

[t_sample_RK,lambda_sample] = runge_kutta_4(dynamics,tspan,dt,lambda_0);

y_sample = C_mat*lambda_sample';

D_sample = zeros(size(t_sample_RK));
for kk = 1:size(lambda_sample,1)
    D_sample(kk) = D_ctrl(t_sample_RK(kk),lambda_sample(kk,:)');
end

%% plot results RK4 manually - P-controller stabilizing setpoint - DIAGNE plot
fig_handle = figure;
sgtitle('Population Dynamics with P-controller -- RK4 manual','Interpreter','Latex')

% nexttile
% plot(t_sample_RK,lambda_sample)
% title('discretized states $\lambda(t)$')
% xlabel('time t')
% grid on

ax1 = subplot(2,6,1:3);
hold on
plot(t_sample_RK,y_des(t_sample_RK),'--k','Linewidth',1.5)
plot(t_sample_RK,y_sample)
title('output $y(t)$')
legend('desired output $y_\mathrm{des}(t)$','output $y(t)$')
xlabel('time $t$')
grid on

ax2 = subplot(2,6,4:6);
hold on
plot(t_sample_RK,ones(size(D_sample))*D_star,'--k','Linewidth',1.5)
plot(t_sample_RK,D_sample)
title('control input $D(t)$')
legend('steady state input $D^\ast$','input $D(t)$')
xlabel('time $t$')
grid on

% plot the PDE state where x(t,a) = lambda(t)'*phi(a);

% time sample from above
% define domain sample
a_sample = 0:0.1:A;

[a_mesh,t_mesh] = meshgrid(a_sample,t_sample_RK);
x_mesh = zeros(size(a_mesh));

for ii = 1:length(t_sample_RK)
    for jj = 1:length(a_sample)
        x_mesh(ii,jj) = lambda_sample(ii,:)*eval_phi(phi,a_sample(jj));
    end
end

axes_handle = subplot(2,6,8:11);
% surf_plot = surf(a_mesh,t_mesh,x_mesh,'FaceColor',[0 0.4470 0.7410]);
surf_plot = surf(a_mesh,t_mesh,x_mesh,'FaceColor','none');
LessEdgeSurf(surf_plot,20,20);
axes_handle.CameraPosition = [15.8393  -65.0215   37.0943];

xlabel('age $a$')
ylabel('time $t$')
title('population density $f(t,a)$')


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

function [t_vec,x_vec]=runge_kutta_4(f,t_range,dt,x_0)
% % inputs:
% dynamics as fcn f(t,x)
% initial time t_0
% final time t_1
% (fixed) step size dt
% initial condition x_0
% % outputs:
% state trajectory x_vec 

% Initilization
t_0 = min(t_range);
t_1 = max(t_range);
if ~isrow(x_0) % turn initial state into row, also modify function
    x_0 = x_0';
    f = @(t,x) f(t,x')';
end
n_steps = ceil((t_1-t_0)/dt)+1;
x_vec = zeros(n_steps,length(x_0));
x_vec(1,:)=x_0;
t_vec = t_0:dt:t_1;
if length(t_vec) ~= n_steps
    error('length error')
end
i = 1;

% Runge-Kutta-Scheme 4th order
for t = t_vec(2:end)
    i = i+1;
    
    k_1 = dt * f((t-dt),x_vec(i-1,:));
    k_2 = dt * f((t-dt) + .5*dt,x_vec(i-1,:) + .5*k_1);
    k_3 = dt * f((t-dt) + .5*dt,x_vec(i-1,:) + .5*k_2);
    k_4 = dt * f((t-dt) + dt,x_vec(i-1,:) + k_3);
    
    x_vec(i,:) = x_vec(i-1,:) + (k_1 + k_4)/6 + (k_2 + k_3)/3   ;
end
end

%% stash of outdated sections
%% plot results - P-controller stabilizing setpoint
% figure
% tiles_handle = tiledlayout(2,2);
% title(tiles_handle,'Population Dynamics with P-controller stabilizing setpoint','Interpreter','Latex')
% nexttile
% plot(t_sample,lambda_sample)
% title('discretized states $\lambda(t)$')
% xlabel('time t')
% grid on
% 
% nexttile
% hold on
% plot(t_sample,ones(size(y_sample))*y_des,'--k','Linewidth',1.5)
% plot(t_sample,y_sample)
% title('output $y(t)$')
% legend('desired output $y_\mathrm{des}$','output $y(t)$')
% xlabel('time $t$')
% grid on
% 
% nexttile
% hold on
% plot(t_sample,ones(size(D_sample))*D_star,'--k','Linewidth',1.5)
% plot(t_sample,D_sample)
% title('control input $D(t)$')
% legend('steady state input $D^\ast$','input $D(t)$')
% xlabel('time $t$')
% grid on
% 
% % plot the PDE state where x(t,a) = lambda(t)'*phi(a);
% 
% % time sample from above
% % define domain sample
% a_sample = 0:0.1:A;
% 
% [a_mesh,t_mesh] = meshgrid(a_sample,t_sample);
% x_mesh = zeros(size(a_mesh));
% 
% for ii = 1:length(t_sample)
%     for jj = 1:length(a_sample)
%         x_mesh(ii,jj) = lambda_sample(ii,:)*eval_phi(phi,a_sample(jj));
%     end
% end
% 
% axes_handle = nexttile;
% % surf_plot = surf(a_mesh,t_mesh,x_mesh,'FaceColor',[0 0.4470 0.7410]);
% surf_plot = surf(a_mesh,t_mesh,x_mesh,'FaceColor','none');
% LessEdgeSurf(surf_plot,10,10);
% axes_handle.CameraPosition = [15.8393  -65.0215   37.0943];
% 
% xlabel('age $a$')
% ylabel('time $t$')
% title('population density $x(t,a)$')
%% simulate linear system - Steady State Input

% % here, with steady state input u(t) == D_star
% % denote the simulation state by lambda
% 
% D_ctrl = @(t,lambda) D_star*1.5;
% 
% dynamics = @(t,lambda) (A_mat-eye(size(A_mat))*D_ctrl(t,lambda))*lambda;
% 
% lambda_0 = zeros(size(A_mat,1),1);
% lambda_0(end) = 1;
% tspan = [0 20];
% 
% [t_sample,lambda_sample] = ode45(dynamics,tspan,lambda_0);
% 
% y_sample = C_mat*lambda_sample';
%% plot results - Steady State Input

% figure
% tiledlayout(2,2)
% nexttile
% plot(t_sample,lambda_sample)
% title('discretized states $\lambda(t)$ - steady state input $D(t) = D^\ast$')
% xlabel('time $t$')
% grid on
% 
% nexttile
% plot(t_sample,y_sample)
% title('output $y(t)$ - steady state input $D(t) = D^\ast$')
% xlabel('time $t$')
% grid on
% 
% % plot the PDE state where x(t,a) = lambda(t)'*phi(a);
% 
% % time sample from above
% % define domain sample
% a_sample = 0:0.1:A;
% 
% [a_mesh,t_mesh] = meshgrid(a_sample,t_sample);
% x_mesh = zeros(size(a_mesh));
% 
% for ii = 1:length(t_sample)
%     for jj = 1:length(a_sample)
%         x_mesh(ii,jj) = lambda_sample(ii,:)*eval_phi(phi,a_sample(jj));
%     end
% end
% 
% axes_handle = nexttile;
% surf_plot = surf(a_mesh,t_mesh,x_mesh);
% LessEdgeSurf(surf_plot,10,10);
% 
% xlabel('age $a$')
% ylabel('time $t$')
% title('population density $x(t,a)$ - steady state input $D(t) = D^\ast$')
% 
% axes_handle.CameraPosition = [15.8393  -65.0215   37.0943];
