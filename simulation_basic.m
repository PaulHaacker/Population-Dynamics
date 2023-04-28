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

par_system.A = A; % max age - double
par_system.mu = mu; % mortality rate - function
par_system.mu_int = mu_int; % mortality rate integral - function
par_system.k = k; % birth kernel - function handle
par_system.p = p; % output kernel - double
par_system.D_star = D_star; % steady-state dilution rate - double

% parameters for IC
par_system.x0 = x0; % function handle

% eigenvalues of the form lambda = -sigma/A+-j*omega/(2*pi*A)

par_system.sigma(1) = sigma(1);
par_system.omega(1) = omega(1);
par_system.sigma(2) = sigma(2);
par_system.omega(2) = omega(2);

[A_mat, C_mat, phi] = getDiscretization(par_system);

par_disc.A_mat = A_mat;
par_disc.C_mat = C_mat;
par_disc.phi = phi;

%% simulate system ODE45 - P-controller stabilizing setpoint
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
dt = .1; % only relevant for RK4 simulation

sim_par.y_des = y_des;
sim_par.y_des_d = y_des_d;
sim_par.D_ctrl = D_ctrl;
sim_par.dynamics = dynamics;
sim_par.lambda_0 = lambda_0;
sim_par.tspan = tspan;
sim_par.dt = dt;

sim_method = 'ODE45';
resultsODE45 = sim_system(par_system, par_disc, sim_par, sim_method);

%% plot results ODE45 - P-controller stabilizing setpoint - DIAGNE plot
plot_results(par_system, par_disc, resultsODE45,sim_method)

%% simulate system RK4 manually - P-controller stabilizing setpoint
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

sim_par.y_des = y_des;
sim_par.y_des_d = y_des_d;
sim_par.D_ctrl = D_ctrl;
sim_par.dynamics = dynamics;
sim_par.lambda_0 = lambda_0;
sim_par.tspan = tspan;
sim_par.dt = dt;

sim_method = 'RK4';
resultsRK4 = sim_system(par_system, par_disc, sim_par, sim_method);

%% plot results RK4 manually - P-controller stabilizing setpoint - DIAGNE plot
% plot_results(par_system, par_disc, resultsRK4,sim_method)

%% simulate system, euler scheme manually - P-controller stabilizing setpoint
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
    
% input delay
delay = .5;

dynamics = @(t,lambda) (A_mat-eye(size(A_mat))*D_ctrl(t,lambda))*lambda;

lambda_0 = zeros(size(A_mat,1),1);
lambda_0(end) = 1;
tspan = [0 15];
dt = .01;

sim_par.y_des = y_des;
sim_par.y_des_d = y_des_d;
sim_par.D_ctrl = D_ctrl;
sim_par.dynamics = dynamics;
sim_par.lambda_0 = lambda_0;
sim_par.tspan = tspan;
sim_par.dt = dt;
sim_par.delay = delay;

sim_method = 'euler';
resultsEuler = sim_system(par_system, par_disc, sim_par, sim_method);

%% plot results euler manually - P-controller stabilizing setpoint - DIAGNE plot
% plot_results(par_system, par_disc, resultsEuler,sim_method)

%% simulate delayed system, euler scheme manually - P-controller stabilizing setpoint
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
% D_ctrl = @(t,lambda) max([D_star - y_des_d(t)./y_des(t)...
%         + log(C_mat*lambda./y_des(t)),0]); % linear P-gain - dynamic FF +
% %       saturation.
D_ctrl = @(t,lambda) D_star - y_des_d(t)./y_des(t)...
        + log(C_mat*lambda./y_des(t)); % linear P-gain - dynamic FF
    
% input delay
delay = 5;

dynamics = @(t,lambda) (A_mat-eye(size(A_mat))*D_ctrl(t,lambda))*lambda;

lambda_0 = zeros(size(A_mat,1),1);
lambda_0(end) = 1;
tspan = [0 15];
dt = .01; % needs to be small for numeric stability

sim_par.y_des = y_des;
sim_par.y_des_d = y_des_d;
sim_par.D_ctrl = D_ctrl;
sim_par.dynamics = dynamics;
sim_par.lambda_0 = lambda_0;
sim_par.tspan = tspan;
sim_par.dt = dt;
sim_par.delay = delay;
sim_par.A_mat = A_mat;

sim_method = 'euler_delay';
resultsEuler_delay = sim_system(par_system, par_disc, sim_par, sim_method);

%% plot results delay euler manually - P-controller stabilizing setpoint - DIAGNE plot
plot_results(par_system, par_disc, resultsEuler_delay,sim_method)


%% debug plot

% figure
% tiledlayout(2,1)
% nexttile % absolute results
% title('absolute comparing ODE45, RK4, euler')
% hold on
% plot(resultsODE45.t_sample,resultsODE45.y_sample,'.-')
% plot(resultsRK4.t_sample,resultsRK4.y_sample,'.-')
% plot(resultsEuler.t_sample,resultsEuler.y_sample,'.-')
% legend(['ODE45 -- \# of steps ',num2str(length(resultsODE45.t_sample))],...
%     ['RK4 -- \# of steps ',num2str(length(resultsRK4.t_sample))],...
%     ['euler -- \# of steps ',num2str(length(resultsEuler.t_sample))])
% grid on
% 
% nexttile
% title('relative of RK4, euler to ODE45')
% hold on
% plot(resultsODE45.t_sample,abs(resultsODE45.y_sample ...
%     - interp1(resultsRK4.t_sample,resultsRK4.y_sample,resultsODE45.t_sample'))...
%     ./abs(resultsODE45.y_sample),'.-')
% plot(resultsODE45.t_sample,abs(resultsODE45.y_sample ...
%     - interp1(resultsEuler.t_sample,resultsEuler.y_sample,resultsODE45.t_sample'))...
%     ./abs(resultsODE45.y_sample),'.-')
% legend(['RK4 to ODE45 -- \# of steps ',num2str(length(resultsRK4.t_sample))],...
%     ['euler to ODE45 -- \# of steps ',num2str(length(resultsEuler.t_sample))])
% grid on

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
t_vec = t_0:dt:t_1;
n_steps = length(t_vec);
x_vec = zeros(n_steps,length(x_0));
x_vec(1,:)=x_0;
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

function [t_vec,x_vec]=euler(sim_par,t_range,dt,x_0)
% % inputs:
% dynamics as fcn f(t,x)
% initial time t_0
% final time t_1
% (fixed) step size dt
% initial condition x_0
% % outputs:
% state trajectory x_vec 

% extract parameters
y_des = sim_par.y_des;
y_des_d = sim_par.y_des_d;
D_ctrl = sim_par.D_ctrl;
f = sim_par.dynamics;
lambda_0 = sim_par.lambda_0;
tspan = sim_par.tspan;
dt = sim_par.dt;

% Initilization
t_0 = min(t_range);
t_1 = max(t_range);
if ~isrow(x_0) % turn initial state into row, also modify function
    x_0 = x_0';
    f = @(t,x) f(t,x')';
end

t_vec = t_0:dt:t_1;
n_steps = length(t_vec);
x_vec = zeros(n_steps,length(x_0));
x_vec(1,:)=x_0;
if length(t_vec) ~= n_steps
    error('length error')
end
i = 1;

% Euler-Scheme 1st order
for t = t_vec(2:end)
    i = i+1;
    
%     k_1 = dt * f((t-dt),x_vec(i-1,:));
%     k_2 = dt * f((t-dt) + .5*dt,x_vec(i-1,:) + .5*k_1);
%     k_3 = dt * f((t-dt) + .5*dt,x_vec(i-1,:) + .5*k_2);
%     k_4 = dt * f((t-dt) + dt,x_vec(i-1,:) + k_3);
    
    x_vec(i,:) = x_vec(i-1,:) + dt * f(t-dt,x_vec(i-1,:));
end
end

function [t_vec,x_vec]=euler_delay(sim_par,t_range,dt,x_0)
% % inputs:
% dynamics as fcn f(t,x)
% initial time t_0
% final time t_1
% (fixed) step size dt
% initial condition x_0
% % outputs:
% state trajectory x_vec 

% extract parameters
y_des = sim_par.y_des;
y_des_d = sim_par.y_des_d;
D_ctrl = sim_par.D_ctrl;
f = sim_par.dynamics;
lambda_0 = sim_par.lambda_0;
tspan = sim_par.tspan;
dt = sim_par.dt;
delay = sim_par.delay;
A_mat = sim_par.A_mat;

% Initilization
t_0 = min(t_range);
t_1 = max(t_range);
if ~isrow(x_0) % turn initial state into row, also modify function
    x_0 = x_0';
    f = @(t,x) f(t,x')';
end

t_vec = t_0:dt:t_1;
n_steps = length(t_vec);
x_vec = zeros(n_steps,length(x_0));
x_vec(1,:)=x_0;
if length(t_vec) ~= n_steps
    error('length error')
end
i = 1;

% Euler-Scheme 1st order
for t = t_vec(2:end)
    i = i+1;
    
    x_stage = x_vec(i-1,:);
    t_stage = t-dt;
    if t_stage>=t_0+delay % delayed input
        t_delayed = t_stage-delay;
        x_delayed = interp1(t_vec,x_vec,t_delayed); % interpolate, since time steps dont neccessarily match up.
        f_stage = x_stage*(A_mat-eye(size(A_mat))*D_ctrl(t_delayed,x_delayed'))';
    else
        f_stage = x_stage*A_mat'; % D_ctrl = 0 for t_stage <= t0 + delay
    end
    
%     k_1 = dt * f((t-dt),x_vec(i-1,:));
%     k_2 = dt * f((t-dt) + .5*dt,x_vec(i-1,:) + .5*k_1);
%     k_3 = dt * f((t-dt) + .5*dt,x_vec(i-1,:) + .5*k_2);
%     k_4 = dt * f((t-dt) + dt,x_vec(i-1,:) + k_3);
    
    x_vec(i,:) = x_stage + dt * f_stage;
%     x_vec(i,:) = x_vec(i-1,:) + dt * f(t-dt,x_vec(i-1,:));
end
end

function [] = plot_results(par_system, par_disc, results, sim_method)
% functionalized plot.

% % extract data
A = par_system.A; % max age - double
mu = par_system.mu; % mortality rate - function
mu_int = par_system.mu_int; % mortality rate integral - function
k = par_system.k; % birth kernel - function handle
p = par_system.p; % output kernel - double
D_star = par_system.D_star; % steady-state dilution rate - double

% % parameters for IC
% x0 = par_system.x0; % function handle

% % eigenvalues of the form lambda = -sigma/A+-j*omega/(2*pi*A)
% 
% par_system.sigma(1) = sigma(1);
% par_system.omega(1) = omega(1);
% par_system.sigma(2) = sigma(2);
% par_system.omega(2) = omega(2);

% par_disc.A_mat = A_mat;
% par_disc.C_mat = C_mat;
phi = par_disc.phi;

t_sample = results.t_sample;
lambda_sample = results.lambda_sample;
y_sample = results.y_sample;
y_des = results.y_des;
D_sample = results.D_sample;

% % plot results ODE54 - P-controller stabilizing setpoint - DIAGNE plot
fig_handle = figure;
switch sim_method
    case 'ODE45'
        sgtitle('Population Dynamics with P-controller -- ODE45','Interpreter','Latex')
    case 'RK4'
        sgtitle('Population Dynamics with P-controller -- RK4','Interpreter','Latex')
    case 'euler'
        sgtitle('Population Dynamics with P-controller -- euler','Interpreter','Latex')
    case 'euler w/ delay'
        sgtitle('Population Dynamics with P-controller and input delay -- euler','Interpreter','Latex')
end

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
end

function results = sim_system(par_system, par_disc, sim_par, sim_method)
% % extract data
A = par_system.A; % max age - double
mu = par_system.mu; % mortality rate - function
mu_int = par_system.mu_int; % mortality rate integral - function
k = par_system.k; % birth kernel - function handle
p = par_system.p; % output kernel - double
D_star = par_system.D_star; % steady-state dilution rate - double

% % parameters for IC
% x0 = par_system.x0; % function handle

A_mat = par_disc.A_mat;
C_mat = par_disc.C_mat;
phi = par_disc.phi;

y_des = sim_par.y_des;
y_des_d = sim_par.y_des_d;
D_ctrl = sim_par.D_ctrl;
dynamics = sim_par.dynamics;
lambda_0 = sim_par.lambda_0;
tspan = sim_par.tspan;
dt = sim_par.dt;

% % simulate
switch sim_method
    case 'ODE45'
        [t_sample,lambda_sample] = ode45(dynamics,tspan,lambda_0);
    case 'RK4'
        [t_sample,lambda_sample] = runge_kutta_4(dynamics,tspan,dt,lambda_0);
    case 'euler'
        [t_sample,lambda_sample] = euler(sim_par,tspan,dt,lambda_0);
    case 'euler_delay'
        [t_sample,lambda_sample] = euler_delay(sim_par,tspan,dt,lambda_0);
end

y_sample = C_mat*lambda_sample';

D_sample = zeros(size(t_sample));
for kk = 1:size(lambda_sample,1)
    D_sample(kk) = D_ctrl(t_sample(kk),lambda_sample(kk,:)');
end

results.t_sample = t_sample;
results.lambda_sample = lambda_sample;
results.y_sample = y_sample;
results.D_sample = D_sample;
results.y_des = y_des;
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
%% simulate system - Steady State Input

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
