close all
clear
%% Population Dynamics with Actuator Dynamics
% this script extends the Population Dynamics system with integrator
% dynamics for the dilution rate D(t)
%       \dot D(t) = u(t),
% where u(t) is the new input. Simulate the new dynamics with the
% backstepping-based controller.

%% booleans for plots
debug_plot = false; % plot for debug
standard_plot = true; % standard plot
transformed_plot = true; % plot of transformed coordinates

%% switch for control
ctrl_mode = 'Karafyllis';'Backstepping';  % options: 'Backstepping', 'Karafyllis'

%% ------ parameters

% [Schmidt17]
A = 2; % max age
mu = @(a) .1; % mortality rate fcn
k = @(a) 2*a.*(A-a); % birth kernel
p = @(a) 1; % output kernel
manuallyProvideMuINT = false; % boolean, that switches integral of mu on or off.
D_star = 1; % steady-state dilution rate
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
% mu = @(a) 1./(20-5*a); % mortality rate function - problem: matlab cannot find the correct integral...
% k = @(a) a; % birth kernel
% p = @(a) 1+.1*a.^2; % output kernel
% manuallyProvideMuINT = true; % boolean, that switches integral of mu on or off.
% mu_int = @(a) -log((4-a)/4)/5; % = int_0^a mu(s) ds for a \in [0,2]
% D_star = 0.4837;
% y0 = 1; % initial output
% x0 = @(a) 8-3*a; % IC
% sigma(1) = -1.8224; % eigenvalues of the form lambda = -sigma/A+-j*omega/(2*pi*A)
% omega(1) = 48.0574;
% sigma(2) = -2.3838;
% omega(2) = 87.8539;

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

% N_EV = length(sigma); % number of nonzero eigenvalues considered
% 
% sign_ImaginaryPart = 1; % only works for +1
% EV = -sigma/A + 1i*omega/(2*pi*A)*sign_ImaginaryPart;

%% ------ symbolic treatment of parameters

% integral of mortality rate needed for eigenfunctions in discretization
% find integral of mortality rate symbolically
if ~manuallyProvideMuINT
    syms a
    mu_sym = mu(a);
    mu_int_sym = int(mu_sym,0,a);
    mu_int = matlabFunction(mu_int_sym);
end

% derivative of birth kernel needed for controller
syms p_sym(a)
p_sym(a) = p(a);
p_prime_sym = diff(p_sym,a);
p_prime = matlabFunction(p_prime_sym);

%% get discretization

par_sys.A = A; % max age - double
par_sys.mu = mu; % mortality rate - function
par_sys.mu_int = mu_int; % mortality rate integral - function
par_sys.k = k; % birth kernel - function handle
par_sys.p = p; % output kernel - function handle
par_sys.D_star = D_star; % steady-state dilution rate - double

% parameters for IC - paper [Schmidt17]
par_sys.x0 = x0; % function handle

% eigenvalues of the form lambda = -sigma/A+-j*omega/(2*pi*A)

par_sys.sigma(1) = sigma(1);
par_sys.omega(1) = omega(1);
par_sys.sigma(2) = sigma(2);
par_sys.omega(2) = omega(2);

[A_mat, C_mat, phi] = getDiscretization(par_sys);

discretization.A_mat = A_mat;
discretization.C_mat = C_mat;
discretization.phi = phi;

%% define P-controller stabilizing setpoint
% here, with controller u(t) = u_cancel(t) + u_stabilize(t)
% notice that y(t) = C*lambda(t)
% denote the simulation state by rho(t) = [lambda'(t),D(t)]

% choose desired setpoint for output - equivalent to choosing a desired
% equilibrium profile x^\ast(a), or better its family parameter.
% y_des = 1.5;
y_des = 20;

% allowed interval of dilution rate
D_min = 0.5; % minimum Dilution rate constraint for Safety-Filter
D_max = 1.5; % maximum Dilution rate constraint for Safety-Filter

switch ctrl_mode
    case 'Karafyllis'
        % --- Karafyllis' Transformation Controller (TC)
        % parameter
        Q_TC = 1; % TC gain parameter
        sigma_TC = 1; % TC gain parameter
        k_TC = 1; % TC gain parameter
        D_min_TC = D_min; % minimum Dilution rate constraint for TC
        D_max_TC = D_max; % maximum Dilution rate constraint for TC
        A_TC = D_max_TC - D_min_TC; % TC parameter
        B_TC = (D_max_TC - D_star)/(D_star - D_min_TC); % TC parameter

        zeta_fcn = @(D) log(B_TC*(D -D_min_TC)/(D_max_TC - D));
        exp_zeta = @(D) exp(zeta_fcn(D));

        if D_star < D_min_TC || D_star > D_max_TC
            error('equilibrium dilution not within constraints')
        end

        u_ctrl = @(rho) u_ctrl_TC_fcn(rho, zeta_fcn, exp_zeta, A_TC, ...
                                B_TC, k_TC, sigma_TC, Q_TC, y_des, C_mat);

    case 'Backstepping'
        % --- Backstepping Controller
        % controller parameters
        c = 2; % control gain c > 0
        k_safety = 1; % 2*(c+1); % safety gain "rate of allowed safety dissipation"
        D_min_safe = D_min; % minimum Dilution rate constraint for Safety-Filter
        D_max_safe = D_max; % maximum Dilution rate constraint for Safety-Filter
        h_fcn = @(D) -(D-D_min_safe).*(D-D_max_safe); % safety function for D_min_safe <= D(t) <= D_max_safe
        L_g_h = @(D) -2*D + D_max_safe + D_min_safe; % lie derivative of h(D) = -(D-D_min_safe)(D-D_max_safe) along g(D) == 1;
        int_par = zeros(size(phi)); % numerical representation of the integral part of the cancelling controller
        for ii = 1:length(phi)
            int_par(ii) = integral(@(a)(p_prime(a)-p(a).*mu(a)).*phi{ii}(a),0,A);
        end

        if D_star < D_min_safe || D_star > D_max_safe
            error('equilibrium dilution not within constraints')
        end

        % --- define controller - backstepping type

        u_cancel = @(rho) -rho(end)-1/(C_mat*rho(1:end-1))...
                    *(p(A)*eval_phi(phi,A)-p(0)*eval_phi(phi,0)-int_par)'*rho(1:end-1); % cancelling terms - exact cancellation
        % u_cancel = @(rho) -rho(end)+D_star; % cancelling terms - super late botching (seems to work)
        % u_cancel = @(rho) 0; % cancelling terms - early botching (seems to work)
        % u_cancel = @(rho) -rho(end)+D_star*y_des/(C_mat*rho(1:end-1)); % cancelling terms - late botching (unstable)
        u_stabilize = @(rho) -c*(rho(end)-D_star-log(C_mat*rho(1:end-1)/y_des)); % stabilizing terms

        % --- safety override
%         u_constraint = @(rho) 0; % ignore constraints on D(t)
        % u_constraint =  @(rho) -log(rho(end)/D_star); % logarithmic penalty of D(t)->0
        % u_constraint =  @(rho) (y_des)/4*(- rho(end) + D_star); % linear penalty of D(t)->0
%         u_constraint =  @(rho) max(0,- u_cancel(rho) - u_stabilize(rho) ...
%                         +k_safety*(-rho(end)+D_min_safe)); % Safety-Filter for D(t) > D_min_safe
        % u_constraint =  @(rho) max(0,-(u_cancel(rho) + u_stabilize(rho))*L_g_h(rho(end))...
        %                         - h_fcn(rho(end)))/L_g_h(rho(end)); % Safety-Filter for D(t) \in [D_min_safe,D_max_safe]

        u_ctrl = @(rho) u_cancel(rho) + u_stabilize(rho) + u_constraint(rho);
        
        % additional plot parameters
        par_ctrl.c = c;
        par_ctrl.k_safety = k_safety;
end

par_ctrl.u_ctrl = u_ctrl;
par_ctrl.D_min = D_min;
par_ctrl.D_max = D_max;
par_ctrl.y_des = y_des;
par_ctrl.ctrl_mode = ctrl_mode;

%% simulate system
dynamics = @(t,rho) [(A_mat-eye(size(A_mat))*rho(end))*rho(1:end-1);
                      u_ctrl(rho)];

lambda_0 = zeros(size(A_mat,1),1); % initial conditions
lambda_0(end) = 1; % DO NOT change IC here, but in x0
rho_0 = [lambda_0;D_star];
tspan = [0 20]; % simulation horizon

[t_sample,rho_sample] = ode45(dynamics,tspan,rho_0); % run simulation

% get results
lambda_sample = rho_sample(:,1:end-1);
D_sample = rho_sample(:,end);

y_sample = C_mat*lambda_sample';

u_ctrl_sample = zeros(size(t_sample));
u_stabilize_sample = zeros(size(t_sample));
% u_cancel_sample = zeros(size(t_sample));
for kk = 1:size(lambda_sample,1)
    u_ctrl_sample(kk) = u_ctrl(rho_sample(kk,:)');
%     u_stabilize_sample(kk) = u_stabilize(rho_sample(kk,:)');
%     u_cancel_sample(kk) = u_cancel(rho_sample(kk,:)');
end

results.t_sample = t_sample;
results.rho_sample = rho_sample;
results.lambda_sample = lambda_sample;
results.u_ctrl_sample = u_ctrl_sample;
results.y_sample = y_sample;
results.D_sample = D_sample;

%% plot results - P-controller stabilizing setpoint DEBUG PLOT

if debug_plot
    plot_debug(par_sys,discretization,par_ctrl,results)
end

%% plot results - P-controller stabilizing setpoint KRSTIC plot

if standard_plot
    plot_std(par_sys,discretization,par_ctrl,results)
end

%% plot results - transformed coordinates
if transformed_plot
    plot_transformed(par_sys,discretization,par_ctrl,results)
end
%% functions

function val = u_ctrl_TC_fcn(rho, zeta_fcn, exp_zeta, A_TC, B_TC, k_TC, sigma_TC, Q_TC, y_des, C_mat)
% controller proposed by Karafyllis
D = rho(end);
% val = A_TC^2*B_TC^2 *(sigma_TC^2*Q_TC + 1)*exp_zeta(D).*(1-exp_zeta(D))./...
%     sigma_TC / Q_TC / (B_TC +1) ./ (B_TC+ exp_zeta(D))^3 - k_TC*A_TC*B_TC*exp_zeta(D)...
%     ./(B_TC+ exp_zeta(D))^2 .*(zeta_fcn(D)-sigma_TC *log(C_mat*rho(1:end-1)/y_des));
val = - k_TC*A_TC*B_TC*exp_zeta(D)...
    ./(B_TC+ exp_zeta(D))^2 .*(zeta_fcn(D)-sigma_TC *log(C_mat*rho(1:end-1)/y_des));
end

function val = eval_phi(phi,a)
% takes N-by-1 cell of functions phi and evaluates them at a, returns
% values as N-by-1 matrix
N = length(phi);
val = zeros(size(phi));
for kk = 1:N
    val(kk) = phi{kk}(a);
end
end
