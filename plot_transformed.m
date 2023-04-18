function output = plot_transformed(par_sys,discretization,par_ctrl,results)
%plot_std

%% extract input parameters:

% par_sys
A = par_sys.A; % max age - double
mu = par_sys.mu; % constant mortality rate - double
mu_int = par_sys.mu_int; % mortality rate integral - function
k = par_sys.k; % birth kernel - function handle
p = par_sys.p; % output kernel - function handle
D_star = par_sys.D_star; % steady-state dilution rate - double

x0 = par_sys.x0; % function handle

sigma(1) = par_sys.sigma(1); % eigenvalues
omega(1) = par_sys.omega(1);
sigma(2) = par_sys.sigma(2);
omega(2) = par_sys.omega(2);

%discretization
A_mat = discretization.A_mat;
C_mat = discretization.C_mat;
phi = discretization.phi;

%par_ctrl
u_ctrl = par_ctrl.u_ctrl;
D_min = par_ctrl.D_min;
D_max = par_ctrl.D_max;
y_des = par_ctrl.y_des;
y_des_d = par_ctrl.y_des_d;
k_nom = par_ctrl.k_nom;
ctrl_mode = par_ctrl.ctrl_mode;

%results
t_sample = results.t_sample;
rho_sample = results.rho_sample;
lambda_sample = results.lambda_sample;
u_ctrl_sample = results.u_ctrl_sample;
y_sample = results.y_sample;
D_sample = results.D_sample;

%% find transformed data
% notice that the transformation includes the desired equilibrium profile,
% so find the desired equilibrium boundary value from y_des: 
p_star = integral(@(a) p(a).*phi{1}(a),0,A);
% f_star_0 = y_des/integral(@(a) p(a).*phi{1}(a),0,A);

% first, find the adjoint eigenfunction (of the zero eigenvalue of the
% differential age operator) as a lookup table
integrand_pi = @(s) k(s).*exp(-mu_int(s)-D_star*s);
a_vec_lookup = linspace(0,A,20);
pi_lookup = zeros(size(a_vec_lookup));
integral_pi = @(a) integral(integrand_pi,a,A);
for kk = 1:length(a_vec_lookup)
    a_sample = a_vec_lookup(kk);
    pi_lookup(kk) = exp(mu_int(a_sample)+D_star*a_sample).*integral_pi(a_sample);
end

pi_fcn = @(a) interp1(a_vec_lookup,pi_lookup,a);

% now find pi_vec
% f_star_fcn = @(a) f_star_0*phi{1}(a); % desired equilibrium profile
gamma_fcn = @(a) phi{1}(a)/p_star; % quasi-static equilibrium profile

pi_vec_denominator = integral(@(a)pi_fcn(a).*gamma_fcn(a),0,A);

pi_vec = zeros(size(phi));
for kk = 1:length(pi_vec)
    pi_vec(kk) = integral(@(a) pi_fcn(a).*phi{kk}(a),0,A)/pi_vec_denominator;
end

% extract transformed states
eta_sample = log(pi_vec'*lambda_sample'./y_des(t_sample)');
psi_sample_posTime = (eval_phi(phi,0)'*lambda_sample')./(pi_vec'*lambda_sample'/p_star)-1;
t_sample_negTime = -2:.01:0;
psi_sample_negTime = phi{end}(-t_sample_negTime)./gamma_fcn(-t_sample_negTime)/pi_vec(end) -1;

% Lyapunov Functional
par_sigma = .1; % suff small parameter
par_M_hat = 2*exp(2*par_sigma*A)/par_sigma; % suff large parameter
delta_sample = D_sample' - D_star - k_nom* log(y_sample./y_des(t_sample)')...
                + y_des_d(t_sample)'/y_des(t_sample)'; % dilution error
G_Lyap_Sample = zeros(size(t_sample)); % Lyap Functional G wrt psi
t_sample_ext = [t_sample_negTime, t_sample']; % extended time sample
psi_sample_ext = [psi_sample_negTime,psi_sample_posTime]; % extended psi sample
g_1_Sample = zeros(size(t_sample)); % non-decreasing functional
g_2_Sample = zeros(size(t_sample)); % non-increasing functional
for kk = length(t_sample_negTime)+1:length(t_sample_ext)
    indx_start = find(t_sample_ext>=t_sample_ext(kk)-A,1);
    psi_stage = psi_sample_ext(indx_start:kk);
    t_sample_stage = t_sample_ext(indx_start:kk);
    G_num = max(abs(psi_stage).*exp(par_sigma*(t_sample_stage-t_sample_ext(kk))));
    G_den = 1 + min(0,min(psi_stage));
    G_Lyap_Sample(kk-length(t_sample_negTime)) = G_num/G_den;
    g_1_Sample(kk-length(t_sample_negTime)) = min(psi_stage);
    g_2_Sample(kk-length(t_sample_negTime)) = max(psi_stage);
end
C_Lyap_Sample = .5*eta_sample.^2+.5*delta_sample.^2+.5*par_M_hat*G_Lyap_Sample'.^2;

% % (quasistatic) active filter set.
% eta_ASF_0 = (c+1-k_safety)/c*D_sample - (c+1)*D_star/c;

%% plot transformed data
% plotting
figure('units','normalized','outerposition',[0 0 1 1])
tiles_handle = tiledlayout(2,2);
title(tiles_handle,'transformed states','Interpreter','Latex')

nexttile
plot(t_sample,eta_sample)
title('1-dim. state $\eta(t)$')
xlim([-A,t_sample(end)])
xlabel('time $t$')
grid on

psi_plots_handle = nexttile;
hold on
plot_handle = plot(t_sample_ext,psi_sample_ext);
plot(t_sample,g_1_Sample)
plot(t_sample,g_2_Sample)
xlabel('time $t$')
legend('infinite-dim. state $\psi(t)$',...
    'functional $g_1(\psi_t)= \min_{a\in[0,A]}\psi(t-a)$ non-decreasing',...
    'functional $g_2(\psi_t)= \max_{a\in[0,A]}\psi(t-a)$ non-increasing',...
    'Location','Southeast')
uistack(plot_handle,'top')
title('infinite-dim. state $\psi(t)$')
xlabel('time $t$')
grid on

% phase portrait plot
nexttile
hold on
switch ctrl_mode
    case 'Backstepping'
        % additional parameters
        c = par_ctrl.c;
        k_safety = par_ctrl.k_safety;
        
        % plot
        traj_plot = plot(D_sample,eta_sample);
        setpoint_pl = plot(D_star, 0, 'k.','MarkerSize', 20);
        D_lim = xlim;
        eta_ASF_0 = (c+1-k_safety)/c*D_lim - (c+1)*D_star/c;
        plot(D_lim,eta_ASF_0,'r','HandleVisibility','off')
        help1 = ylim;
        area_plot = area(D_lim,eta_ASF_0,help1(1),'FaceColor','#ffcccb');
        area_plot.EdgeColor = 'none';
        area_plot.FaceAlpha = .5;
        legend('trajectory $(D,\eta)(t)$',...
            'setpoint $(D,\eta,\psi) = (D^\star,0,0)$',...
            'quasistatic active safety filter set $\mathcal{X}_{\mathrm{ASF},\psi = 0}$',...
            'Location', 'best')
        uistack(traj_plot,'top')
        uistack(setpoint_pl,'top')
        title('phase portrait projected to $D$-$\eta$ plane')
        xlabel('Dilution rate $D$')
        ylabel('1-dim. state $\eta$')
        grid on
    case 'Karafyllis'
        traj_plot = plot(D_sample,eta_sample);
        setpoint_pl = plot(D_star, 0, 'k.','MarkerSize', 20);
        legend('trajectory $(D,\eta)(t)$',...
            'setpoint $(D,\eta,\psi) = (D^\star,0,0)$',...
            'Location', 'best')
        uistack(traj_plot,'top')
        uistack(setpoint_pl,'top')
        title('phase portrait projected to $D$-$\eta$ plane')
        xlabel('Dilution rate $D$')
        ylabel('1-dim. state $\eta$')
        grid on
end


% Lyap functional plot
nexttile
hold on
plot(t_sample,C_Lyap_Sample)
title('Lyapunov Functional $\tilde C(\eta(t),\delta(t),\psi_t)$')
xlabel('time $t$')
grid on

% OUTPUT - data

output.eta_sample = eta_sample;
output.C_Lyap_Sample = C_Lyap_Sample;
output.t_sample_ext = t_sample_ext;
output.psi_sample_ext = psi_sample_ext;

%% quick and dirty: plot of local simultaneous stability and safety
% figure
% hold on
% c = par_ctrl.c;
% k_safety = par_ctrl.k_safety;
% 
% % plot
% traj_plot = plot(D_sample,eta_sample);
% setpoint_pl = plot(D_star, 0, 'k.','MarkerSize', 20);
% D_lim = xlim;
% eta_ASF_0 = (c+1-k_safety)/c*D_lim - (c+1)*D_star/c;
% plot(D_lim,eta_ASF_0,'r','HandleVisibility','off')
% help1 = ylim;
% area_plot = area(D_lim,eta_ASF_0,help1(1),'FaceColor','#ffcccb');
% area_plot.EdgeColor = 'none';
% area_plot.FaceAlpha = .5;
% legend('trajectory $(D,\eta)(t)$',...
%     'setpoint $(D,\eta,\psi) = (D^\star,0,0)$',...
%     'quasistatic active safety filter set $\mathcal{X}_{\mathrm{ASF},\psi = 0}$',...
%     'Location', 'best')
% uistack(traj_plot,'top')
% uistack(setpoint_pl,'top')
% title('phase portrait projected to $D$-$\eta$ plane')
% xlabel('Dilution rate $D$')
% ylabel('1-dim. state $\eta$')
% grid on

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