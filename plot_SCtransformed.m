function output = plot_SCtransformed(par_sys,discretization,par_ctrl,results)
%plot_std

%% extract input parameters:

A = par_sys.A; % max age - double
mu = par_sys.mu; % constant mortality rate - double
mu_int = par_sys.mu_int; % mortality rate integral - function
k = par_sys.k; % birth kernel - function handle
p = par_sys.p; % output kernel - function handle
gamma = par_sys.gamma; % steady-state dilution rate - double
b = par_sys.b; % SelfCompetitionKernel - fcn handle

% parameters for IC - paper [Schmidt17]
x0 = par_sys.x0; % function handle

% eigenvalues of the form lambda = -sigma/A+-j*omega/(2*pi*A)

sigma(1) = par_sys.sigma(1);
omega(1) = par_sys.omega(1);
sigma(2) = par_sys.sigma(2);
omega(2) = par_sys.omega(2);

%discretization
A_mat = discretization.A_mat;
C_mat = discretization.C_mat;
phi = discretization.phi;
b_star = discretization.b_star;
p_star = discretization.p_star;

%par_ctrl
D_ctrl = par_ctrl.D_ctrl;
D_des = par_ctrl.D_des;
% D_min = par_ctrl.D_min;
% D_max = par_ctrl.D_max;
y_des = par_ctrl.y_des;
% ctrl_mode = par_ctrl.ctrl_mode; % unneeded here

%results
t_sample = results.t_sample;
lambda_sample = results.lambda_sample;
y_sample = results.y_sample;
D_sample = results.D_sample;

%% find transformed data
% notice that the transformation includes the desired equilibrium profile,
% so find the desired equilibrium boundary value from y_des: 
f_star_0 = y_des/integral(@(a) p(a).*phi{1}(a),0,A); % c = y_des/p_star

% first, find the adjoint eigenfunction (of the zero eigenvalue of the
% differential age operator) as a lookup table
k_tilde = @(s) k(s).*exp(-mu_int(s)-gamma*s);
a_vec_lookup = linspace(0,A,20);
pi_lookup = zeros(size(a_vec_lookup));
int_k_tilde_lookup = zeros(size(a_vec_lookup));
int_k_tilde = @(a) integral(k_tilde,a,A);
for kk = 1:length(a_vec_lookup)
    a_sample = a_vec_lookup(kk);
    pi_lookup(kk) = exp(mu_int(a_sample)+gamma*a_sample).*int_k_tilde(a_sample);
    int_k_tilde_lookup(kk) = int_k_tilde(a_sample);
end

pi_fcn = @(a) interp1(a_vec_lookup,pi_lookup,a); % eigenfunction
int_k_tilde_fcn = @(a) interp1(a_vec_lookup,int_k_tilde_lookup,a); % for finding par_lambda, par_sigma
par_r = integral(@(a) k_tilde(a).*a,0,A);

% now find pi_vec
f_star_fcn = @(a) f_star_0*phi{1}(a); % desired equilibrium profile

pi_vec_denominator = integral(@(a)pi_fcn(a).*f_star_fcn(a),0,A);

pi_vec = zeros(size(phi));
for kk = 1:length(pi_vec)
    pi_vec(kk) = integral(@(a) pi_fcn(a).*phi{kk}(a),0,A)/pi_vec_denominator;
end

% extract transformed states
eta_sample = log(pi_vec'*lambda_sample');
psi_sample_posTime = (eval_phi(phi,0)'*lambda_sample')./(f_star_0*pi_vec'*lambda_sample')-1;
t_sample_negTime = -2:.01:0;
psi_sample_negTime = phi{end}(-t_sample_negTime)./f_star_fcn(-t_sample_negTime)/pi_vec(end) -1;
t_sample_ext = [t_sample_negTime(1:end-1),t_sample']; % extended time sample
psi_sample_ext = [psi_sample_negTime(1:end-1),psi_sample_posTime]; % extended psi sample

% finding par_lambda
lambda_vec = 0:.05:2;
lambda_res_vec = zeros(size(lambda_vec)); % residual for finding par_lambda

for kk = 1:length(lambda_vec)
    lambda_res_vec(kk) = 1 - integral(@(a) abs(k_tilde(a) - par_r*lambda_vec(kk)*int_k_tilde_fcn(a)),0,A);
end

[max_res_lambda,lambda_indx] = max(lambda_res_vec);
par_lambda = lambda_vec(lambda_indx);

if par_lambda <= 0 || max_res_lambda <= 0
    error('Error finding par_lambda.')
end

% finding par_sigma
sigma_vec = 0:.01:2;
sigma_res_vec = zeros(size(sigma_vec)); % residual for finding par_sigma

for kk = 1:length(sigma_vec)
    sigma_res_vec(kk) = 1 - integral(@(a) abs(k_tilde(a) - par_r*par_lambda*int_k_tilde_fcn(a))...
        .*exp(sigma_vec(kk)*a),0,A);
end

sigma_admissable = sigma_vec(sigma_res_vec>0);
par_sigma = sigma_admissable(end); % assuming sigma_res strictly decreasing in sigma

error_help = 1 - integral(@(a) abs(k_tilde(a) - par_r*par_lambda*int_k_tilde_fcn(a))...
        .*exp(par_sigma*a),0,A);
if error_help <= 0
    error('Error finding par_sigma.')
end

% Lyapunov Functional
% par_sigma = .05; % suff small parameter
par_b_1 = 10*b_star*exp(2*par_sigma*A)/par_sigma; % suff large parameter
G_Lyap_Sample = zeros(size(t_sample)); % Lyap Functional G wrt psi
g_1_Sample = zeros(size(t_sample)); % non-decreasing functional
g_2_Sample = zeros(size(t_sample)); % non-increasing functional
for kk = length(t_sample_negTime)+1:length(t_sample_ext)
    indx_start = find(t_sample_ext>=t_sample_ext(kk)-A,1);
    psi_stage = psi_sample_ext(indx_start:kk);
    t_sample_stage = t_sample_ext(indx_start:kk);
    G_Lyap_Sample(kk-length(t_sample_negTime)) = ...
     max(abs(psi_stage).*exp(par_sigma*(t_sample_stage-t_sample_ext(kk))));
    g_1_Sample(kk-length(t_sample_negTime)) = min(psi_stage);
    g_2_Sample(kk-length(t_sample_negTime)) = max(psi_stage);
end
V_Lyap_Sample = .5*(1-exp(-eta_sample)).^2+.5*par_b_1*G_Lyap_Sample'.^2;

% % (quasistatic) active filter set.
% eta_ASF_0 = (c+1-k_safety)/c*D_sample - (c+1)*D_star/c;

% find functional v2(psi_t)
v2_psi_vec = zeros(size(t_sample));

for kk = 1:length(t_sample)
    a_vec_lookup = t_sample_ext(1:length(t_sample_negTime)+kk-1)...
                    - t_sample(kk);
    psi_lookup = psi_sample_ext(1:length(t_sample_negTime)+kk-1);
    psi_stage = @(a) interp1(a_vec_lookup,psi_lookup,a);
    v2_psi_vec(kk)=integral(@(a)b(a).*f_star_fcn(a).*psi_stage(-a),0,A);
end

%% plotting the (evil) functional v2 of psi over time...
figure
plot(t_sample,v2_psi_vec)
xlim([0 10])
xlabel('time $t$')
ylabel('$v_2(\psi_t)$')
grid on

%% plotting the lyapunov function of psi over time
% figure
% plot(t_sample,G_Lyap_Sample)
% hold on
% plot(t_sample,G_Lyap_Sample(1)*exp(-par_sigma*t_sample))

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

traj_plot = plot(eta_sample,psi_sample_posTime);
setpoint_pl = plot(0, 0, 'k.','MarkerSize', 20);
legend('trajectory $(\psi_t,\eta(t))$',...
        'Location', 'best')
title('phase portrait projected to $\psi$-$\eta$ plane')
ylabel('internal state $\psi$')
xlabel('1-dim. state $\eta$')
grid on

% Lyap functional plot
nexttile
hold on
plot(t_sample,V_Lyap_Sample)
plot(t_sample,V_Lyap_Sample(1)*exp(-min(par_sigma,b_star*f_star_0)*t_sample))
title('Lyapunov Functional $V(\eta(t),\delta(t),\psi_t)$')
legend('Lyapunov Functional $V(\eta(t),\delta(t),\psi_t)$',...
        'exponential bound $V_0 e^{-\mu t}$')
xlabel('time $t$')
grid on

% OUTPUT - data

output.eta_sample = eta_sample;
output.C_Lyap_Sample = V_Lyap_Sample;
output.t_sample_ext = t_sample_ext;
output.psi_sample_ext = psi_sample_ext;

%% plot of new lyap-terms

% figure('units','normalized','outerposition',[0 0 1 1])
% tiles_handle = tiledlayout(2,2);
% title(tiles_handle,'intuition of the functional $1-\frac{1}{\Pi(f_\rho)}$','Interpreter','Latex')
% 
% % nexttile
% % hold on
% % plot(t_sample,(1-exp(-eta_sample)))
% 
% % with a basis function
% nexttile 
% chebychev = @(t,a) cos(t.*acos(1-a/A)).^2;
% t_vec = 0:.1:20;
% a_vec = 0:.01:A;
% [t_mesh,a_mesh] = meshgrid(t_vec,a_vec);
% basis_fcn = @(t,a) f_star_fcn(a) + (1 - chebychev(t,a)).*exp(1-a/A);
% Z = basis_fcn(t_mesh,a_mesh);
% s1_handle = surf(t_mesh,a_mesh,Z);
% % s1_handle.FaceColor = 'none';
% LessEdgeSurf(s1_handle,20,20)
% s1_handle.EdgeColor = 'none';
% xlabel('family parameter $\rho$')
% ylabel('age $a$')
% zlabel('population density $f_\rho(a)$')
% 
% % %debug
% % test_vec = 0:.01:1;
% % figure
% % plot(test_vec, cos(2*acos(test_vec)))
% % hold on
% % plot(test_vec, chebychev(2,test_vec*A))
% 
% N_t = size(t_mesh,2);
% basisfcn_lyap = zeros(N_t,1);
% 
% for kk = 1:N_t
%     basisfcn_lyap(kk) = 1-integral(@(a) pi_fcn(a).*basis_fcn(t_vec(kk),a),0,A)^-1 ...
%     * pi_vec_denominator;
% end
% 
% nexttile
% plot(t_vec,basisfcn_lyap)
% xlabel('family parameter $\rho$')
% ylabel('functional $1-\frac{1}{\Pi(f_\rho)}$')
% grid on
% 
% % another basis function
% nexttile 
% t_vec = 0:.1:20;
% a_vec = 0:.01:A;
% [t_mesh,a_mesh] = meshgrid(t_vec,a_vec);
% basis_fcn = @(t,a) f_star_fcn(a).*(.5+t).^2;
% Z = basis_fcn(t_mesh,a_mesh);
% s1_handle = surf(t_mesh,a_mesh,Z);
% % s1_handle.FaceColor = 'none';
% LessEdgeSurf(s1_handle,20,20)
% s1_handle.EdgeColor = 'none';
% xlabel('family parameter $\rho$')
% ylabel('age $a$')
% zlabel('population density $f_\rho(a)$')
% 
% % %debug
% % test_vec = 0:.01:1;
% % figure
% % plot(test_vec, cos(2*acos(test_vec)))
% % hold on
% % plot(test_vec, chebychev(2,test_vec*A))
% 
% N_t = size(t_mesh,2);
% basisfcn_lyap = zeros(N_t,1);
% 
% for kk = 1:N_t
%     basisfcn_lyap(kk) = 1-integral(@(a) pi_fcn(a).*basis_fcn(t_vec(kk),a),0,A)^-1 ...
%     * pi_vec_denominator;
% end
% 
% nexttile
% plot(t_vec,basisfcn_lyap)
% xlabel('family parameter $\rho$')
% ylabel('functional $1-\frac{1}{\Pi(f_\rho)}$')
% grid on

%% plot of transformed states - KRSTIC plot

figure('units','normalized','outerposition',[0 0 1 1])
tiles_handle = tiledlayout(2,2);
title(tiles_handle,'functionals $(\zeta,\psi)$','Interpreter','Latex')

% nexttile
% hold on
% plot(t_sample,(1-exp(-eta_sample)))

% with a basis function
rho_par_vec = 0:.01:10; %4:.0001:4.2; %0.85:.0001:1; %
a_vec = 0:.01:2;
chebychev = @(rho,a) cos(rho.*acos(a/a_vec(end)));
[rho_par_mesh,a_mesh] = meshgrid(rho_par_vec,a_vec);
% basis_fcn = @(t,a) f_star_fcn(a) + (1 - chebychev(t,a)).*exp(1-a/A);
basis_fcn = @(rho,a) f_star_fcn(a).*(1+chebychev(rho,a));
Z = basis_fcn(rho_par_mesh,a_mesh);

% plot basis function
nexttile 
s1_handle = surf(rho_par_mesh,a_mesh,Z);
% s1_handle.FaceColor = 'none';
LessEdgeSurf(s1_handle,20,20)
s1_handle.EdgeColor = 'none';
xlabel('family parameter $\rho$')
ylabel('age $a$')
zlabel('population density $f_\rho(a)$')

% %debug
% test_vec = 0:.01:1;
% figure
% plot(test_vec, cos(2*acos(test_vec)))
% hold on
% plot(test_vec, chebychev(2,test_vec*A))

N_rho = size(rho_par_mesh,2);
basisfcn_Pi_c = zeros(N_rho,1);

for kk = 1:N_rho
    basisfcn_Pi_c(kk) = integral(@(a) pi_fcn(a).*basis_fcn(rho_par_vec(kk),a),0,A) ...
    / pi_vec_denominator;
end

basisfcn_zeta = 1-basisfcn_Pi_c.^(-1);

% finding psi
nexttile

basisfcn_Pi_c_fcn = @(rho) interp1(rho_par_vec,basisfcn_Pi_c,rho);

Psi_mesh = basis_fcn(rho_par_mesh,a_mesh)./f_star_fcn(a_mesh)...
            ./basisfcn_Pi_c_fcn(rho_par_mesh)-1;
s1_handle = surf(rho_par_mesh,a_mesh,Psi_mesh);
% s1_handle.FaceColor = 'none';
LessEdgeSurf(s1_handle,20,20)
s1_handle.EdgeColor = 'none';
xlabel('family parameter $\rho$')
ylabel('age $a$')
zlabel('functional $\psi = \frac{f_\rho(a)}{f^\ast_c(a) \Pi_c(f_\rho)}-1$')

% plot zeta-Pi
nexttile
Pi_vec = .01:.02:5;
zeta_vec = 1-1./Pi_vec;
plot(Pi_vec,zeta_vec)
xlabel('functional $\Pi$')
ylabel('functional $\zeta = 1-\Pi^{-1}$')
ylim([-10,2])
grid on

% plot zeta-rho
nexttile
plot(rho_par_vec,basisfcn_zeta)
xlabel('family parameter $\rho$')
ylabel('functional $\zeta = 1-\Pi(f_\rho)^{-1}$')
grid on
% %debug
% test_vec = 0:.01:1;
% figure
% plot(test_vec, cos(2*acos(test_vec)))
% hold on
% plot(test_vec, chebychev(2,test_vec*A))


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