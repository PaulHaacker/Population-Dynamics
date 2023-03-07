function plot_debug(par_sys,discretization,par_ctrl,results)
%plot_debug 

%% extract input parameters:

%par_sys
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

%results
t_sample = results.t_sample;
rho_sample = results.rho_sample;
lambda_sample = results.lambda_sample;
u_ctrl_sample = results.u_ctrl_sample;
y_sample = results.y_sample;
D_sample = results.D_sample;

%% debug plot
figure('units','normalized','outerposition',[0 0 1 1])
tiles_handle = tiledlayout(2,2);
title(tiles_handle,'Debug Plot','Interpreter','Latex')

nexttile
plot(t_sample,lambda_sample)
title('discretized states $\lambda(t)$ - backstepping controller')
legend('1','2','3','4','5','6')
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
% str = {'steady-state error $\Delta y = y_\mathrm{ss} - y_\mathrm{des}$ = '...
%     ,num2str(y_sample(end)-y_des)};
% text(max(xlim), min(ylim),str, 'Horiz','right', 'Vert','top')

% REMARK: notation in 
% - documentation is Dilution rate D(t), voltage input u(t)

nexttile
hold on
plot(t_sample,ones(size(D_sample))*D_star,'--k','Linewidth',1.5)
plot(t_sample,D_sample)
plot(t_sample,u_ctrl_sample)
% plot(t_sample,u_cancel_sample)
% plot(t_sample,u_stabilize_sample)
title('dilution rate $D(t)$ and input $u(t)$ - backstepping controller')
legend('steady state dilution $D^\ast$','dilution rate $D(t)$',...
    'input $u(t) = u_\mathrm c(t) +u_\mathrm s(t)$'...
    ); % ,'cancelling terms $u_\mathrm c(t)$','stabilizing terms $u_\mathrm s(t)$')
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
LessEdgeSurf(surf_plot,20,20);
axes_handle.CameraPosition = [15.7853   91.8902    2.8718];

xlabel('age $a$')
ylabel('time $t$')
title('population density $x(t,a)$ - backstepping controller')
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