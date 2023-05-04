function plot_std(par_sys,discretization,par_ctrl,results)
%plot_std

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

%% std plot
figure('units','normalized','outerposition',[0 0 1 1])
tiles_handle = tiledlayout(2,2);
title(tiles_handle,'Print Plot','Interpreter','Latex')

nexttile
plot(t_sample,u_ctrl_sample)
title('input $u(t)$')
xlabel('time t')
grid on

output_ax_handle = nexttile;
hold on
plot(t_sample,y_des(t_sample),'--k','Linewidth',1.5)
plot(t_sample,y_sample)
title('output $y(t)$')
legend('desired output $y_\mathrm{des}$','output $y(t)$','Location', 'best')
xlabel('time $t$')
grid on

nexttile
hold on
plot(t_sample,ones(size(D_sample))*D_star,'--k','Linewidth',1.5)
plot(t_sample,ones(size(D_sample))*D_min,'--b','Linewidth',1.5)
plot(t_sample,ones(size(D_sample))*D_max,'--b','Linewidth',1.5,'HandleVisibility','off')
area(t_sample(D_sample<= D_min),D_sample(D_sample<= D_min), D_min,'FaceColor',[0.8500 0.3250 0.0980],'HandleVisibility','off')
plot(t_sample,D_sample,'b')
title('dilution rate $D(t)$')
legend('steady state dilution $D^\ast$','$D_\mathrm{min}$, $D_\mathrm{max}$','dilution rate $D(t)$','Location', 'best')
xlabel('time $t$')
grid on

% plot the PDE state where x(t,a) = lambda(t)'*phi(a);

% time sample from above
% define domain sample
a_sample = 0:0.2:A;
t_plot = linspace(t_sample(1),t_sample(end),10);
sep_help = 3;
t_plot = [linspace(t_plot(1),t_plot(sep_help),20),t_plot(sep_help+1:end)];

[a_mesh,t_mesh] = meshgrid(a_sample,t_plot);
x_mesh = zeros(size(a_mesh));

for ii = 1:length(t_plot)
    for jj = 1:length(a_sample)
%         x_mesh(ii,jj) = lambda_sample(ii,:)*eval_phi(phi,a_sample(jj));
        x_mesh(ii,jj) = interp1(t_sample,lambda_sample,t_plot(ii))*eval_phi(phi,a_sample(jj));
    end
end

axes_handle = nexttile;
% surf_plot = surf(a_mesh,t_mesh,x_mesh,'FaceColor',[0 0.4470 0.7410]); %
% matlab blue
surf_plot = surf(a_mesh,t_mesh,x_mesh,'FaceColor','none');
LessEdgeSurf(surf_plot,20,10);
axes_handle.CameraPosition = [16.7896   57.3334    3.7910];
xlabel('age $a$')
ylabel('time $t$')
title('population density $f(t,a)$')

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