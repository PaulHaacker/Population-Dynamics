close all
clear
%% Population Dynamics with Actuator Dynamics
% this script extends the Population Dynamics system with integrator
% dynamics for the dilution rate u(t)
%       \dot u(t) = w(t),
% where w(t) is the new input. Simulate the new dynamics with the
% backstepping-based controller.


%% finding the matrices for the Galerkin-based simulation

% choose the basis phi(a) = [zero eigenfunction; two first eigenfunction
% pairs and function for the IC]
% --> requires two pairs of nonzero eigenvalues, computed by [Schmidt17]


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
sigma(1) = 4.0335;
omega(1) = 55.4606;
sigma(2) = 4.9866;
omega(2) = 95.7048;

N_EV = 2; % number of nonzero eigenvalues considered

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

%% ------ basis of trial functions

phi = cell(6,1);
phi{1} = @(a) exp(-(u_star+mu).*a);
for kk = 1:N_EV
    phi{2*kk} = @(a) sin(omega(kk).*a/(2*pi*A)).*exp(-sigma(kk).*a/A).*phi{1}(a);
    phi{2*kk+1} = @(a) cos(omega(kk).*a/(2*pi*A)).*exp(-sigma(kk).*a/A).*phi{1}(a);
end
phi{2*N_EV+2} = @(a) x0(a); % initial condition fcn

% apply differential operator to basis:

D_phi = cell(6,1);
D_phi{1} = @(a) 0;

Lambda_mat = @(sigma,omega) [   sigma/A,    omega/2/pi/A;
                            -omega/2/pi/A,  sigma/A];

for kk = 1:N_EV
    Lambda_k = Lambda_mat(sigma(kk),omega(kk));
    
    D_phi{2*kk} = @(a) Lambda_k(1,1)*phi{2*kk}(a) + Lambda_k(1,2)*phi{2*kk+1}(a);
    D_phi{2*kk+1} = @(a) Lambda_k(2,1)*phi{2*kk}(a) + Lambda_k(2,2)*phi{2*kk+1}(a);
end


% differential operator applied to IC
syms a_sym
x0_sym = x0(a_sym);
D_x0_sym = diff(x0_sym) + (mu+u_star)*x0_sym;

D_phi{2*N_EV+2} = matlabFunction(D_x0_sym);

%% Plot Basis functions - works, but with dirty tweaks...
% differences:
% - needed to flip the sign of the eigenvalues real part sigma
%
% Q: Does the sign of an eigenfunction matter?
% A: let f(a) satisfy the ODE 
%           Df(a) = \lambda f(a)
% where D is a differential operator. Let 
%           g(a) := - f(a).
% Whenever D is a linear operator, g also satisfies the above ODE and thus
% even the function with the sign flipped is an eigenfunction..
%  --> This is verified by switching the signs in the sine-argument, or
%  equivalently, flipping the signs of the whole eigenfucntions.

% a_sample = linspace(0,A,100);
% 
% figure
% plot(a_sample,phi{1}(a_sample))
% grid on
% hold on
% 
% for kk = 1:N_EV
%     plot(a_sample,phi{2*kk}(a_sample))
%     plot(a_sample,phi{2*kk+1}(a_sample))
% end
% plot(a_sample,x0(a_sample))
% 
% legend('$\varphi_0$','$\varphi_{1,1}$','$\varphi_{1,2}$','$\varphi_{2,1}$',...
%     '$\varphi_{2,2}$','$x_0$','Interpreter','latex','FontSize',12)

%% ------ find matrices

N = length(phi);
Phi_1 = zeros(N);
Phi_2 = zeros(N);

for ii = 1:N
    for jj = 1:N
        Phi_1(ii,jj) = integral(@(a) phi{ii}(a).* phi{jj}(a),0,A);
        Phi_2(ii,jj) = integral(@(a) phi{ii}(a).* D_phi{jj}(a),0,A);
    end
end

% system matrix
A_mat = eye(N)*u_star - Phi_1\Phi_2;

% output matrix, where y(t) = C*lambda(t)
C = zeros(size(phi))';
for ii = 1: length(phi)
   C(ii) = p*integral(@(a) phi{ii}(a),0,A);
end

% output value at equilibrium lambda = [1 0 ... 0]
y_eq = C(1);

%% simulate linear system - P-controller stabilizing setpoint
% here, with controller w(t) = w_cancel(t) + w_stabilize(t)
% notice that y(t) = C*lambda(t)
% denote the simulation state by rho(t) = [lambda'(t),u(t)]

% choose desired setpoint for output - equivalent to choosing a desired
% equilibrium profile x^\ast(a), or better its family parameter.
% y_des = 1.5;
y_des = .5;

% --- define controller - backstepping type
c = 1; % control gain c > 0

w_cancel = @(rho) -rho(end)-p/(C*rho(1:end-1))...
            *(rho(1:end-1)'*(eval_phi(phi,A)-eval_phi(phi,0)))-mu;
w_stabilize = @(rho) -c*(rho(end)-u_star-log(C*rho(1:end-1)/y_des));
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

y_sample = C*lambda_sample';

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
