close all
clear

%% ------ parameters

A = 2; % max age
mu = .1; % constant mortality rate
k = @(a) 2*a.*(A-a); % birth kernel
p = 1; % output kernel

% %% ------ find equilibrium dilution rate numerically w/ analytic solution
% % of integral
% 
% res_bc = @(c) 2./c.^3.*((A*c+2).*exp(-A*c)+A*c-2)-1; % analytic expression for 
%  % residual of boundary condition with the above parameters where c =
%  % mu+u_star
% 
% c_0 = 1.5; % initial guess
% c_star = fzero(res_bc,c_0);
% u_star = c_star - mu; % equilibrium dilution rate
 
%% ------ find equilibrium dilution rate numerically w/o analytic solution
% of integral - provides a higher level of automation

% integrand of the Lotka-Sharpe condition
res_LS_int = @(a,u_star) k(a).*exp(-(u_star+mu).*a);

% Lotka-Sharpe in residual form as fcn of u_star
res_LotkaSharpe_fcn = @(u_star) integral(@(a) res_LS_int(a,u_star),0,A)-1;

% verify the numerical integration by plugging in the known value of u_star
u_star_test = .9998;
res_LS_test = res_LotkaSharpe_fcn(u_star_test);

u_0 = .5; % initial guess
u_star_num = fzero(res_LotkaSharpe_fcn,u_0);

%% plot Lotka-Sharpe residual fcn wrt. u_star

% u_star_sample = 0:.1:10;
% res_LotkaSharpe_sample = zeros(size(u_star_sample));
% for ii = 1: length(u_star_sample)
%     res_LotkaSharpe_sample(ii) = res_LotkaSharpe_fcn(u_star_sample(ii));
% end
% 
% 
% figure
% plot(u_star_sample,res_LotkaSharpe_sample)
% title('Lotka-Sharpe condition')
% xlabel('$u^\star$','Interpreter','Latex')
% xlabel('residual')
% grid on

%% ----- define eigenfunction of zero eigenvalue

phi_0 = @(a) exp(-(u_star_num+mu)*a); % eigenfunction of steady-state operator 
 % with eigenvalue = 0

%% ----- finding eigenvalues using newton-iteration - convergence to nonzero roots unsuccessful...
% ---> failed to find initial conditions with pos. def. jacobians

% % integrands of the residuals of the implicit(to-be solved) eigenvalue
% % equations.
% res_1_eig_int = @(a,sigma,omega) k(a).*cos(omega.*a/2/pi/A).*exp(-sigma.*a/A).*phi_0(a);
% res_2_eig_int = @(a,sigma,omega) k(a).*sin(omega.*a/2/pi/A).*exp(-sigma.*a/A).*phi_0(a);
% 
% % define 2-dim residual - vec(1) = sigma, vec(2) = omega.
% res_ev = @(vec) [integral(@(a) res_1_eig_int(a,vec(1),vec(2)),0,A)-1;
%                  integral(@(a) res_2_eig_int(a,vec(1),vec(2)),0,A)];   
%              
% % integrands of partial derivatives of residual wrt. (sigma,omega)
% dsigma_R1 = @(a,sigma,omega) -k(a).*cos(omega*a/2/pi/A).*sigma.*exp(-sigma*a/A).*phi_0(a);
% dsigma_R2 = @(a,sigma,omega) -k(a).*sin(omega*a/2/pi/A).*sigma.*exp(-sigma*a/A).*phi_0(a);
% domega_R1 = @(a,sigma,omega) -k(a).*sin(omega*a/2/pi/A).*a/(2*pi*A).*exp(-sigma*a/A).*phi_0(a);
% domega_R2 = @(a,sigma,omega) k(a).*cos(omega*a/2/pi/A).*a/(2*pi*A).*exp(-sigma*a/A).*phi_0(a);
% 
% % jacobian of residual wrt. (sigma, omega)
% resJac_ev = @(vec) [integral(@(a) dsigma_R1(a,vec(1),vec(2)),0,A), integral(@(a) domega_R1(a,vec(1),vec(2)),0,A);
%                     integral(@(a) dsigma_R2(a,vec(1),vec(2)),0,A), integral(@(a) domega_R2(a,vec(1),vec(2)),0,A);];   

%%  run the newton algorithm

% steps = 1000;
% x_01 = [-3.6;56]; % initial guess for k = 1 from surf plot below
% tol = .001;
% [y_1,yfinal_1,act_steps_1,res_history_1] = NewtonIteration(res_ev, resJac_ev, x_01, tol, steps);
% 
% figure
% plot(y_1','.-')
% grid on
% title('convergence of newton-iteration run01 - root variables')
% figure
% plot(res_history_1','.-')
% grid on
% title('convergence of newton-iteration run01 - residuals')
% 
% steps = 1000;
% x_02 = [-4.2;98]; % initial guess for k = 2 from surf plot below
% tol = .001;
% [y_2,yfinal_2,act_steps_2,res_history_2] = NewtonIteration(res_ev, resJac_ev, x_01, tol, steps);
% 
% figure
% plot(y_2','.-')
% grid on
% title('convergence of newton-iteration run02 - root variables')
% figure
% plot(res_history_2','.-')
% grid on
% title('convergence of newton-iteration run02 - residuals')

%% ----- finding eigenvalues using newton-iteration - using symbolic toolbox

% integrands of the residuals of the implicit(to-be solved) eigenvalue
% equations.
res_1_eig_int = @(a,sigma,omega) k(a).*cos(omega.*a/2/pi/A).*exp(-sigma.*a/A).*phi_0(a);
res_2_eig_int = @(a,sigma,omega) k(a).*sin(omega.*a/2/pi/A).*exp(-sigma.*a/A).*phi_0(a);

% define symbolic fcns
syms a sigma omega
res_1_eig_int_sym(a,sigma,omega) = res_1_eig_int(a,sigma,omega);
res_2_eig_int_sym(a,sigma,omega) = res_2_eig_int(a,sigma,omega);

% compute symbolic integral
res_ev_sym = [int(res_1_eig_int_sym,a,[0 A]) - 1;
             int(res_2_eig_int_sym,a,[0 A])];

% find symbolic jacobian
resJac_ev_sym = jacobian(res_ev_sym,[sigma, omega]);

% convert into functions compatible with below newton fcn
help = matlabFunction(res_ev_sym);
res_ev_fcn = @(vec) help(vec(1),vec(2));
help = matlabFunction(resJac_ev_sym);
resJac_ev_fcn = @(vec) help(vec(1),vec(2));

%%  run the newton algorithm

steps = 1000; % max number of steps
x_01 = [-3.6;56]; % initial guess for k = 1 from surf plot below
tol = .001; % absolute errror tolerance
[y_1,yfinal_1,act_steps_1,res_history_1] = ...
    NewtonIteration(res_ev_fcn, resJac_ev_fcn,...
    x_01, tol, steps);

figure
plot(y_1','.-')
grid on
title('convergence of newton-iteration run01 - root variables')
figure
plot(res_history_1','.-')
grid on
title('convergence of newton-iteration run01 - residuals')

% steps = 1000;
% x_02 = [-4.2;98]; % initial guess for k = 2 from surf plot below
% tol = .001;
% [y_2,yfinal_2,act_steps_2,res_history_2] = NewtonIteration(res_ev, resJac_ev, x_01, tol, steps);
% 
% figure
% plot(y_2','.-')
% grid on
% title('convergence of newton-iteration run02 - root variables')
% figure
% plot(res_history_2','.-')
% grid on
% title('convergence of newton-iteration run02 - residuals')


%% test the residuals with the known solution

% % values suggested by paper in [Schmidt17] eq. (49)-(50) ... residual
% % test FAILED
% sigma_test = 1.250;1.010;
% omega_test = 0.601;.351;

% % values suggested by paper in [Schmidt17] figure 2 ... residual test
% % FAILED
% sigma_test = 4.1;
% omega_test = 37.69;
% % sigma_test = 5;
% % omega_test = 62.83;

% values I suspect by paper in [Schmidt17] eq. (49)-(50)
% -> they match the plots exactly, with a flipped sign of argument in the
% sine - residual test successful !
% sigma_test = -4.04;
% omega_test = 55.43;
sigma_test = -5;
omega_test = 94.92;

res_test = res_ev_fcn([sigma_test, omega_test]);

%% plot residuals surface over meshgrid


% define 2-dim residual - vec(1) = sigma, vec(2) = omega.
res_ev = @(vec) [integral(@(a) res_1_eig_int(a,vec(1),vec(2)),0,A)-1;
                 integral(@(a) res_2_eig_int(a,vec(1),vec(2)),0,A)];    

sigma_sample = -5:.2:5;
omega_sample = 1:2:100;

% sigma varies with columns, omega with rows:
[sigma_grid,omega_grid] = meshgrid(sigma_sample,omega_sample);

res_norm_grid = zeros(size(sigma_grid));

for ii = 1: length(omega_sample)
    for jj = 1: length(sigma_sample)
        res_norm_grid(ii,jj) = norm(res_ev([sigma_sample(jj), omega_sample(ii)]));
    end
end

% add colormap to indicate local minima
% colormap_grid = res_norm_grid < .75;

figure
% surface1 = surf(sigma_grid,omega_grid,res_norm_grid,double(colormap_grid));
surface1 = surf(sigma_grid,omega_grid,log(res_norm_grid));
xlabel('$\sigma$','Interpreter','Latex')
ylabel('$\omega$','Interpreter','Latex')
zlabel('logarithmic norm of residual')
LessEdgeSurf(surface1)
c = parula;
c = flipud(c);
colormap(c);

figure
contour1 = contour(sigma_grid,omega_grid,res_norm_grid,50);
xlabel('$\sigma$','Interpreter','Latex')
ylabel('$\omega$','Interpreter','Latex')
title('Contour plot of norm of residual')


%% ----- finding eigenvalues using fzero - succeeds to converge to nonzero EVs

% define 2-dim residual - vec(1) = sigma, vec(2) = omega.
res_ev = @(vec) [integral(@(a) res_1_eig_int(a,vec(1),vec(2)),0,A)-1;
                 integral(@(a) res_2_eig_int(a,vec(1),vec(2)),0,A)];    

% --- k = 1:
x_01 = [-3.6;56]; % initial guess for k = 1 from surf plot

options = optimoptions('fsolve','Display','iter-detailed','Diagnostics','on','FunctionTolerance',.01);

[EV_01,~,~,output01] = fsolve(res_ev,x_01,options);

% --- k = 2:
x_02 = [-4.2;98]; % initial guess for k = 2 from surf plot

options = optimoptions('fsolve','Display','iter-detailed','Diagnostics','on','FunctionTolerance',.01);

[EV_02,~,~,output] = fsolve(res_ev,x_02,options);

%% ----- finding eigenvalues using fzero with deflation at zero EV - does not seem to give an advantage.

% % define 2-dim residual with deflation at origin - vec(1) = sigma, vec(2) = omega.
% % use 2-norm of (sigma,omega) to deflate.
% res_ev_deflation = @(vec) [(integral(@(a) res_1_eig_int(a,vec(1),vec(2)),0,A)-1)/norm(vec);
%                  integral(@(a) res_2_eig_int(a,vec(1),vec(2)),0,A)/norm(vec)];    
% 
% options = optimoptions('fsolve','Display','iter-detailed','Diagnostics','on','FunctionTolerance',.1);
% x_0 = [-3.6;56];
% 
% [test_ev,a,b,output] = fsolve(res_ev_deflation,x_0,options);

%% functions

function [y,yfinal,act_steps,res_history] = NewtonIteration(R_handle, RJ_handle, y_0, tol, steps)
% inputs:
% R_handle - function handle of residual fcn of output size Nx1
% RJ_handle - function handle of residual jacobian of output size NxN
% x_0 - initial value for iteration of size Nx1
% steps - maximum number of steps
% tol - tolerance of absolute error in norm of residual

% outputs:
% y - history of iteration
% yfinal - final value of estimated root
% act_steps - number of steps taken
% res_history - history of residual values

disp('starting newton iteration...')

y = zeros(length(y_0),steps);
res_history = zeros(length(R_handle([0,0])),steps);
y(:,1) = y_0;
res_history(:,1) = R_handle(y_0);

for k = 2:steps+1
    y(:,k) = y(:,k-1) - RJ_handle(y(:,k-1))\R_handle(y(:,k-1));
    res_history(:,k) = R_handle(y(:,k));
    if norm(R_handle(y(:,k)))<tol
        yfinal = y(:,k);
        y = y(:,1:k);
        res_history = res_history(:,1:k);
        act_steps = k-1;
        disp(['tolerance reached! # of steps taken: ',num2str(act_steps)])
        disp(['final values: (sigma,omega) = ',num2str(yfinal')])
        disp('stopping newton iteration...')
        return
    end
end
yfinal = y(:,end);
act_steps = steps;
disp(['maximum number of steps reached! # of steps taken: ',num2str(act_steps)])
disp(['final values: (sigma,omega) = ',num2str(yfinal')])
disp('stopping newton iteration...')
end
