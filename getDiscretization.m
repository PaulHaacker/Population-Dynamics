function [A_mat, C_mat, phi] = getDiscretization(parameter)
%% finding the matrices for the Galerkin-based simulation

% choose the basis phi(a) = [zero eigenfunction; two first eigenfunction
% pairs and function for the IC]
% --> requires two pairs of nonzero eigenvalues, computed by [Schmidt17]

%% extract input parameters:

A = parameter.A; % max age - double
mu = parameter.mu; % constant mortality rate - double
mu_int = parameter.mu_int; % mortality rate integral - function
k = parameter.k; % birth kernel - function handle
p = parameter.p; % output kernel - function handle
D_star = parameter.D_star; % steady-state dilution rate - double

% parameters for IC - paper [Schmidt17]
x0 = parameter.x0; % function handle

% eigenvalues of the form lambda = -sigma/A+-j*omega/(2*pi*A)

sigma(1) = parameter.sigma(1);
omega(1) = parameter.omega(1);
sigma(2) = parameter.sigma(2);
omega(2) = parameter.omega(2);

N_EV = length(sigma); % number of nonzero eigenvalues considered

sign_ImaginaryPart = 1; % only works for +1
EV = -sigma/A + 1i*omega/(2*pi*A)*sign_ImaginaryPart;

%% ------ verify compatibility of parameters

% verify steady-state input D_star:
phi_0 = @(a) exp(-D_star*a-mu_int(a)); % eigenfunction of steady-state operator with eigenvalue = 0
res_LS_int = @(a,D_star) k(a).*phi_0(a); % integrand of the Lotka-Sharpe condition

res_LotkaSharpe_fcn = @(D_star) integral(@(a) res_LS_int(a,D_star),0,A)-1; % Lotka-Sharpe in residual form as fcn of D_star

res_LS_test = res_LotkaSharpe_fcn(D_star);

if abs(res_LS_test)>.1
    error('steady-state input D_star is incompatible with Lotka-Sharpe!')
end

% verify eigenvalues:

res_1_eig_int = @(a,sigma,omega) k(a).*cos(omega.*a/2/pi/A).*exp(-sigma.*a/A).*phi_0(a);
res_2_eig_int = @(a,sigma,omega) k(a).*sin(omega.*a/2/pi/A).*exp(-sigma.*a/A).*phi_0(a);
res_ev_fcn = @(vec) [integral(@(a) res_1_eig_int(a,vec(1),vec(2)),0,A)-1;
                 integral(@(a) res_2_eig_int(a,vec(1),vec(2)),0,A)];  % define 2-dim residual - vec(1) = sigma, vec(2) = omega.

for kk = 1:N_EV
    res_ev(kk) = norm(res_ev_fcn([sigma(kk),omega(kk)]));
end

if max(res_ev)>.1
    error('eigenvalues are incompatible with Lotka-Sharpe!')
end

%% ------ basis of trial functions

phi = cell(6,1);
phi{1} = @(a) phi_0(a);
for kk = 1:N_EV
    phi{2*kk} = @(a) sign_ImaginaryPart*sin(omega(kk).*a/(2*pi*A)).*exp(-sigma(kk).*a/A).*phi{1}(a);
    phi{2*kk+1} = @(a) cos(omega(kk).*a/(2*pi*A)).*exp(-sigma(kk).*a/A).*phi{1}(a);
end
phi{2*N_EV+2} = @(a) x0(a); % initial condition fcn

% apply differential operator to basis:

D_phi = cell(6,1);
D_phi{1} = @(a) 0;

Lambda_mat = @(sigma,omega) [   -sign_ImaginaryPart*sigma/A,    sign_ImaginaryPart*omega/2/pi/A;
                                -omega/2/pi/A,                  -sigma/A];

for kk = 1:N_EV
    Lambda_k = Lambda_mat(sigma(kk),omega(kk));
    
    D_phi{2*kk} = @(a) Lambda_k(1,1)*phi{2*kk}(a) + Lambda_k(1,2)*phi{2*kk+1}(a);
    D_phi{2*kk+1} = @(a) Lambda_k(2,1)*phi{2*kk}(a) + Lambda_k(2,2)*phi{2*kk+1}(a);
end


% differential operator applied to IC
syms a_sym
x0_sym = x0(a_sym);
D_x0_sym = diff(x0_sym) + (mu(a_sym)+D_star)*x0_sym;

D_phi{2*N_EV+2} = matlabFunction(D_x0_sym);

%% Plot Basis functions
% % 
% % Q: Does the sign of an eigenfunction matter?
% % A: let f(a) satisfy the ODE 
% %           Df(a) = \lambda f(a)
% % where D is a differential operator. Let 
% %           g(a) := - f(a).
% % Whenever D is a linear operator, g also satisfies the above ODE and thus
% % even the function with the sign flipped is an eigenfunction..
% %  --> This is verified by switching the signs in the sine-argument, or
% %  equivalently, flipping the signs of the whole eigenfucntions.
% 
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
A_mat = eye(N)*D_star - Phi_1\Phi_2;

% output matrix, where y(t) = C*lambda(t)
C_mat = zeros(size(phi))';
for ii = 1: length(phi)
   C_mat(ii) = integral(@(a)p(a).*phi{ii}(a),0,A);
end

end