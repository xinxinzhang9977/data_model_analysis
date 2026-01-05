clear; clc;

%% =========================================================
%  Physical parameters
%% =========================================================
kB   = 1.380649e-23;
T    = 300;
beta = 1/(kB*T);
e0   = 1.602176634e-19;

eps0 = 8.8541878128e-12;
epsr = 10;
eps  = eps0*epsr;

R  = 0.5e-9;
v  = 4/3*pi*R^3;
c0 = 5/v;
qs = 0.6;

%% =========================================================
%  Spatial grid
%% =========================================================
L  = 2e-9;
Nx = 200;
x  = linspace(0,L,Nx)';
dx = x(2)-x(1);

%% =========================================================
%  Nonlocal weights
%% =========================================================
[Ws, Wv] = build_weights(x, R);

%% =========================================================
%  Poisson operator
%% =========================================================
A = poisson_matrix(Nx, dx, eps);

%% =========================================================
%  Initial log-densities
%% =========================================================
u_p = zeros(Nx,1);
u_m = zeros(Nx,1);

%% =========================================================
%  Time stepping
%% =========================================================
dt      = 1e-16;
maxIter = 50000;
tol     = 1e-8;

F_hist  = nan(maxIter,1);

%% =========================================================
%  Main gradient flow loop
%% =========================================================
for it = 1:maxIter

    % ---- densities ----
    c_p = c0 * exp(u_p);
    c_m = c0 * exp(u_m);

    % ---- packing fraction ----
    p_bar = v * (Wv * (c_p + c_m));
    p_bar = min(p_bar, 0.9);

    % ---- excess chemical potential ----
    mu_ex = (8*p_bar - 9*p_bar.^2 + 3*p_bar.^3) ./ (1 - p_bar).^3;
    mu_ex_bar = Wv * mu_ex;

    % ---- electrostatics ----
    rho = e0 * (c_p - c_m);
    rho_bar = Ws * rho;

    b = -rho_bar;
    b(1)   = -qs;
    b(end) = 0;

    phi = A \ b;
    phi_bar = Ws * phi;

    % ---- log-density gradient flow (semi-implicit) ----
    u_p_new = (u_p - dt*( beta*e0*phi_bar + beta*mu_ex_bar )) ...
                / (1 + dt);
    u_m_new = (u_m - dt*( -beta*e0*phi_bar + beta*mu_ex_bar )) ...
                / (1 + dt);

    % ---- free energy ----
    F = sum( c_p.*(log(c_p/c0)-1) + ...
             c_m.*(log(c_m/c0)-1) ) * dx ...
        + sum( mu_ex .* (c_p+c_m) ) * dx ...
        + 0.5 * sum( phi .* rho ) * dx;

    F_hist(it) = F;

    % ---- convergence ----
    err = max([norm(u_p_new-u_p,inf), norm(u_m_new-u_m,inf)]);

    u_p = u_p_new;
    u_m = u_m_new;

    if err < tol
        fprintf('Converged at iteration %d\n', it);
        F_hist = F_hist(1:it);
        break
    end
end

%% =========================================================
%  Post-processing
%% =========================================================
c_p = c0 * exp(u_p);
c_m = c0 * exp(u_m);
rho = e0 * (c_p - c_m);

%% ---- Potential ----
figure;
plot(x*1e9, beta*e0*phi, 'LineWidth',2)
xlabel('x (nm)')
ylabel('\beta e \phi')
title('Electrostatic potential')
grid on

%% ---- Charge density ----
figure;
plot(x*1e9, rho/max(abs(rho)), 'LineWidth',2)
xlabel('x (nm)')
ylabel('\rho_e (normalized)')
title('Charge density')
grid on

%% ---- Free energy ----
figure;
plot(F_hist,'LineWidth',2)
xlabel('Iteration')
ylabel('Free energy')
title('Free energy decay (log-density gradient flow)')
grid on
