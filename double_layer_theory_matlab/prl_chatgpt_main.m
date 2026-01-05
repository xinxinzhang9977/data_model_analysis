clear; clc;

%% ===== Physical parameters (PRL-like) =====
kB   = 1.380649e-23;
T    = 300;
beta = 1/(kB*T);
e    = 1.602176634e-19;

eps0 = 8.8541878128e-12;
epsr = 10;
eps  = eps0*epsr;

R  = 0.5e-9;                 % ion radius
v  = 4/3*pi*R^3;             % ion volume
c0 = 3/v;                    % bulk conc (~IL)
qs = 20e-6;                  % surface charge (C/m^2)

%% ===== Spatial grid =====
L  = 12e-9;
Nx = 1200;
x  = linspace(0,L,Nx)';
dx = x(2)-x(1);

%% ===== Build nonlocal weights =====
[Ws, Wv] = build_weights(x, R);

%% ===== Initial fields =====
phi        = zeros(Nx,1);
mu_ex_bar  = zeros(Nx,1);

%% ===== Iteration parameters =====
maxIter = 4000;
tol     = 1e-6;
alpha   = 0.2;

%% ===== Poisson matrix =====
A = poisson_matrix(Nx, dx, eps, qs);

%% ===== Main iteration =====
for it = 1:maxIter

    phi_old = phi;

    % --- nonlocal electrostatic potential ---
    phi_bar = Ws * phi;

    % --- Boltzmann distributions (Eq.6) ---
    c_plus  = c0 * exp(-beta*e*phi_bar - beta*mu_ex_bar);
    c_minus = c0 * exp(+beta*e*phi_bar - beta*mu_ex_bar);

    % --- nonlocal packing fraction ---
    p_bar = v * (Wv * (c_plus + c_minus));

    % --- Carnahan–Starling chemical potential ---
    mu_ex = (8*p_bar - 9*p_bar.^2 + 3*p_bar.^3) ./ (1 - p_bar).^3;

    % --- nonlocal excess chemical potential ---
    mu_ex_bar = Wv * mu_ex;

    % --- weighted charge density ---
    rho_bar = e * (c_plus - c_minus);
    rho_eff = Ws * rho_bar;

    % --- solve Poisson equation ---
    phi = A \ (-rho_eff);

    % --- relaxation ---
    phi = (1-alpha)*phi_old + alpha*phi;

    % --- convergence ---
    if max(abs(phi-phi_old)) < tol
        fprintf('Converged in %d iterations\n',it);
        break
    end
end

%% ===== Derived quantities =====
phi_bar = Ws * phi;
c_plus  = c0 * exp(-beta*e*phi_bar - beta*mu_ex_bar);
c_minus = c0 * exp(+beta*e*phi_bar - beta*mu_ex_bar);
rho     = e*(c_plus - c_minus);

%% ===== Plot: potential =====
figure;
plot(x*1e9, beta*e*phi,'LineWidth',2)
xlabel('x (nm)')
ylabel('\beta e \phi')
title('Electrostatic potential (PRL Fig.2-like)')
grid on

%% ===== Plot: charge density =====
figure;
plot(x*1e9, rho/max(abs(rho)),'LineWidth',2)
xlabel('x (nm)')
ylabel('\rho_e (normalized)')
title('Charge density oscillations (overscreening)')
grid on
