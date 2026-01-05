
clear; clc;

%% =========================================================
%  Physical parameters (PRL 2020 – de Souza et al.)
%% =========================================================
kB   = 1.380649e-23;
T    = 300;
beta = 1/(kB*T);
e    = 1.602176634e-19;

eps0 = 8.8541878128e-12;
epsr = 2;
eps  = eps0*epsr;

R  = 0.25e-9;                % ion radius (m)
v  = 4/3*R^3;            % ion volume
c0 = 0.5/v;               % bulk concentration (SAFE < close packing)
qs = 0.6;                 % surface charge density (C/m^2)

%% =========================================================
%  Spatial grid
%% =========================================================
L  = 2e-9;
Nx = 200;
x  = linspace(0,L,Nx)';
dx = x(2)-x(1);

%% =========================================================
%  Build nonlocal weights (surface & volume)
%% =========================================================
[Ws, Wv] = build_weights(x, R);

%% =========================================================
%  Initial conditions (BREAK symmetry!)
%% =========================================================
phi        =  0.01 * exp(-x/R);   % small bias near electrode
mu_ex_bar  = zeros(Nx,1);

%% =========================================================
%  Iteration parameters
%% =========================================================
maxIter = 200000;
tol     = 1e-5;
alpha   = 0.00001;
err_hist = nan(maxIter,1);   % record max |phi - phi_old|

%% =========================================================
%  Poisson matrix (Neumann at x=0, Dirichlet at bulk)
%% =========================================================
A = poisson_matrix(Nx, dx, eps);

%% =========================================================
%  Main self-consistent iteration
%% =========================================================
for it = 1:maxIter
    %disp(it);
    if it < 100
        mu_ex_bar = zeros(Nx,1);
    end

    phi_old = phi;

    % ---- nonlocal electrostatic potential ----
    phi_bar = Ws * phi;

    % ---- Boltzmann distributions (Eq. 6) ----
    arg_plus  = -beta*e*phi_bar - beta*mu_ex_bar;
    arg_minus = +beta*e*phi_bar - beta*mu_ex_bar;
    
    arg_plus  = max(min(arg_plus,  50), -50);
    arg_minus = max(min(arg_minus, 50), -50);
    
    c_plus  = c0 * exp(-beta*e*phi_bar - beta*mu_ex_bar);
    c_minus = c0 * exp(+beta*e*phi_bar - beta*mu_ex_bar);

mask = (x >= R);
c_plus(~mask)  = 0;
c_minus(~mask) = 0;

    % ---- nonlocal packing fraction ----
    p_bar = v * (Wv * (c_plus + c_minus));
    p_bar = min(p_bar, 0.95);   % avoid CS divergence

    % ---- Carnahan–Starling excess chemical potential ----
    mu_ex = (8*p_bar - 9*p_bar.^2 + 3*p_bar.^3) ./ (1 - p_bar).^3;

    % ---- nonlocal excess chemical potential ----
    mu_ex_bar = Wv * mu_ex;

    % ---- weighted charge density (Poisson source) ----
    rho_bar = e * (c_plus - c_minus);
    rho_eff = Ws * rho_bar;

    % ---- Poisson equation: A*phi = b ----
    b = -rho_eff;
    b(1)   = -qs;   % surface charge enters RHS !!
    b(end) = 0;     % reference potential in bulk

    phi = A \ b;

    % ---- relaxation ----
    phi = (1-alpha)*phi_old + alpha*phi;
    
    % ---- record convergence metric ----
    err = max(abs(phi - phi_old));
    err_hist(it) = err;

    % ---- convergence check ----
    if err < tol
        fprintf('Converged in %d iterations\n', it);
        err_hist = err_hist(1:it);   % truncate unused entries
        break
    end

end

%% =========================================================
%  Post-processing (Fig.2 / Fig.3 style)
%% =========================================================
phi_bar = Ws * phi;
c_plus  = c0 * exp(-beta*e*phi_bar - beta*mu_ex_bar);
c_minus = c0 * exp(+beta*e*phi_bar - beta*mu_ex_bar);
rho     = e * (c_plus - c_minus);
mask = (x >= R);
c_plus(~mask)  = 0;
c_minus(~mask) = 0;
%% ---- Electrostatic potential (Fig.2-like) ----
figure;
plot(x*1e9, beta*e*phi, 'LineWidth', 2)
xlabel('x (nm)')
ylabel('\beta e \phi')
title('Electrostatic potential (nonlocal PB)')
grid on

%% ---- Charge density oscillations (Fig.2 / Fig.3-like) ----
figure;
plot(x*1e9, c_plus, 'LineWidth', 2)
hold on;
plot(x*1e9, c_minus, 'LineWidth', 2)
xlabel('x (nm)')
ylabel('\rho_e (normalized)')
title('Charge density oscillations and overscreening')
grid on
%% =========================================================
%  Convergence history
%% =========================================================
figure;
semilogy(err_hist, 'LineWidth', 2)
xlabel('Iteration')
ylabel('max |\phi^{(n)} - \phi^{(n-1)}|')
title('Self-consistent iteration convergence')
grid on
%% =========================================================
%  END OF MAIN SCRIPT
%% =========================================================


%% =========================================================
%  Local functions
%% =========================================================

function [Ws, Wv] = build_weights(x, R)
    Nx = length(x);
    dx = x(2)-x(1);

    Ws = zeros(Nx);
    Wv = zeros(Nx);

    for i = 1:Nx
        for j = 1:Nx
            d = abs(x(i)-x(j));

            % --- surface (charge) weight ---
            if d <= R
                Ws(i,j) = 1/(2*R);
            end

            % --- volume (packing) weight ---
            if d <= R
                Wv(i,j) = 1*(R^2-d^2)/(2*R);%Wv(i,j) = 3*(R^2-d^2)/(4*R^3);
            end
        end
    end

    Ws = Ws * dx;
    Wv = Wv * dx;
end


function A = poisson_matrix(Nx, dx, eps)
    A = sparse(Nx, Nx);

    % bulk points
    for i = 2:Nx-1
        A(i,i-1) =  eps/dx^2;
        A(i,i)   = -2*eps/dx^2;
        A(i,i+1) =  eps/dx^2;
    end

    % Neumann BC at x = 0  (surface charge via RHS)
    A(1,1) = -eps/dx;
    A(1,2) =  eps/dx;

    % Dirichlet BC at bulk
    A(Nx,Nx) = 1;
end
