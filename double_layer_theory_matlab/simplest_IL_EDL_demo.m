function simplest_IL_EDL_demo()
% Minimal demonstration of the paper's idea:
% c± from modified Boltzmann (with ls^2 d2phi, CS mu_ex and lv^2 d2mu_ex),
% coupled with differential-form Poisson: eps d2phi = -rho - ls^2 d2rho.

clc; clear; close all;

%% constants
eps0 = 8.8541878128e-12;
e0   = 1.602176634e-19;
kB   = 1.380649e-23;
NA   = 6.02214076e23;

%% parameters (paper-like)
T  = 300; beta = 1/(kB*T);
er = 10;  eps  = eps0*er;

d  = 0.5e-9;                 % ion diameter
R  = d/2;
v  = (4/3)*pi*R^3;            % ion volume
ls = d/sqrt(24);              % ell_s
lv = d/sqrt(40);              % ell_v

c0_M = 5;                     % 5 M
c0   = c0_M*1000;             % mol/m^3
c0n  = c0*NA;                 % number density (1/m^3), each species

% surface charge (choose small for simple Picard to behave)
qs_uCcm2 = 1;                 % try 1 first; you can increase later
qs = qs_uCcm2*1e-6/1e-4;       % C/m^2

%% 1D grid
dx = R/10;
L  = 60*R;
x  = (0:dx:L).';
N  = numel(x);

% exclusion zone: ions cannot approach closer than R
mask = (x >= R);

%% bulk mu_ex
p_bulk = v*(2*c0n);
mu_ex_bulk = mu_ex_hat_CS(p_bulk);

%% initial guess
phi = zeros(N,1);
cp  = c0n*ones(N,1); cm = c0n*ones(N,1);
cp(~mask)=0; cm(~mask)=0;

%% iteration controls (simple Picard + under-relax)
maxIter = 2000;
tol = 1e-6;
alpha_phi = 0.00002;
alpha_c   = 0.00002;

res = zeros(maxIter,1);

for it = 1:maxIter

    % second derivatives (simple central FD; we impose phi(L)=0 as reference)
    d2phi = d2_dirichlet_right(phi, dx, 0);

    % packing fraction p = v*(cp+cm)  (minimal version: no lv^2*d2c in pbar)
    p = v*(cp + cm);
    mu_ex = mu_ex_hat_CS(p);
    d2mu_ex = d2_dirichlet_right(mu_ex, dx, mu_ex_bulk);

    % modified Boltzmann (differential form)
    ueff = beta*e0*(phi + ls^2*d2phi);
    cp_new = c0n .* exp( - (+1)*ueff - mu_ex + mu_ex_bulk - lv^2*d2mu_ex );
    cm_new = c0n .* exp( - (-1)*ueff - mu_ex + mu_ex_bulk - lv^2*d2mu_ex );

    % apply exclusion zone and far-field bulk (simple)
    cp_new(~mask)=0; cm_new(~mask)=0;
    cp_new(end)=c0n; cm_new(end)=c0n;

    % relax concentrations
    cp = (1-alpha_c)*cp + alpha_c*cp_new;
    cm = (1-alpha_c)*cm + alpha_c*cm_new;

    % charge density and its Laplacian
    rho = e0*(cp - cm);
    d2rho = d2_dirichlet_right(rho, dx, 0);

    % differential-form Poisson: eps d2phi = -rho - ls^2 d2rho
    rhs = -(rho + ls^2*d2rho)/eps;         % d2phi = rhs
    phi_new = solve_poisson_mixed(rhs, dx, qs, eps);  % left Neumann, right Dirichlet phi(L)=0

    % relax potential
    phi_upd = (1-alpha_phi)*phi + alpha_phi*phi_new;

    res(it) = max(abs(phi_upd - phi));
    phi = phi_upd;

    if res(it) < tol
        res = res(1:it);
        fprintf('Converged at it=%d, res=%.3e\n', it, res(end));
        break;
    end
end

%% plots
figure; plot(x/d, beta*e0*phi, 'LineWidth',1.6);
xlabel('x/d'); ylabel('\beta e\phi'); grid on;
title(sprintf('Potential (q_s=%.1f \\muC/cm^2)', qs_uCcm2));

figure; plot(x/d, cp/c0n, 'LineWidth',1.6); hold on;
plot(x/d, cm/c0n, 'LineWidth',1.6);
xlabel('x/d'); ylabel('c_{\pm}/c_0'); grid on;
legend('c_+/c_0','c_-/c_0','Location','best');
title('Ion concentrations (modified Boltzmann idea)');

figure; semilogy(res, 'LineWidth',1.6);
xlabel('Iteration'); ylabel('max|\Delta\phi|'); grid on;
title('Picard iteration history');

end

%% --- helpers ---
function muhat = mu_ex_hat_CS(p)
% Carnahan–Starling: (8p - 9p^2 + 3p^3)/(1-p)^3
muhat = (8*p - 9*p.^2 + 3*p.^3) ./ (1 - p).^3;
end

function d2f = d2_dirichlet_right(f, dx, fR)
% simple second derivative, left Neumann(0) ghost, right Dirichlet f(N)=fR
N = numel(f);
d2f = zeros(N,1);

% left Neumann df/dx=0 => f0=f2
f0 = f(2);

% right Dirichlet: enforce f(N)=fR in stencil via ghost reflection
fNp1 = 2*fR - f(N-1);

d2f(1) = (f(2) - 2*f(1) + f0)/dx^2;
for i=2:N-1
    d2f(i) = (f(i+1) - 2*f(i) + f(i-1))/dx^2;
end
d2f(N) = (fNp1 - 2*f(N) + f(N-1))/dx^2;
end

function phi = solve_poisson_mixed(rhs, dx, qs, eps)
% d2phi = rhs
% left: -eps dphi/dx|0 = qs (Neumann)
% right: phi(L)=0 (Dirichlet)
N = numel(rhs);
A = spalloc(N,N,3*N);
b = rhs;

for i=2:N-1
    A(i,i-1)=1/dx^2; A(i,i)=-2/dx^2; A(i,i+1)=1/dx^2;
end

% left Neumann via ghost: phi0 = phi2 + 2dx*qs/eps
A(1,1)=-2/dx^2; A(1,2)=2/dx^2;
b(1)=rhs(1) - (2*qs)/(eps*dx);

% right Dirichlet
A(N,:)=0; A(N,N)=1; b(N)=0;

phi = A\b;
end
