function deSouza2020_differential_form_1D()
% Differential-form (gradient expansion) solver for ionic-liquid EDL
% Based on Supplemental Material equations (S11) etc.
% - Picard iteration with under-relaxation
% - Finite differences in 1D, Neumann BCs (surface charge at x=0, bulk at x=L)
%
% References:
%  - Main text & Supplemental Material of:
%    J. P. de Souza et al., PRL 125, 116001 (2020).

clc; clear; close all;

%% ------------------ Physical constants ------------------
eps0 = 8.8541878128e-12;   % F/m
e0   = 1.602176634e-19;    % C
kB   = 1.380649e-23;       % J/K
NA   = 6.02214076e23;      % 1/mol

%% ------------------ Model parameters (paper-like) ------------------
T    = 300;                % K
beta = 1/(kB*T);
er   = 10;                 % relative permittivity
eps  = eps0*er;

d    = 0.5e-9;             % ion diameter (m)
R    = d/2;
vion = (4/3)*pi*R^3;       % volume of one ion (m^3)

% Gradient-expansion length scales (from Fourier expansion)
ls = d/sqrt(24);           % ell_s
lv = d/sqrt(40);           % ell_v

% Bulk concentration: c0 = 5 M for each species in symmetric 1:1 IL model
c0_M = 5;                  % mol/L
c0   = c0_M*1000;          % mol/m^3
c0n  = c0*NA;              % number density (1/m^3) for each species

% Surface charge density (example)
qs_uCcm2 = 10;             % microC/cm^2
qs = qs_uCcm2 * 1e-6 / 1e-4; % -> C/m^2

%% ------------------ Numerical grid (recommend dx=R/10, L=60R) ------------------
dx = R/10;                 % suggested in SM numerical section for integral form; works well here too
L  = 60*R;
x  = (0:dx:L).';           % include boundaries
N  = numel(x);

% Hard-sphere exclusion zone near wall:
% In differential-form section they place BCs at x=d/2; we set ions = 0 for x < R
maskFluid = (x >= R);      % fluid region with mobile ion centers

%% ------------------ Bulk excess chemical potential (Carnahan-Starling) ------------------
p_bulk = vion*(2*c0n);                         % filling fraction in bulk (symmetric, c+=c-=c0)
mu_ex_bulk_hat = mu_ex_hat_CS(p_bulk);         % dimensionless beta*mu_ex

%% ------------------ Iteration controls ------------------
maxIter = 5000;
tolPhi  = 1e-9;            % V
tolCrel = 1e-8;            % relative concentration tolerance
alphaPhi = 0.05;           % under-relaxation for potential
alphaC   = 0.03;           % under-relaxation for concentrations

%% ------------------ Initial guess ------------------
phi = zeros(N,1);          % V
cplus  = c0n*ones(N,1);
cminus = c0n*ones(N,1);
cplus(~maskFluid)  = 0;
cminus(~maskFluid) = 0;

resHist = nan(maxIter,3);  % [resPhi, resCrel, resEqPoisson]

%% ------------------ Main Picard loop ------------------
for it = 1:maxIter

    % --- second derivative of phi (for the ls^2 * d2phi term in Boltzmann factor)
    % phi(L)=0, left dphi/dx given by surface charge, 但在 d2phi 里只需要一个"左端导数条件"
    d2phi = d2dx2_mixed(phi, dx, 'neumann', -qs/eps, 'dirichlet', 0);

    % --- compute weighted filling fraction pbar ≈ v Σ (ci + lv^2 d2ci)
    d2cplus  = d2dx2_mixed(cplus,  dx, 'neumann', 0, 'dirichlet', c0n);
    d2cminus = d2dx2_mixed(cminus, dx, 'neumann', 0, 'dirichlet', c0n);
    cbar_plus  = cplus  + lv^2 * d2cplus;
    cbar_minus = cminus + lv^2 * d2cminus;

    pbar = vion*(cbar_plus + cbar_minus);
    % keep pbar in (0,1) to avoid blow-up
    pbar = min(max(pbar, 1e-12), 1-1e-8);

    % --- excess chemical potential (dimensionless) and its Laplacian
    mu_ex_hat = mu_ex_hat_CS(pbar);
    d2mu_ex_hat = d2dx2_mixed(mu_ex_hat, dx, 'neumann', 0, 'dirichlet', mu_ex_bulk_hat);

    % --- update concentrations from differential-form Boltzmann (S11)
    % ci = c0 exp( -zi*beta*e*(phi + ls^2 d2phi) - mu_ex_hat + mu_ex_bulk_hat - lv^2 d2mu_ex_hat )
    u_eff = beta*e0*(phi + ls^2*d2phi);  % dimensionless electrostatic term

    arg_plus  = ( -(+1)*u_eff - mu_ex_hat + mu_ex_bulk_hat - lv^2*d2mu_ex_hat );
    arg_minus = ( -(-1)*u_eff - mu_ex_hat + mu_ex_bulk_hat - lv^2*d2mu_ex_hat );

    % 指数截断，避免溢出（±50 对 double 已经很安全）
    arg_plus  = min(max(arg_plus,  -20), 20);
    arg_minus = min(max(arg_minus, -20), 20);

    cplus_new  = c0n .* exp(arg_plus);
    cminus_new = c0n .* exp(arg_minus);
    cplus(end)=c0n; cminus(end)=c0n;

    % enforce exclusion zone
    cplus_new(~maskFluid)  = 0;
    cminus_new(~maskFluid) = 0;

    % under-relax concentrations
    cplus  = (1-alphaC)*cplus  + alphaC*cplus_new;
    cminus = (1-alphaC)*cminus + alphaC*cminus_new;

    % --- charge density rho_e and its Laplacian
    rho = e0*(cplus - cminus);           % C/m^3
    d2rho = d2dx2_mixed(rho, dx, 'neumann', 0, 'dirichlet', 0);

    % --- solve Poisson with smearing term in differential form:
    % eps * d2phi = -rho - ls^2 * d2rho
    rhs = -(rho + ls^2*d2rho)/eps;       % RHS for d2phi = rhs

    phi_new = solve_poisson_mixedBC_rhs(rhs, dx, qs, eps);
    % 不再需要 phi_new = phi_new - phi_new(end);


    % under-relax potential
    phi_updated = (1-alphaPhi)*phi + alphaPhi*phi_new;

    % --- residuals
    resPhi = max(abs(phi_updated - phi));
    cscale = max(c0n, 1); % avoid /0
    resCrel = max( abs([cplus_new-cplus; cminus_new-cminus]) ) / cscale;

    % check Poisson equation residual with updated phi
d2phi_chk = d2dx2_mixed(phi_updated, dx, 'neumann', -qs/eps, 'dirichlet', 0);    resEq = max(abs(eps*d2phi_chk + rho + ls^2*d2rho)); % C/m^3 scale

    resHist(it,:) = [resPhi, resCrel, resEq];

    phi = phi_updated;

    if (resPhi < tolPhi) && (resCrel < tolCrel)
        resHist = resHist(1:it,:);
        fprintf('Converged at iter %d: resPhi=%.3e V, resCrel=%.3e\n', it, resPhi, resCrel);
        break;
    end

    if mod(it,200)==0
        fprintf('Iter %d: resPhi=%.3e V, resCrel=%.3e, resEq=%.3e\n', it, resPhi, resCrel, resEq);
    end

    if it==maxIter
        fprintf('Reached maxIter without full convergence.\n');
    end
end

%% ------------------ Plot results ------------------
xd = x/d;

figure;
plot(xd, beta*e0*phi, 'LineWidth', 1.8);
xlabel('x/d'); ylabel('\beta e \phi');
title(sprintf('Potential profile (qs = %.1f \\muC/cm^2)', qs_uCcm2));
grid on;

figure;
plot(xd, cplus/c0n, 'LineWidth', 1.8); hold on;
plot(xd, cminus/c0n, 'LineWidth', 1.8);
xlabel('x/d'); ylabel('c_{\pm}/c_0');
title('Ion concentration profiles');
legend('c_+/c_0','c_-/c_0','Location','best');
grid on;

figure;
semilogy(resHist(:,1), 'LineWidth', 1.8); hold on;
semilogy(resHist(:,2), 'LineWidth', 1.8);
semilogy(resHist(:,3), 'LineWidth', 1.8);
xlabel('Iteration'); ylabel('Residual (log scale)');
title('Iteration history');
legend('||\Delta\phi||_\infty','max|\Delta c|/c_0','Poisson eq residual','Location','best');
grid on;

end

%% ====================== Helper functions ======================

function muhat = mu_ex_hat_CS(p)
% Carnahan-Starling excess chemical potential (dimensionless beta*mu_ex)
% muhat = (8p - 9p^2 + 3p^3) / (1-p)^3
muhat = (8*p - 9*p.^2 + 3*p.^3) ./ (1 - p).^3;
end

function d2f = d2dx2_mixed(f, dx, leftType, leftVal, rightType, rightVal)
% 1D second derivative with ghost-point BCs
% leftType:  'neumann' => df/dx = leftVal
% rightType: 'dirichlet' => f = rightVal   (at last node)

N = numel(f);
d2f = zeros(N,1);

% ----- left ghost (i=1 uses f0) -----
if strcmpi(leftType,'neumann')
    % (f2 - f0)/(2dx) = leftVal  => f0 = f2 - 2dx*leftVal
    f0 = f(2) - 2*dx*leftVal;
else
    error('leftType not supported');
end

% ----- right ghost (i=N uses fNp1) -----
if strcmpi(rightType,'dirichlet')
    % enforce f(N)=rightVal strongly; for second derivative at N-1 need fNp1:
    % simplest: mirror with fixed boundary: fNp1 = 2*rightVal - f(N-1)
    % also overwrite f(N) outside before calling if you want strict.
    fNp1 = 2*rightVal - f(N-1);
else
    error('rightType not supported');
end

% boundaries
d2f(1) = (f(2) - 2*f(1) + f0)/dx^2;
d2f(N) = (fNp1 - 2*f(N) + f(N-1))/dx^2;

% interior
for i=2:N-1
    d2f(i) = (f(i+1) - 2*f(i) + f(i-1))/dx^2;
end
end

function phi = solve_poisson_mixedBC_rhs(rhs, dx, qs, eps)
% Solve d2phi = rhs with:
%   -eps dphi/dx|_{0} = qs   (Neumann at x=0)
%    phi(L) = 0              (Dirichlet reference at far boundary)

N = numel(rhs);

A = spalloc(N,N,3*N);
b = rhs;

% Interior nodes
for i=2:N-1
    A(i,i-1) = 1/dx^2;
    A(i,i)   = -2/dx^2;
    A(i,i+1) = 1/dx^2;
end

% Left boundary Neumann using ghost:
% dphi/dx|0 = -qs/eps, phi(0)=phi(2)+2dx*qs/eps
A(1,1) = -2/dx^2;
A(1,2) =  2/dx^2;
b(1)   = rhs(1) - (2*qs)/(eps*dx);

% Right boundary Dirichlet: phi(N)=0
A(N,:) = 0;
A(N,N) = 1;
b(N)   = 0;

phi = A\b;
end
