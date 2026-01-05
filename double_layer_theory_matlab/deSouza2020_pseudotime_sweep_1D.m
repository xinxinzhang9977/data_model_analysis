function deSouza2020_pseudotime_sweep_1D()
clc; clear; close all;

%% ---------------- Physical constants ----------------
eps0 = 8.8541878128e-12;
e0   = 1.602176634e-19;
kB   = 1.380649e-23;
NA   = 6.02214076e23;

%% ---------------- Parameters ----------------
T    = 300;
beta = 1/(kB*T);

er   = 10;
eps  = eps0*er;

d    = 0.5e-9;
R    = d/2;
vion = (4/3)*pi*R^3;

ls = d/sqrt(24);
lv = d/sqrt(40);

c0_M = 5;              % mol/L
c0   = c0_M*1000;      % mol/m^3
c0n  = c0*NA;          % number density per species

% ---- sweep list (uC/cm^2) ----
qs_uCcm2_list = [20];  % 你可以改成论文同款序列
qs_list = qs_uCcm2_list * 1e-6 / 1e-4;       % -> C/m^2

%% ---------------- Grid ----------------
dx = R/10;
L  = 60*R;
x  = (0:dx:L).';
N  = numel(x);
xd = x/d;

maskFluid = (x >= R);
idxFluid  = find(maskFluid);

%% ---------------- Bulk CS mu_ex ----------------
p0 = vion*(2*c0n);
mu_ex_bulk_hat = mu_ex_hat_CS(p0);

%% ---------------- Packing caps ----------------
ctot_max = 0.99 / vion;
ci_max   = 0.5 * ctot_max;

%% ---------------- Pseudo-time controls ----------------
M = 1.0;
dt = 0.05 * dx^2;

maxStep_eachQs = 6000;   % 每个qs最多推进步数
tol_mu  = 2e-7;
tol_phi = 2e-9;

% far-field clamped thickness
nbulk = 20;

%% ---------------- Storage for sweep results ----------------
phi_all   = zeros(N, numel(qs_list));
cp_all    = zeros(N, numel(qs_list));
cm_all    = zeros(N, numel(qs_list));
rho_all   = zeros(N, numel(qs_list));
f_all     = zeros(N, numel(qs_list));  % cumulative charge fraction
V0_all    = zeros(numel(qs_list),1);   % V0 = phi(0)-phi(L) (phi(L)=0)
conv_all  = zeros(numel(qs_list),3);   % [steps, final_resMu, final_resPhi]

%% ---------------- Initial fields (start from bulk) ----------------
phi = zeros(N,1);
cplus  = c0n*ones(N,1);
cminus = c0n*ones(N,1);
cplus(~maskFluid)=0; cminus(~maskFluid)=0;

%% ==================== SWEEP LOOP ====================
for k = 1:numel(qs_list)
    qs = qs_list(k);
    fprintf('\n=== Solving for qs = %.2f uC/cm^2 (%.3g C/m^2) ===\n', qs_uCcm2_list(k), qs);

    % pseudo-time march to steady state for this qs
    last_resMu = inf; last_resPhi = inf;
    for step = 1:maxStep_eachQs

        % ---- enforce bulk far-field strongly ----
        cplus(end-nbulk+1:end)  = c0n;
        cminus(end-nbulk+1:end) = c0n;
        phi(end) = 0;

        % ---- Poisson solve: eps d2phi = -rho - ls^2 d2rho ----
        rho = e0*(cplus - cminus);
        d2rho = d2dx2_mixed(rho, dx, 'neumann', 0, 'dirichlet', 0);

        rhs = -(rho + ls^2*d2rho)/eps;
        phi_new = solve_poisson_mixedBC_rhs(rhs, dx, qs, eps);

        % ---- d2phi for electrostatic correction ----
        d2phi = d2dx2_mixed(phi_new, dx, 'neumann', -qs/eps, 'dirichlet', 0);

        % ---- packing fraction p and CS mu_ex ----
        d2cplus  = d2dx2_mixed(cplus,  dx, 'neumann', 0, 'dirichlet', c0n);
        d2cminus = d2dx2_mixed(cminus, dx, 'neumann', 0, 'dirichlet', c0n);

        cbar_plus  = cplus  + lv^2 * d2cplus;
        cbar_minus = cminus + lv^2 * d2cminus;

        p = vion*(cbar_plus + cbar_minus);
        p = min(max(p, 1e-12), 1-1e-8);

        mu_ex_hat = mu_ex_hat_CS(p);
        d2mu_ex_hat = d2dx2_mixed(mu_ex_hat, dx, 'neumann', 0, 'dirichlet', mu_ex_bulk_hat);

        % ---- chemical potentials (dimensionless beta*mu) ----
        u_eff = beta*e0*(phi_new + ls^2*d2phi);

        mu_plus  = safe_log(cplus./c0n)  + (+1)*u_eff + (mu_ex_hat - mu_ex_bulk_hat) + lv^2*d2mu_ex_hat;
        mu_minus = safe_log(cminus./c0n) + (-1)*u_eff + (mu_ex_hat - mu_ex_bulk_hat) + lv^2*d2mu_ex_hat;

        % ---- gradient flow: ct = d/dx ( M*c * d(mu)/dx ) ----
        dmu_plus_dx  = ddx_central_mixed(mu_plus,  dx, 'neumann', 0, 'dirichlet', 0);
        dmu_minus_dx = ddx_central_mixed(mu_minus, dx, 'neumann', 0, 'dirichlet', 0);

        flux_plus  = -M * cplus  .* dmu_plus_dx;
        flux_minus = -M * cminus .* dmu_minus_dx;

        % no-flux at wall
        flux_plus(1)  = 0;
        flux_minus(1) = 0;

        dflux_plus_dx  = ddx_central_mixed(flux_plus,  dx, 'neumann', 0, 'dirichlet', 0);
        dflux_minus_dx = ddx_central_mixed(flux_minus, dx, 'neumann', 0, 'dirichlet', 0);

        cplus_new  = cplus  + dt * (-dflux_plus_dx);
        cminus_new = cminus + dt * (-dflux_minus_dx);

        % ---- enforce exclusion, positivity, packing ----
        cplus_new(~maskFluid)=0; cminus_new(~maskFluid)=0;
        cplus_new  = max(cplus_new,  0);
        cminus_new = max(cminus_new, 0);

        % cap each species
        cplus_new  = min(cplus_new,  ci_max);
        cminus_new = min(cminus_new, ci_max);

        % cap total packing
        ctot = cplus_new + cminus_new;
        over = ctot > ctot_max;
        if any(over)
            scale = ctot_max ./ ctot(over);
            cplus_new(over)  = cplus_new(over).*scale;
            cminus_new(over) = cminus_new(over).*scale;
        end

        % ---- residuals ----
        resPhi = max(abs(phi_new - phi));
        resMu  = max( abs(dmu_plus_dx(idxFluid)) + abs(dmu_minus_dx(idxFluid)) );

        % accept
        phi = phi_new;
        cplus = cplus_new;
        cminus = cminus_new;

        last_resMu = resMu;
        last_resPhi = resPhi;

        if (resMu < tol_mu) && (resPhi < tol_phi)
            fprintf('Converged: step=%d, resMu=%.3e, resPhi=%.3e\n', step, resMu, resPhi);
            break;
        end

        if mod(step,5000)==0
            fprintf(' step=%d, resMu=%.3e, resPhi=%.3e\n', step, resMu, resPhi);
        end
    end

    % ---- store results for this qs ----
    rho = e0*(cplus - cminus);
    phi_all(:,k) = phi;
    cp_all(:,k)  = cplus;
    cm_all(:,k)  = cminus;
    rho_all(:,k) = rho;

    % cumulative charge fraction: f(x) = -(1/qs) * ∫_0^x rho(x') dx'
    if abs(qs) > 0
        qcum = cumtrapz(x, rho);     % C/m^2
        f_all(:,k) = -qcum / qs;
    else
        f_all(:,k) = 0;
    end

    V0_all(k)   = phi(1) - phi(end);      % phi(end)=0 so V0=phi(1)
    conv_all(k,:) = [step, last_resMu, last_resPhi];
end

%% ==================== PLOTS ====================

% 1) Potential profiles (beta e phi)
figure;
hold on;
for k = 1:numel(qs_list)
    plot(xd, beta*e0*phi_all(:,k), 'LineWidth', 1.4);
end
xlabel('x/d'); ylabel('\beta e \phi'); grid on;
title('Potential profiles for different q_s');
legend(compose('%.0f', qs_uCcm2_list), 'Location','best');
set(gca,'Box','on');

% 2) Concentration profiles (show a subset to avoid clutter)
figure;
kshow = unique(round(linspace(1,numel(qs_list),min(numel(qs_list),6))));
tiledlayout(2,1);

nexttile;
hold on;
for kk = kshow
    plot(xd, cp_all(:,kk)/c0n, 'LineWidth', 1.2);
end
xlabel('x/d'); ylabel('c_+/c_0'); grid on;
title('Cation concentration profiles (subset)');
legend(compose('qs=%.0f', qs_uCcm2_list(kshow)), 'Location','best');

nexttile;
hold on;
for kk = kshow
    plot(xd, cm_all(:,kk)/c0n, 'LineWidth', 1.2);
end
xlabel('x/d'); ylabel('c_-/c_0'); grid on;
title('Anion concentration profiles (subset)');
legend(compose('qs=%.0f', qs_uCcm2_list(kshow)), 'Location','best');

% 3) Cumulative charge fraction f(x)
figure;
hold on;
for k = 2:numel(qs_list) % skip qs=0
    plot(xd, f_all(:,k), 'LineWidth', 1.4);
end
xlabel('x/d'); ylabel('f(x) = -\int_0^x \rho dx / q_s'); grid on;
title('Cumulative charge fraction');
legend(compose('%.0f', qs_uCcm2_list(2:end)), 'Location','best');

% 4) Differential capacitance from sweep: C = dqs/dV0
% Use central finite difference in (qs, V0)
[Vsort, idx] = sort(V0_all);
qsort = qs_list(idx);

C = zeros(size(qsort));
for i = 2:numel(qsort)-1
    dqdV = (qsort(i+1)-qsort(i-1)) / (Vsort(i+1)-Vsort(i-1));
    C(i) = dqdV; % F/m^2
end
% one-sided ends
C(1) = (qsort(2)-qsort(1)) / (Vsort(2)-Vsort(1));
C(end) = (qsort(end)-qsort(end-1)) / (Vsort(end)-Vsort(end-1));

figure;
plot(Vsort, C, 'LineWidth', 1.8);
xlabel('V_0 = \phi(0) - \phi(L) (V)'); ylabel('C = d q_s / dV_0 (F/m^2)');
grid on; title('Differential capacitance (from sweep)');

% 5) Print convergence summary
fprintf('\n=== Convergence summary ===\n');
fprintf('qs(uC/cm2)   steps    resMu      resPhi\n');
for k = 1:numel(qs_list)
    fprintf('%8.1f   %6.0f  %9.2e  %9.2e\n', qs_uCcm2_list(k), conv_all(k,1), conv_all(k,2), conv_all(k,3));
end

end

%% ================= Helper functions =================

function muhat = mu_ex_hat_CS(p)
muhat = (8*p - 9*p.^2 + 3*p.^3) ./ (1 - p).^3;
end

function y = safe_log(x)
y = log(max(x, 1e-300));
end

function d2f = d2dx2_mixed(f, dx, leftType, leftVal, rightType, rightVal)
N = numel(f);
d2f = zeros(N,1);

if strcmpi(leftType,'neumann')
    f0 = f(2) - 2*dx*leftVal;
else
    error('leftType not supported');
end

if strcmpi(rightType,'dirichlet')
    fNp1 = 2*rightVal - f(N-1);
else
    error('rightType not supported');
end

d2f(1) = (f(2) - 2*f(1) + f0)/dx^2;
for i=2:N-1
    d2f(i) = (f(i+1) - 2*f(i) + f(i-1))/dx^2;
end
d2f(N) = (fNp1 - 2*f(N) + f(N-1))/dx^2;
end

function df = ddx_central_mixed(f, dx, leftType, leftVal, rightType, rightVal)
N = numel(f);
df = zeros(N,1);

if strcmpi(leftType,'neumann')
    f0 = f(2) - 2*dx*leftVal;
else
    error('leftType not supported');
end

if strcmpi(rightType,'dirichlet')
    fNp1 = 2*rightVal - f(N-1);
else
    error('rightType not supported');
end

df(1) = (f(2) - f0)/(2*dx);
for i=2:N-1
    df(i) = (f(i+1) - f(i-1))/(2*dx);
end
df(N) = (fNp1 - f(N-1))/(2*dx);
end

function phi = solve_poisson_mixedBC_rhs(rhs, dx, qs, eps)
N = numel(rhs);
A = spalloc(N,N,3*N);
b = rhs;

for i=2:N-1
    A(i,i-1) = 1/dx^2;
    A(i,i)   = -2/dx^2;
    A(i,i+1) = 1/dx^2;
end

% left Neumann: -eps dphi/dx|0 = qs
A(1,1) = -2/dx^2;
A(1,2) =  2/dx^2;
b(1)   = rhs(1) - (2*qs)/(eps*dx);

% right Dirichlet: phi(L)=0
A(N,:) = 0;
A(N,N) = 1;
b(N)   = 0;

phi = A\b;
end
