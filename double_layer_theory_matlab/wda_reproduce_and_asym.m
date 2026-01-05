function wda_reproduce_and_asym()
clc; clear;

%% ============================================================
%  Choose a case
%  "symmetric"  : reproduce (Rp=Rm, zp=-zm, vp=vm)
%  "asymmetric" : Rp!=Rm, zp and zm arbitrary, vp!=vm (fully asymmetric)
%% ============================================================
case_name = "symmetric";  % <-- change: "symmetric" or "asymmetric"

%% ------------------------
% Physical constants
%% ------------------------
kB   = 1.380649e-23;
T    = 300;
beta = 1/(kB*T);
e0   = 1.602176634e-19;
NA   = 6.02214076e23;

eps0 = 8.8541878128e-12;
epsr = 10;
eps  = eps0*epsr;

%% ============================================================
%  Model parameters (edit here)
%% ============================================================

switch case_name
    case "symmetric"
        % ---- symmetric ions ----
        Rp = 0.25e-9;     Rm = 0.25e-9;     % radii
        zp = +1;          zm = -1;          % valences
        p0 = 0.30;                            % bulk packing fraction (0.1~0.5 typical)
        u0_target = -5;                       % target dimensionless electrode potential u0=beta e phi0
        L_factor = 60;                         % domain length L = L_factor * R
        dz_factor = 10;                        % dz = R/dz_factor

    case "asymmetric"
        % ---- fully asymmetric example ----
        Rp = 0.40e-9;     Rm = 0.80e-9;      % radii (different -> different volumes & kernels)
        zp = +2;          zm = -1;           % valences (asymmetric)
        p0 = 0.35;                             % IMPORTANT: keep <=~0.5 for stability; can ramp later
        u0_target = -5;                        % target u0
        L_factor = 40;                          % for large Rm, choose e.g. 30~80
        dz_factor = 12;                         % dz = Rmin/dz_factor (resolve smaller ion)
end

% Ion volumes
vp = 4*pi*Rp^3/3;
vm = 4*pi*Rm^3/3;

%% ------------------------
% Bulk concentrations from (packing fraction p0) + (electroneutrality)
% electroneutrality: zp*c0p + zm*c0m = 0  (zm should be negative)
% p0 = vp*c0p + vm*c0m
%% ------------------------
[c0p, c0m] = bulk_from_p0_and_neutrality(p0, vp, vm, zp, zm);

fprintf('Bulk check:\n');
fprintf('  Rp=%.3f nm, Rm=%.3f nm\n', Rp*1e9, Rm*1e9);
fprintf('  zp=%.3g, zm=%.3g\n', zp, zm);
fprintf('  c0p=%.3f M, c0m=%.3f M\n', c0p/(1000*NA), c0m/(1000*NA));
fprintf('  neutrality zp*c0p+zm*c0m = %.3e (should be ~0)\n', zp*c0p+zm*c0m);
fprintf('  p_bulk = %.3f (target p0=%.3f)\n\n', vp*c0p+vm*c0m, p0);

%% ------------------------
% Grid
%% ------------------------
Rmin = min(Rp, Rm);
Rmax = max(Rp, Rm);

L  = L_factor * Rmax;
dz = Rmin / dz_factor;
N  = max(50, round(L/dz));
z  = ((1:N) - 0.5) * dz;   % cell centers

% Exclusion zones (ion centers cannot enter z < Ri)
mask_p = (z >= Rp);
mask_m = (z >= Rm);

%% ------------------------
% Build kernels for each species
% ws: surface kernel, wv: volume kernel (1D forms used in SI)
%% ------------------------
nRp = round(Rp/dz);
nRm = round(Rm/dz);

[ws_p, wv_p] = build_1d_weights(Rp, vp, dz, nRp);
[ws_m, wv_m] = build_1d_weights(Rm, vm, dz, nRm);

%% ------------------------
% Carnahan–Starling bulk excess chemical potential (dimensionless beta*mu)
%% ------------------------
p_bulk = vp*c0p + vm*c0m;
beta_mu_bulk = carnahan_starling_beta_mu(p_bulk);

%% ============================================================
%  Continuation in u0 (highly recommended for stability)
%% ============================================================
u0_list = ramp_list(0, u0_target, 60);   % 15 steps; increase if needed

% initial guess for first step
phi        = zeros(N,1);                 % V
beta_mu_ex = beta_mu_bulk*ones(N,1);     % dimensionless beta*mu_ex

for ku = 1:numel(u0_list)
    u0 = u0_list(ku);
    phi0 = u0/(beta*e0);

    fprintf('=== Continuation step %d/%d: u0 = %.4g ===\n', ku, numel(u0_list), u0);

    [phi, beta_mu_ex, out] = solve_wda_step( ...
        phi, beta_mu_ex, ...
        z, dz, eps, phi0, ...
        beta, e0, ...                      % <<< add
        Rp, Rm, nRp, nRm, ...
        ws_p, ws_m, wv_p, wv_m, ...
        mask_p, mask_m, ...
        zp, zm, c0p, c0m, vp, vm, ...
        beta_mu_bulk);

    fprintf('  converged=%d, it=%d, qs=%.3e C/m^2\n\n', out.converged, out.it, out.qs);
end

%% ============================================================
%  Final recompute for plots
%% ============================================================
phi_bar_p = conv_band_norm_trap(phi, ws_p, dz, nRp);
phi_bar_m = conv_band_norm_trap(phi, ws_m, dz, nRm);
mu_bar_p  = conv_band_norm_trap(beta_mu_ex, wv_p, dz, nRp);
mu_bar_m  = conv_band_norm_trap(beta_mu_ex, wv_m, dz, nRm);

arg_p = -zp*beta*e0*phi_bar_p - mu_bar_p + beta_mu_bulk;
arg_m = -zm*beta*e0*phi_bar_m - mu_bar_m + beta_mu_bulk;

[arg_p, arg_m] = clip_args(arg_p, arg_m);

c_p = c0p * exp(arg_p);
c_m = c0m * exp(arg_m);

% enforce exclusion
c_p(~mask_p) = 0;
c_m(~mask_m) = 0;

% optional physical ceiling (avoid local packing>~1)
c_p = min(c_p, 0.95/vp);
c_m = min(c_m, 0.95/vm);

u = beta*e0*phi;

qs_out = eps*( (u0_target/(beta*e0)) - phi(1) ) / dz;

%% ------------------------
% Plots
%% ------------------------
figure('Name','Potential');
subplot(2,1,1);
plot(z*1e9, phi, 'LineWidth', 1.5); grid on
xlabel('z (nm)'); ylabel('\phi(z) (V)');
title(sprintf('Potential, final u_0=%.3g, q_s=%.3e C/m^2', u0_target, qs_out));

subplot(2,1,2);
plot(z*1e9, u, 'LineWidth', 1.5); grid on
xlabel('z (nm)'); ylabel('u(z)=\\beta e \\phi(z)');

figure('Name','Concentrations');
subplot(2,1,1);
plot(z*1e9, c_p, 'LineWidth', 1.5); hold on
plot(z*1e9, c_m, 'LineWidth', 1.5); grid on
xlabel('z (nm)'); ylabel('c_i(z) (m^{-3})');
legend('c_+(z)','c_-(z)','Location','best');
title('Number density');

subplot(2,1,2);
semilogy(z*1e9, max(c_p./c0p, 1e-20), 'LineWidth', 1.5); hold on
semilogy(z*1e9, max(c_m./c0m, 1e-20), 'LineWidth', 1.5); grid on
xlabel('z (nm)'); ylabel('c_i(z)/c_{i,bulk} (log)');
legend('c_+/c_{+,bulk}','c_-/c_{-,bulk}','Location','best');
title('Normalized (log scale)');

end

%% ============================================================
%  One continuation step solver (Picard + under-relaxation)
%% ============================================================
function [phi, beta_mu_ex, out] = solve_wda_step( ...
    phi, beta_mu_ex, ...
    z, dz, eps, phi0, ...
    beta, e0, ...                      % <<< add
    Rp, Rm, nRp, nRm, ...
    ws_p, ws_m, wv_p, wv_m, ...
    mask_p, mask_m, ...
    zp, zm, c0p, c0m, vp, vm, ...
    beta_mu_bulk)


N = numel(z);

% iteration params (tune here)
maxIter = 30000;
tol     = 1e-6;

% under-relaxation (these usually work; adjust if needed)
w_phi = 0.02;
w_mu  = 0.05;

% packing cap (too close to 1 => CS blows up)
p_cap = min(0.90, 0.98*(vp*c0p+vm*c0m));  % you can also set fixed 0.6~0.9

% store error history for debug/plot if you want
err_phi_hist = zeros(maxIter,1);
err_mu_hist  = zeros(maxIter,1);

for it = 1:maxIter
    % weighted potentials
    phi_bar_p = conv_band_norm_trap(phi, ws_p, dz, nRp);
    phi_bar_m = conv_band_norm_trap(phi, ws_m, dz, nRm);

    % weighted beta*mu_ex (dimensionless)  <-- IMPORTANT: do NOT divide by beta
    mu_bar_p  = conv_band_norm_trap(beta_mu_ex, wv_p, dz, nRp);
    mu_bar_m  = conv_band_norm_trap(beta_mu_ex, wv_m, dz, nRm);

    % keep finite inside exclusion (they don't matter because c=0 there)
    if any(~mask_p)
        i0 = find(mask_p,1,'first');
        phi_bar_p(~mask_p) = phi_bar_p(i0);
        mu_bar_p(~mask_p)  = mu_bar_p(i0);
    end
    if any(~mask_m)
        i0 = find(mask_m,1,'first');
        phi_bar_m(~mask_m) = phi_bar_m(i0);
        mu_bar_m(~mask_m)  = mu_bar_m(i0);
    end

    % concentrations (Eq-like update)
    arg_p = -zp*beta*e0*phi_bar_p - mu_bar_p + beta_mu_bulk;
    arg_m = -zm*beta*e0*phi_bar_m - mu_bar_m + beta_mu_bulk;

    [arg_p, arg_m] = clip_args(arg_p, arg_m);

    c_p = c0p * exp(arg_p);
    c_m = c0m * exp(arg_m);

    % exclusion zones
    c_p(~mask_p) = 0;
    c_m(~mask_m) = 0;

    % optional local ceilings (prevents runaway)
    c_p = min(c_p, 0.95/vp);
    c_m = min(c_m, 0.95/vm);

    % weighted concentrations -> packing fraction
    cbar_p = conv_band_norm_trap(c_p, wv_p, dz, nRp);
    cbar_m = conv_band_norm_trap(c_m, wv_m, dz, nRm);

    p_bar = vp*cbar_p + vm*cbar_m;
    p_bar = min(max(p_bar, 0), p_cap);

    % CS beta*mu_ex(p_bar)
    beta_mu_new = carnahan_starling_beta_mu(p_bar);
    beta_mu_new = min(max(beta_mu_new, 0), 60);

    % weighted charge density (species-dependent smearing)
    % rho_bar = e*( zp * (c_p * ws_p) + zm * (c_m * ws_m) )
    rho_bar = e0*( zp*conv_band_norm_trap(c_p, ws_p, dz, nRp) + ...
                   zm*conv_band_norm_trap(c_m, ws_m, dz, nRm) );

    % solve Poisson with left ghost Dirichlet phi0, right gauge phi(N)=0
    phi_new = solve_poisson_dirichlet_neumann_dense(rho_bar, eps, dz, phi0);


    % under-relax
    phi_next      = (1-w_phi)*phi + w_phi*phi_new;
    beta_mu_next  = (1-w_mu )*beta_mu_ex + w_mu*beta_mu_new;

    err_phi = max(abs(phi_next - phi));
    err_mu  = max(abs(beta_mu_next - beta_mu_ex));

    err_phi_hist(it) = err_phi;
    err_mu_hist(it)  = err_mu;

    phi = phi_next;
    beta_mu_ex = beta_mu_next;

    if mod(it,200)==0 || it==1
        qs_out = eps*(phi0 - phi(1))/dz;
        fprintf('  it=%d  err_phi=%.3e  err_mu=%.3e  qs=%.3e\n', it, err_phi, err_mu, qs_out);
    end

    if max(err_phi, err_mu) < tol
        out.converged = true;
        out.it = it;
        out.qs = eps*(phi0 - phi(1))/dz;
        out.err_phi_hist = err_phi_hist(1:it);
        out.err_mu_hist  = err_mu_hist(1:it);
        return;
    end
end

out.converged = false;
out.it = maxIter;
out.qs = eps*(phi0 - phi(1))/dz;
out.err_phi_hist = err_phi_hist;
out.err_mu_hist  = err_mu_hist;

end

%% ============================================================
%  Helper: build weights
%% ============================================================
function [ws, wv] = build_1d_weights(R, vion, dz, nR)
offset = (-nR:nR)';
zz = offset*dz;
inside = (abs(zz) <= R);

ws = zeros(size(zz));
wv = zeros(size(zz));

% surface kernel: ws = 1/(2R) within |z|<=R
ws(inside) = 1/(2*R);

% volume kernel: wv = pi*(R^2 - z^2)/v within |z|<=R
wv(inside) = pi*(R^2 - zz(inside).^2)/vion;
end

%% ============================================================
%  Helper: normalized band convolution using trapezoidal endpoints
%  ybar_i = (sum_j w_ij y_j dz) / (sum_j w_ij dz)
%% ============================================================
function ybar = conv_band_norm_trap(y, w, dz, nR)
N = numel(y);
ybar = zeros(N,1);

for i = 1:N
    j1 = max(1, i-nR);
    j2 = min(N, i+nR);

    off = (j1:j2) - i;
    ww  = w(off + nR + 1);

    % trapezoid weights
    trap = ones(size(ww));
    if numel(trap) >= 2
        trap(1)   = 0.5;
        trap(end) = 0.5;
    end

    denom = sum(trap .* ww) * dz;
    if denom == 0
        ybar(i) = 0;
    else
        ybar(i) = sum(trap .* ww .* y(j1:j2)) * dz / denom;
    end
end
end

%% ============================================================
%  Helper: Poisson solver (dense, readable, no spdiags)
%  eps * (phi_{i+1}-2phi_i+phi_{i-1})/dz^2 = rho_i
%  Left boundary via ghost phi0, right gauge phi(N)=phiL
%% ============================================================
function phi = solve_poisson_dirichlet_neumann_dense(rho_bar, eps, dz, phi0)
% Solve:  eps * phi'' = -rho_bar
% Left:   phi(0)=phi0 (ghost point)
% Right:  phi'(L)=0  -> phi_{N+1} = phi_N

N = numel(rho_bar);
rhs = -(dz^2/eps) * rho_bar;   % IMPORTANT sign

A = zeros(N,N);
b = rhs;

% i=1: phi_2 - 2phi_1 + phi0 = rhs1  -> (-2)phi1 + (1)phi2 = rhs1 - phi0
A(1,1) = -2;
if N>=2, A(1,2)=1; end
b(1) = b(1) - phi0;

% interior: i=2..N-1
for i=2:N-1
    A(i,i-1)=1; A(i,i)=-2; A(i,i+1)=1;
end

% i=N with Neumann: phi_{N+1}=phi_N
% equation: phi_{N+1} - 2phi_N + phi_{N-1} = rhs_N
% -> phi_N - 2phi_N + phi_{N-1} = rhs_N  -> (1)*phi_{N-1} + (-1)*phi_N = rhs_N
A(N,:) = 0;
A(N,N-1) = 1;
A(N,N)   = -1;
b(N)     = rhs(N);

phi = A\b;
end

%% ============================================================
%  Helper: CS beta*mu_ex
%% ============================================================
function beta_mu = carnahan_starling_beta_mu(p)
beta_mu = (8*p - 9*p.^2 + 3*p.^3) ./ (1 - p).^3;
end

%% ============================================================
%  Helper: stable exp argument clipping
%% ============================================================
function [ap, am] = clip_args(ap, am)
% Keep exponentials in a safe numeric range.
% If you care about co-ion visibility, do NOT clamp too low (e.g. -80).
ap = max(min(ap,  80), -60);
am = max(min(am,  80), -60);
end

%% ============================================================
%  Helper: compute bulk concentrations from p0 + electroneutrality
%  zp*c0p + zm*c0m = 0
%  p0 = vp*c0p + vm*c0m
%% ============================================================
function [c0p, c0m] = bulk_from_p0_and_neutrality(p0, vp, vm, zp, zm)
if zm >= 0
    error('Please pass a negative zm for an anion (e.g., zm=-1, -2).');
end
ratio = -zp/zm;   % c0m = ratio * c0p
denom = vp + vm*ratio;
if denom <= 0
    error('Invalid parameters: vp + vm*(-zp/zm) must be > 0.');
end
c0p = p0 / denom;
c0m = ratio * c0p;
end

%% ============================================================
%  Helper: build a ramp list from a to b
%% ============================================================
function arr = ramp_list(a, b, n)
if n <= 1
    arr = b;
    return;
end
arr = linspace(a, b, n);
end
