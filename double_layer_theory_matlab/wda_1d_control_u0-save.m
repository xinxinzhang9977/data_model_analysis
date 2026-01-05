function wda_1d_control_u0_asymR()
clc; clear;

%% ------------------------
% Physical parameters
%% ------------------------
kB   = 1.380649e-23;
T    = 300;
beta = 1/(kB*T);
e0   = 1.602176634e-19;
NA   = 6.02214076e23;

eps0 = 8.8541878128e-12;
epsr = 2;
eps  = eps0*epsr;

% --------- Asymmetric ion radii ---------
Rp = 0.40e-9;   % R_+  (m)
Rm = 0.80e-9;   % R_-  (m)

vp = 4*pi*Rp^3/3;
vm = 4*pi*Rm^3/3;

zval = [+1, -1];

% bulk concentration (per species)  (you can set different if needed)
c0p_M = 10;  % mol/L
c0m_M = 10;  % mol/L
c0p   = c0p_M*1000*NA;  % 1/m^3
c0m   = c0m_M*1000*NA;  % 1/m^3

%% ------------------------
% Grid
% Choose dz based on the smaller radius for resolution
%% ------------------------
Rmin = min(Rp, Rm);
Rmax = max(Rp, Rm);

L  = 10*Rmax;        % you used 10*R before; now use max radius
dz = Rmin/10;        % resolve the smaller ion
N  = round(L/dz);
z  = ((1:N) - 0.5)*dz;

% species-dependent exclusion (no ion centers closer than their radius)
mask_p = (z >= Rp);
mask_m = (z >= Rm);

%% ------------------------
% Control u0 = beta*e*phi0
%% ------------------------
u0   = -5;                 % dimensionless electrode potential
phi0 = u0/(beta*e0);      % volts

%% ------------------------
% Build species-dependent weights (1D SI kernels)
%% ------------------------
nRp = round(Rp/dz);
nRm = round(Rm/dz);

[ws_p, wv_p] = build_1d_weights(Rp, vp, dz, nRp);
[ws_m, wv_m] = build_1d_weights(Rm, vm, dz, nRm);

%% ------------------------
% Bulk mu_ex from total packing fraction (simple CS closure)
% p_bulk = v+ c+ + v- c-
%% ------------------------
p_bulk = vp*c0p + vm*c0m;
beta_mu_bulk = carnahan_starling_beta_mu(p_bulk);

%% ------------------------
% Initial guess
%% ------------------------
phi        = linspace(phi0, 0, N)';   % V
beta_mu_ex = beta_mu_bulk*ones(N,1);  % still one field (simple closure)

%% ------------------------
% Iteration controls
%% ------------------------
maxIter = 8000;
tol     = 1e-4;
w_phi   = 0.00001;
w_mu    = 0.005;

err_phi_hist = nan(maxIter,1);
err_mu_hist  = nan(maxIter,1);

%% ------------------------
% Main iteration
%% ------------------------
for it = 1:maxIter

    % ---- weighted potentials for each species ----
    phi_bar_p = conv_band_norm(phi, ws_p, dz, nRp);
    phi_bar_m = conv_band_norm(phi, ws_m, dz, nRm);

    % ---- weighted mu_ex (use volume kernel; could also choose ws, but keep consistent with your current) ----
    mu_bar_p  = conv_band_norm(beta_mu_ex, wv_p, dz, nRp);
    mu_bar_m  = conv_band_norm(beta_mu_ex, wv_m, dz, nRm);

    % keep finite values inside exclusion zones (not physically used since c=0 there)
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

    % ---- concentration update (Eq.6), species uses its own phi_bar and mu_bar ----
    arg_p = -zval(1)*beta*e0*phi_bar_p - mu_bar_p + beta_mu_bulk;
    arg_m = -zval(2)*beta*e0*phi_bar_m - mu_bar_m + beta_mu_bulk;

    arg_p = max(min(arg_p,  80), -80);
    arg_m = max(min(arg_m,  80), -80);

    c_p = c0p * exp(arg_p);
    c_m = c0m * exp(arg_m);

    % species-dependent exclusion
    c_p(~mask_p) = 0;
    c_m(~mask_m) = 0;

    % ---- weighted concentrations for packing fraction ----
    cbar_p = conv_band_norm(c_p, wv_p, dz, nRp);
    cbar_m = conv_band_norm(c_m, wv_m, dz, nRm);

    p_bar = vp*cbar_p + vm*cbar_m;

    % cap for stability
    p_cap = 0.40;
    p_bar = min(max(p_bar, 0), p_cap);

    beta_mu_new = carnahan_starling_beta_mu(p_bar);
    beta_mu_new = min(max(beta_mu_new, 0), 60);

    % ---- weighted charge density (species-dependent smearing) ----
    rho_bar = e0*( zval(1)*conv_band_norm(c_p, ws_p, dz, nRp) + ...
                   zval(2)*conv_band_norm(c_m, ws_m, dz, nRm) );

    % ---- Poisson solve ----
    phi_new = solve_poisson_dirichlet2(rho_bar, eps, dz, phi0, 0.0);

    % ---- relax ----
    phi_next      = (1-w_phi)*phi + w_phi*phi_new;
    mu_next       = (1-w_mu)*beta_mu_ex + w_mu*beta_mu_new;

    err_phi = max(abs(phi_next - phi));
    err_mu  = max(abs(mu_next  - beta_mu_ex));

    if ~isfinite(err_phi), err_phi = inf; end
    if ~isfinite(err_mu),  err_mu  = inf; end

    err_phi_hist(it) = err_phi;
    err_mu_hist(it)  = err_mu;

    phi = phi_next;
    beta_mu_ex = mu_next;

    if mod(it,50)==0 || it==1
        qs_out = eps*(phi0 - phi(1))/dz;
        fprintf('it=%d, err_phi=%.3e, err_mu=%.3e, qs=%.3e C/m^2\n', it, err_phi, err_mu, qs_out);
    end

    if max(err_phi, err_mu) < tol
        qs_out = eps*(phi0 - phi(1))/dz;
        fprintf('Converged at it=%d, qs=%.6e C/m^2\n', it, qs_out);
        break;
    end
end

err_phi_hist = err_phi_hist(1:it);
err_mu_hist  = err_mu_hist(1:it);

%% ------------------------
% Recompute final profiles for plotting
%% ------------------------
phi_bar_p = conv_band_norm(phi, ws_p, dz, nRp);
phi_bar_m = conv_band_norm(phi, ws_m, dz, nRm);
mu_bar_p  = conv_band_norm(beta_mu_ex, wv_p, dz, nRp);
mu_bar_m  = conv_band_norm(beta_mu_ex, wv_m, dz, nRm);

arg_p = -zval(1)*beta*e0*phi_bar_p - mu_bar_p + beta_mu_bulk;
arg_m = -zval(2)*beta*e0*phi_bar_m - mu_bar_m + beta_mu_bulk;
arg_p = max(min(arg_p,80),-80);
arg_m = max(min(arg_m,80),-80);

c_p = c0p * exp(arg_p);  c_p(~mask_p)=0;
c_m = c0m * exp(arg_m);  c_m(~mask_m)=0;

u = beta*e0*phi;
qs_out = eps*(phi0 - phi(1))/dz;

%% ------------------------
% Plots
%% ------------------------
figure('Name','Iteration errors');
semilogy(1:numel(err_phi_hist), err_phi_hist, 'LineWidth', 1.5); hold on
semilogy(1:numel(err_mu_hist),  err_mu_hist,  'LineWidth', 1.5);
grid on
xlabel('Iteration'); ylabel('max norm error');
legend('err\_phi','err\_\mu','Location','best');
title('Convergence history');

figure('Name','Potential vs z');
subplot(2,1,1);
plot(z*1e9, phi, 'LineWidth', 1.5); grid on
xlabel('z (nm)'); ylabel('\phi(z) (V)');
title(sprintf('Potential (u_0=%.3g), q_s=%.3e C/m^2', u0, qs_out));
subplot(2,1,2);
plot(z*1e9, u, 'LineWidth', 1.5); grid on
xlabel('z (nm)'); ylabel('u(z)=\beta e \phi(z)');

figure('Name','Ion concentrations vs z');
subplot(2,1,1);
plot(z*1e9, c_p, 'LineWidth', 1.5); hold on
plot(z*1e9, c_m, 'LineWidth', 1.5); grid on
xlabel('z (nm)'); ylabel('c_i(z) (m^{-3})');
legend('c_+(z)','c_-(z)','Location','best');
title(sprintf('Ion number densities (R_+=%.2f nm, R_-=%.2f nm)', Rp*1e9, Rm*1e9));

subplot(2,1,2);
plot(z*1e9, c_p./c0p, 'LineWidth', 1.5); hold on
plot(z*1e9, c_m./c0m, 'LineWidth', 1.5); grid on
xlabel('z (nm)'); ylabel('c_i(z)/c_{i,bulk}');
legend('c_+/c_{+,bulk}','c_-/c_{-,bulk}','Location','best');
title('Normalized densities');

end

%% ========================= Helper functions =========================

function [ws, wv] = build_1d_weights(R, vion, dz, nR)
offset = (-nR:nR)';
zz = offset*dz;
Theta = @(x) double(x >= 0);
ws = (1/(2*R)) * Theta(R - abs(zz));
wv = (pi*(R^2 - zz.^2)/vion) .* Theta(R - abs(zz));
end

function ybar = conv_band_norm(y, w, dz, nR)
N = numel(y);
ybar = zeros(N,1);
for i = 1:N
    j1 = max(1, i-nR);
    j2 = min(N, i+nR);
    off = (j1:j2) - i;
    ww  = w(off + nR + 1);
    denom = sum(ww)*dz;
    if denom == 0
        ybar(i) = 0;
    else
        ybar(i) = sum(ww .* y(j1:j2))*dz / denom;
    end
end
end

function phi = solve_poisson_dirichlet2(rho_bar, eps, dz, phi0, phiL)
N = numel(rho_bar);
rhs = (dz^2/eps) * rho_bar;

main = 2*ones(N,1);
sub  = -ones(N-1,1);
sup  = -ones(N-1,1);

b = rhs;
b(1) = b(1) + phi0;

% enforce phi(N)=phiL (gauge)
main(N) = 1;
sub(N-1) = 0;
b(N) = phiL;

A = spdiags([[sub;0], main, [0;sup]], [-1,0,1], N, N);
phi = A\b;
end

function beta_mu = carnahan_starling_beta_mu(p)
beta_mu = (8*p - 9*p.^2 + 3*p.^3) ./ (1 - p).^3;
end
