% ============================================================
% Linear BSK-alpha (de Souza gradient expansion) + Goodwin short-range
% Single-run profile solver (NO capacitance scan)
%
% Adjustable parameters: a, b, gamma, u0
% (delta_c is kept as a fixed model parameter here; you can still change it.)
%
% alpha_eff = 1 / (1 + (gamma/2)*(a-b))
%
% Dimensionless coordinate: Z = z / lambda_D
% Dimensionless potential:  u(Z) = beta e phi(Z)
%
% Linear 4th-order ODE:
%   delta_c^4 u'''' + (2*delta_c^2 - 1/alpha_eff) u'' + u = 0
%
% BCs (potential-controlled):
%   u(0)     = u0
%   u'''(0)  = 0
%   u(Zmax)  = 0
%   u'(Zmax) = 0
%
% Linear concentration (charge-mode only):
%   c+/c0 = 1 - alpha_eff*u(Z)
%   c-/c0 = 1 + alpha_eff*u(Z)
%
% ============================================================

clear; clc; close all;

%% -------- Adjustable parameters (YOU change these) ----------
par.a     = 1.0;     % a (=a+=a-)  dimensionless
par.b     = 3;     % b           dimensionless
par.gamma = 0.90;    % bulk packing fraction gamma = n+_inf + n-_inf

par.u0    = 8.0;     % surface potential u(0) = u0 = beta e phi0  (YOUR "u")

%% -------- Fixed model parameters (can also change if needed) ----------
par.d = 4;        % nm
par.lambda_D = 5; % nm
par.delta_c = par.d / sqrt(24) / par.lambda_D;

par.Z0      = 0.0;   % left boundary (ion-center accessible boundary); set >0 if you want excluded layer
par.Zmax    = 30;    % far-field boundary

%% -------- Solver options ----------
opts = bvpset('RelTol', 1e-9, 'AbsTol', 1e-11, 'NMax', 40000);

%% -------- Derived alpha_eff(a,b,gamma) ----------
den = 1 + 0.5 * par.gamma * (par.a - par.b);  % 1 + (gamma/2)(a-b)
if den <= 0
    error('alpha_eff non-positive: 1 + (gamma/2)*(a-b)=%.6g <= 0. Adjust (a,b,gamma).', den);
end
par.alpha_eff = 1/den;

fprintf('alpha_eff = %.10f\n', par.alpha_eff);

%% -------- Solve BVP (single run) ----------
zmesh  = linspace(par.Z0, par.Zmax, 300);
solinit = bvpinit(zmesh, @(z) guess_u(z, par));

sol = bvp4c(@(z,y) odefun_bsk_alpha(z,y,par), ...
            @(ya,yb) bcfun_profile(ya,yb,par), ...
            solinit, opts);

%% -------- Evaluate solution ----------
Z = linspace(par.Z0, par.Zmax, 1200);
Y = deval(sol, Z);

u  = Y(1,:);   % u(Z)
uz = Y(2,:);   % du/dZ

%% -------- Linear concentrations (normalized) ----------
aeff = par.alpha_eff;

cplus_over_c0  = 1 - aeff * u;
cminus_over_c0 = 1 + aeff * u;

%% -------- Plot u(z) ----------
figure('Color','w');
plot(Z, u, 'LineWidth', 2.0);
xlabel('z/\lambda_D');
ylabel('u(z)=\beta e\phi(z)');
title(sprintf('u(z) profile: a=%.3g, b=%.3g, \\gamma=%.3g, u0=%.3g, \\delta_c=%.3g, \\alpha=%.3g', ...
    par.a, par.b, par.gamma, par.u0, par.delta_c, par.alpha_eff));
grid on;

%% -------- Plot c+/c0 and c-/c0 ----------
figure('Color','w');
plot(Z, cplus_over_c0, 'LineWidth', 2.0); hold on;
plot(Z, cminus_over_c0, 'LineWidth', 2.0);
yline(1,'--');
xlabel('z/\lambda_D');
ylabel('c_{\pm}(z)/c_0');
legend('c_+/c_0','c_-/c_0','bulk','Location','best');
title('Ion concentration profiles (linear response from u)');
grid on;

%% -------- Optional: warn if linear approximation becomes unphysical ----------
if any(cplus_over_c0 < 0) || any(cminus_over_c0 < 0)
    warning('Linear model gives negative concentration somewhere. Reduce |u0| or use nonlinear WDA+HS model.');
end

%% ============================================================
% ODE system: y = [u; uZ; uZZ; uZZZ]
% delta_c^4 u'''' + (2*dc^2 - 1/alpha) u'' + u = 0
% ============================================================
function dydz = odefun_bsk_alpha(~, y, par)
    dc   = par.delta_c;
    aeff = par.alpha_eff;

    if dc == 0
        error('delta_c=0 makes the 4th-order model singular. Use small delta_c or switch to 2nd-order DH.');
    end

    dydz = zeros(4,1);
    dydz(1) = y(2);
    dydz(2) = y(3);
    dydz(3) = y(4);

    coef = (2*dc^2 - 1/aeff);
    dydz(4) = ( -coef*y(3) - y(1) ) / (dc^4);
end

%% ============================================================
% Boundary conditions (profile, potential-controlled):
% u(Z0)=u0, u'''(Z0)=0, u(Zmax)=0, u'(Zmax)=0
% ============================================================
function res = bcfun_profile(ya, yb, par)
    res = zeros(4,1);
    res(1) = ya(1) - par.u0; % u(Z0)=u0
    res(2) = ya(4);          % u'''(Z0)=0
    res(3) = yb(1);          % u(Zmax)=0
    res(4) = yb(2);          % u'(Zmax)=0
end

%% ============================================================
% Initial guess
% ============================================================
function y0 = guess_u(z, par)
    % decaying exponential guess
    L = 2.0;
    A = par.u0;
    u    = A*exp(-(z-par.Z0)/L);
    uz   = -(A/L)*exp(-(z-par.Z0)/L);
    uzz  = (A/L^2)*exp(-(z-par.Z0)/L);
    uzzz = -(A/L^3)*exp(-(z-par.Z0)/L);
    y0 = [u; uz; uzz; uzzz];
end
