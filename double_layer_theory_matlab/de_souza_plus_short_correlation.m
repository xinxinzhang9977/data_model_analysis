% ============================================================
% de Souza (gradient expansion) + Goodwin short-range (a,b)
% LINEAR BSK-alpha model in 1D planar geometry
%
% What you get:
% 1) sigma_tilde(u0) by scanning surface potential u0
% 2) C_tilde(u0) = d sigma_tilde / d u0
% 3) At max voltage (u0_max): plots of
%    - potential profile u(z)
%    - ion concentration profiles c_+(z)/c0 and c_-(z)/c0 (linear response)
%
% Controls: a (=a+=a-), b, gamma, delta_c
% alpha_eff = 1 / (1 + (gamma/2)*(a-b))
%
% Dimensionless:
%   X = x / lambda_D
%   u = beta e phi
%   sigma_tilde = beta e lambda_D sigma / epsilon
%
% ODE (in X):
%   delta_c^4 u'''' + (2*delta_c^2 - 1/alpha_eff) u'' + u = 0
%
% BCs (potential scan):
%   u(XR)    = u0
%   u'''(XR) = 0
%   u(Xmax)  = 0
%   u'(Xmax) = 0
%
% Notes:
% - This is a LINEAR model. For too large |u0|, linear concentration may go negative.
% - XR is the ion-center accessible boundary. If you want excluded layer, set XR>0.
% ============================================================

clear; clc; close all;

%% ---- User controls ----
par.a       = 2.0;     % short-range self (dimensionless)
par.b       = -1.5;     % short-range cross (dimensionless)
par.gamma   = 0.30;    % bulk packing fraction gamma = n+_inf + n-_inf (dimensionless)
par.delta_c = 0.40;    % delta_c = l_s / lambda_D

par.XR      = 0.0;     % XR = R/lambda_D (set >0 if you enforce excluded layer)
par.Xmax    = 30;      % far-field boundary (increase if decay is slow / oscillatory)

% Surface potential scan: u0 = beta e phi0
u0_list = linspace(-8, 8, 81);

% bvp4c options
opts = bvpset('RelTol',1e-8,'AbsTol',1e-10,'NMax',30000);

%% ---- Derived alpha_eff(a,b,gamma) ----
den = 1 + 0.5*par.gamma*(par.a - par.b);
if den <= 0
    error('alpha_eff non-positive: 1 + (gamma/2)*(a-b)=%.6g <= 0. Adjust (a,b,gamma).', den);
end
par.alpha_eff = 1/den;
fprintf('alpha_eff = %.8f\n', par.alpha_eff);

%% ---- Solve for each u0 (continuation) ----
sigma_list = nan(size(u0_list));

% initial mesh
xmesh = linspace(par.XR, par.Xmax, 250);

sol_prev = [];
u0_prev  = NaN;

% store solutions
[~, idx0] = min(abs(u0_list - 0));
sol_at0    = [];
sol_at_u0min = [];
sol_at_u0max = [];

for k = 1:numel(u0_list)
    par.u0 = u0_list(k);

    if isempty(sol_prev)
        % first run guess
        solinit = bvpinit(xmesh, @(x) guess_u(x, par));
    else
        du0 = par.u0 - u0_prev;
        % robust continuation: guess from previous solution + shift in u only
        solinit = bvpinit(xmesh, @(x) guess_from_prev(x, sol_prev, du0));
    end

    sol = bvp4c(@(x,y) odefun_bsk_alpha(x,y,par), ...
                @(ya,yb) bcfun_potential_scan(ya,yb,par), ...
                solinit, opts);

    % sigma_tilde = -u'(XR)
    ya = deval(sol, par.XR);
    sigma_list(k) = -ya(2);

    % store special solutions
    if k == 1
        sol_at_u0min = sol;
    end
    if k == numel(u0_list)
        sol_at_u0max = sol;
    end
    if k == idx0
        sol_at0 = sol;
    end

    % update continuation state
    sol_prev = sol;
    u0_prev  = par.u0;

    if mod(k,10)==0 || k==1 || k==numel(u0_list)
        fprintf('k=%d/%d, u0=%+.3f, sigma_tilde=%+.6f\n', ...
            k, numel(u0_list), par.u0, sigma_list(k));
    end
end

%% ---- Differential capacitance (dimensionless) ----
Ctilde_list = gradient(sigma_list, u0_list);

%% ---- Plot sigma(u0) ----
figure('Color','w');
plot(u0_list, sigma_list, 'LineWidth', 1.8);
xlabel('u_0 = \beta e \phi_0');
ylabel('\sigma~ = \beta e \lambda_D \sigma / \epsilon');
title(sprintf('\\sigma~(u_0): a=%.3g, b=%.3g, \\gamma=%.3g, \\delta_c=%.3g, \\alpha=%.3g', ...
    par.a, par.b, par.gamma, par.delta_c, par.alpha_eff));
grid on;

%% ---- Plot C(u0) ----
figure('Color','w');
plot(u0_list, Ctilde_list, 'LineWidth', 1.8);
xlabel('u_0 = \beta e \phi_0');
ylabel('C~ = d\sigma~/du_0');
title('Dimensionless differential capacitance');
grid on;

%% ---- Profile at u0 ~ 0 (optional) ----
if ~isempty(sol_at0)
    X = linspace(par.XR, par.Xmax, 800);
    Y = deval(sol_at0, X);
    figure('Color','w');
    plot(X, Y(1,:), 'LineWidth', 1.8);
    xlabel('z/\lambda_D'); ylabel('u(z)');
    title(sprintf('Potential profile at u0=%.3g', u0_list(idx0)));
    grid on;
end

%% ---- Profiles at max voltage: potential and concentrations ----
if ~isempty(sol_at_u0max)
    X = linspace(par.XR, par.Xmax, 1000);
    Y = deval(sol_at_u0max, X);
    u = Y(1,:);

    % Linear concentration profiles (charge-mode only):
    % n0 = gamma/2, and c±/c0 = n±/n0 = 1 ∓ alpha_eff * u
    aeff = par.alpha_eff;

    cplus_over_c0  = 1 - aeff*u;
    cminus_over_c0 = 1 + aeff*u;

    figure('Color','w');
    plot(X, u, 'LineWidth', 1.8);
    xlabel('z/\lambda_D'); ylabel('u(z)=\beta e\phi(z)');
    title(sprintf('Potential profile at u0_{max}=%.3g', u0_list(end)));
    grid on;

    figure('Color','w');
    plot(X, cplus_over_c0, 'LineWidth', 1.8); hold on;
    plot(X, cminus_over_c0, 'LineWidth', 1.8);
    yline(1,'--');
    xlabel('z/\lambda_D'); ylabel('c_{\pm}(z)/c_0');
    legend('c_+/c_0','c_-/c_0','bulk','Location','best');
    title(sprintf('Ion concentration profiles (linear) at u0_{max}=%.3g', u0_list(end)));
    grid on;

    % Warn if linear approximation becomes unphysical
    if any(cplus_over_c0 < 0) || any(cminus_over_c0 < 0)
        warning('Linear concentration went negative somewhere. Reduce |u0| or use nonlinear WDA+HS model.');
    end
end

%% ---- (Optional) Profiles at min voltage for symmetry check ----
if ~isempty(sol_at_u0min)
    X = linspace(par.XR, par.Xmax, 1000);
    Y = deval(sol_at_u0min, X);
    u = Y(1,:);

    aeff = par.alpha_eff;
    cplus_over_c0  = 1 - aeff*u;
    cminus_over_c0 = 1 + aeff*u;

    figure('Color','w');
    plot(X, u, 'LineWidth', 1.8);
    xlabel('z/\lambda_D'); ylabel('u(z)=\beta e\phi(z)');
    title(sprintf('Potential profile at u0_{min}=%.3g', u0_list(1)));
    grid on;

    figure('Color','w');
    plot(X, cplus_over_c0, 'LineWidth', 1.8); hold on;
    plot(X, cminus_over_c0, 'LineWidth', 1.8);
    yline(1,'--');
    xlabel('z/\lambda_D'); ylabel('c_{\pm}(z)/c_0');
    legend('c_+/c_0','c_-/c_0','bulk','Location','best');
    title(sprintf('Ion concentration profiles (linear) at u0_{min}=%.3g', u0_list(1)));
    grid on;
end

%% ============================================================
% --- ODE system: y = [u; uX; uXX; uXXX]
% delta_c^4 u'''' + (2*dc^2 - 1/alpha_eff) u'' + u = 0
% ============================================================
function dydx = odefun_bsk_alpha(~, y, par)
    dc   = par.delta_c;
    aeff = par.alpha_eff;

    if dc == 0
        error('delta_c=0 makes the 4th-order model singular.');
    end

    dydx = zeros(4,1);
    dydx(1) = y(2);
    dydx(2) = y(3);
    dydx(3) = y(4);

    coef = (2*dc^2 - 1/aeff);
    dydx(4) = ( -coef*y(3) - y(1) ) / (dc^4);
end

%% ============================================================
% --- Boundary conditions for potential scan:
% u(XR)=u0, u'''(XR)=0, u(Xmax)=0, u'(Xmax)=0
% ============================================================
function res = bcfun_potential_scan(ya, yb, par)
    res = zeros(4,1);
    res(1) = ya(1) - par.u0; % u(XR)=u0
    res(2) = ya(4);          % u'''(XR)=0
    res(3) = yb(1);          % u(Xmax)=0
    res(4) = yb(2);          % u'(Xmax)=0
end

%% ============================================================
% --- Initial guess for the first u0
% ============================================================
function y0 = guess_u(x, par)
    L = 2.0;         % decay length guess (dimensionless)
    A = par.u0;
    u    = A*exp(-(x-par.XR)/L);
    ux   = -(A/L)*exp(-(x-par.XR)/L);
    uxx  = (A/L^2)*exp(-(x-par.XR)/L);
    uxxx = -(A/L^3)*exp(-(x-par.XR)/L);
    y0 = [u; ux; uxx; uxxx];
end

%% ============================================================
% --- Continuation guess from previous solution (version-safe)
% Must return a vector for scalar x, and 4xN for vector x
% ============================================================
function y = guess_from_prev(x, sol_prev, du0)
    y = deval(sol_prev, x);   % 4xN
    y(1,:) = y(1,:) + du0;    % shift only u
end
