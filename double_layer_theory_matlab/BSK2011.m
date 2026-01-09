%% Bazant-Storey-Kornyshev (PRL 2011) 1D model: overscreening + crowding
% Solve: (1 - δc^2 d^2/dx^2) d^2φ/dx^2 = ρ(φ)
% i.e.  φ'' - δc^2 φ'''' = ρ(φ)
% with  ρ(φ) = sinh(φ) / (1 + 2γ sinh^2(φ/2))
%
% BCs (half-space truncated to [0,L]):
%   φ(0) = Vtilde
%   φ'''(0) = 0          (flat mean-field charge at surface)
%   φ(L) = 0
%   φ'(L) = 0
%
% After φ(x) is found, concentrations (relative to bulk) are:
%   c+/c0 = exp(-φ) / (1 + 2γ sinh^2(φ/2))
%   c-/c0 = exp(+φ) / (1 + 2γ sinh^2(φ/2))

clear; clc; close all;

%% ---- Parameters (you can tweak) ----
delta_c = 10;     % δc = lc / lambda_D  (PRL uses ~10 to match MD)
gamma   = 0.5;    % γ = 2 v c_ref  (bulk volume fraction parameter, 0<γ<1)
Vtilde  = 10;     % V~ = zeV/(kBT) (try 1, 10, 100)
L       = 30;     % domain size in x~ units (increase if decay not finished)

% Numerical controls
Nmesh   = 300;    % initial mesh points
opts = bvpset('RelTol',1e-6,'AbsTol',1e-8,'NMax',20000);

%% ---- Solve BVP using bvp4c ----
xmesh = linspace(0, L, Nmesh);
solinit = bvpinit(xmesh, @(x) guess_fun(x, Vtilde)); % initial guess
sol = bvp4c(@(x,y) odefun(x,y,delta_c,gamma), ...
            @(ya,yb) bcfun(ya,yb,Vtilde), solinit, opts);

% Evaluate solution on a fine grid
x = linspace(0, L, 1200);
Y = deval(sol, x);

phi  = Y(1,:);  % φ
dphi = Y(2,:);  % φ'
ddphi= Y(3,:);  % φ''
dddphi=Y(4,:);  % φ'''

%% ---- Post-processing: charge density and concentrations ----
den = 1 + 2*gamma*sinh(phi/2).^2;
rho = sinh(phi)./den;                       % dimensionless charge source term (Eq. 6 RHS)

cplus_rel  = exp(-phi)./den;                % c+/c0
cminus_rel = exp(+phi)./den;                % c-/c0
ctot_rel   = cplus_rel + cminus_rel;        % (c+ + c-)/c0

%% ---- Plots ----
figure;
plot(x, phi, 'LineWidth', 1.5);
xlabel('x~  (distance / \lambda_D)');
ylabel('\phi~');
title(sprintf('Potential profile: \\delta_c=%.2g, \\gamma=%.2g, V~=%.2g', delta_c, gamma, Vtilde));
grid on;

figure;
plot(x, cplus_rel, 'LineWidth', 1.5); hold on;
plot(x, cminus_rel,'LineWidth', 1.5);
plot(x, ctot_rel, 'LineWidth', 1.5);
xlabel('x~  (distance / \lambda_D)');
ylabel('c/c_0  (relative to bulk)');
legend('c_+ / c_0','c_- / c_0','(c_+ + c_-)/c_0','Location','best');
title('Ion concentration profiles (lattice-gas / Bikerman crowding)');
grid on;

figure;
plot(x, rho, 'LineWidth', 1.5);
xlabel('x~  (distance / \lambda_D)');
ylabel('\rho~(\phi~)');
title('Dimensionless charge density source term');
grid on;

%% ---- Print a quick sanity check at far field ----
fprintf('At x=L: phi=%.3e, phi''=%.3e, c+/c0=%.6f, c-/c0=%.6f\n', ...
    phi(end), dphi(end), cplus_rel(end), cminus_rel(end));

%% ================== Local functions ==================

function dy = odefun(~, y, delta_c, gamma)
% State vector:
% y(1)=phi, y(2)=phi', y(3)=phi'', y(4)=phi'''
phi = y(1);
dphi = y(2);
ddphi = y(3);
dddphi = y(4);

den = 1 + 2*gamma*sinh(phi/2)^2;
rho = sinh(phi)/den;  % RHS of Eq. (6) in 1D

% Governing ODE: phi'' - delta_c^2 * phi'''' = rho
% => phi'''' = (phi'' - rho)/delta_c^2
phi4 = (ddphi - rho)/(delta_c^2);

dy = zeros(4,1);
dy(1) = dphi;
dy(2) = ddphi;
dy(3) = dddphi;
dy(4) = phi4;
end

function res = bcfun(ya, yb, Vtilde)
% Boundary conditions:
% at x=0: phi = Vtilde, phi''' = 0
% at x=L: phi = 0, phi' = 0
res = zeros(4,1);
res(1) = ya(1) - Vtilde;  % phi(0)=V~
res(2) = ya(4);           % phi'''(0)=0
res(3) = yb(1);           % phi(L)=0
res(4) = yb(2);           % phi'(L)=0
end

function yguess = guess_fun(x, Vtilde)
% A simple decaying exponential initial guess
phi  = Vtilde*exp(-x);
dphi = -Vtilde*exp(-x);
ddphi=  Vtilde*exp(-x);
dddphi=-Vtilde*exp(-x);
yguess = [phi; dphi; ddphi; dddphi];
end
