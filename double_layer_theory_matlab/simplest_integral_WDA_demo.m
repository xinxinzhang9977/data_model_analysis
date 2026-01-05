function integral_WDA_concentration_only()
% Minimal demo of the INTEGRAL (convolution) idea in the paper:
% given an external potential phi(x),
% compute c+ and c- using weighted-density approximation.

clc; clear; close all;

%% constants
e0 = 1.602176634e-19;
kB = 1.380649e-23;
NA = 6.02214076e23;

%% parameters (paper-like)
T    = 300;
beta = 1/(kB*T);

d  = 0.5e-9;
R  = d/2;
v  = (4/3)*pi*R^3;

c0_M = 5;
c0   = c0_M*1000;
c0n  = c0*NA;

%% grid
dx = R/10;
L  = 60*R;
x  = (dx/2:dx:L-dx/2).';
N  = numel(x);

mask = (x >= R);

%% weighting kernels (1D, from SM)
mmax = round(R/dx);
offset = (-mmax:mmax)*dx;

ws = (abs(offset)<=R)/(2*R);
wv = zeros(size(offset));
inside = abs(offset)<=R;
wv(inside) = pi*(R^2-offset(inside).^2)/v;

% normalize discretely
ws = ws/(sum(ws)*dx);
wv = wv/(sum(wv)*dx);

%% bulk CS chemical potential
p_bulk = v*(2*c0n);
mu_ex_bulk = mu_ex_hat_CS(p_bulk);

%% ---- external potential (chosen, not solved) ----
% this is the KEY simplification
% ---- external oscillatory potential (effective EDL-like) ----
phi0 = 0.05;           % 50 mV
lambda = 2*d;          % decay length
kosc = 2*pi/d;         % oscillation ~ ion diameter

phi = phi0 * exp(-x/lambda) .* cos(kosc*x);

%% weighted potential
phi_bar = conv(phi, ws, 'same')*dx;
ubar = beta*e0*phi_bar;

%% initial concentrations
cp = c0n*ones(N,1);
cm = c0n*ones(N,1);
cp(~mask)=0; cm(~mask)=0;

%% simple Picard iteration for consistency of mu_ex
for it = 1:5000

    % weighted packing
    ctot = cp + cm;
    pbar = v * (conv(ctot, wv, 'same')*dx);

    mu_ex = mu_ex_hat_CS(pbar);
    mu_ex_bar = conv(mu_ex, wv, 'same')*dx;

    cp_new = c0n .* exp( -(+1)*ubar - mu_ex_bar + mu_ex_bulk );
    cm_new = c0n .* exp( -(-1)*ubar - mu_ex_bar + mu_ex_bulk );

    cp_new(~mask)=0; cm_new(~mask)=0;

    % mild relaxation
    cp = 0.7*cp + 0.3*cp_new;
    cm = 0.7*cm + 0.3*cm_new;
end

%% plots
figure;
plot(x/d, beta*e0*phi, 'LineWidth',1.8);
xlabel('x/d'); ylabel('\beta e\phi'); grid on;
title('Prescribed electrode potential');

figure;
plot(x/d, cp/c0n, 'LineWidth',1.8); hold on;
plot(x/d, cm/c0n, 'LineWidth',1.8);
xlabel('x/d'); ylabel('c_{\pm}/c_0'); grid on;
legend('c_+','c_-','Location','best');
title('Integral WDA concentration response');

end

function muhat = mu_ex_hat_CS(p)
muhat = (8*p - 9*p.^2 + 3*p.^3) ./ (1 - p).^3;
end
