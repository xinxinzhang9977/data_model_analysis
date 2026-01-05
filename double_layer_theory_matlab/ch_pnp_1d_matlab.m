function ch_pnp_1d_matlab()
% CH–PNP (Onsager) in 1D, MATLAB template
% Variables: phase-field phi (0..1), cP, cM, potential psi
% BC: no-flux for phi and ions; electrical insulation (Neumann) for psi.
%
% Author: (template)

%% -------------------- parameters --------------------
kB  = 1.380649e-23;      % J/K
T   = 298.15;            % K
e0  = 1.602176634e-19;   % C

% Domain
L  = 1e-6;               % m
Nx = 400;
x  = linspace(0, L, Nx).';
dx = x(2)-x(1);

% Time
dt    = 2e-6;            % s  (adjust for stability; ions explicit-ish)
tEnd  = 5e-3;            % s
nStep = ceil(tEnd/dt);

% Phase-field (Cahn-Hilliard)
W        = 2e7;          % J/m^3  (double-well barrier strength)
kappaPhi = 2e-12;        % J/m    (gradient coeff; sets interface thickness)
Mphi     = 5e-16;        % (m^5)/(J*s)  mobility (CH)

% Electrolyte dielectric (two phases)
eps0 = 8.8541878128e-12; % F/m
epsA = 40*eps0;
epsB = 10*eps0;

% Ions
Dplus  = 2e-10;          % m^2/s
Dminus = 2e-10;          % m^2/s
beta   = 1/(kB*T);

% Solvation preference (partition) term: g_i(phi)=DeltaMu_i*phi
DeltaMu_plus  = 2.0*kB*T;   % J (per particle); tune for partition
DeltaMu_minus = 0.0*kB*T;

% Fixed charge (optional)
rho_f = zeros(Nx,1);     % C/m^3

% Numerical safety
cFloor = 1e-6;           % mol/m^3 or number density floor (consistent units!)
% NOTE: In this template, c is number density [1/m^3] for consistency with Poisson.
% If you prefer mol/m^3, convert with NA and e accordingly.

%% -------------------- initial conditions --------------------
% Phase-field: one interface around L/2
phi = 0.5*(1 - tanh((x - 0.5*L)/(0.03*L)));  % smooth step

% Ions: start uniform (number density). Example: 1000 mol/m^3 ~ 6e26 1/m^3
NA  = 6.02214076e23;
c0_molm3 = 1000;         % 1 M ~ 1000 mol/m^3
c0 = c0_molm3*NA;        % 1/m^3

cP = c0*ones(Nx,1);
cM = c0*ones(Nx,1);

% Potential
psi = zeros(Nx,1);

%% -------------------- operators (Neumann Laplacian) --------------------
D2 = laplacian_neumann_1d(Nx, dx);   % second derivative matrix
D4 = D2*D2;                          % biharmonic operator (Neumann-compatible)

I  = speye(Nx);

% Precompute CH implicit left-hand matrix (constant Mphi,kappaPhi here)
A_phi = (I - dt*Mphi*kappaPhi*D4);

%% -------------------- time loop --------------------
for n = 1:nStep
    t = n*dt;

    % ----- dielectric field -----
    eps = epsA*phi + epsB*(1-phi);

    % ----- solve variable-coefficient Poisson: -d/dx(eps dpsi/dx) = rho -----
    rho = e0*(cP - cM) + rho_f;  % C/m^3
    psi = solve_poisson_var_eps_neumann(eps, rho, dx);

    % ----- compute CH chemical potential pieces (explicit except laplacian term) -----
    % f_mix = W*phi^2(1-phi)^2  => f'(phi)=2W*phi(1-phi)(1-2phi)
    fprime = 2*W*phi.*(1-phi).*(1-2*phi);

    % coupling from solvation: sum_i c_i d g_i/dphi = cP*DeltaMu_plus + cM*DeltaMu_minus
    coupling_sol = cP*DeltaMu_plus + cM*DeltaMu_minus;

    % dielectric coupling: + 0.5 * eps'(phi) * |dpsi/dx|^2, eps' = epsA-epsB
    dpsi = grad_central_1d(psi, dx);
    coupling_die = 0.5*(epsA - epsB)*(dpsi.^2);

    % RHS for phi update:
    RHS_phi = phi + dt*Mphi * (D2*(fprime + coupling_sol + coupling_die));

    % Semi-implicit solve for phi^{n+1}
    phi_new = A_phi \ RHS_phi;

    % Optional: clamp (not strictly variational, but prevents overshoot early on)
    phi = min(max(phi_new, 0), 1);

    % ----- ion chemical potentials (no gradient penalty for ions in this minimal model) -----
    % mu_i = kBT ln(c_i/c_ref) + DeltaMu_i*phi +/- e psi
    % Choose c_ref = c0 for nondimensional log
    cP = max(cP, cFloor);
    cM = max(cM, cFloor);

    muP = (1/beta)*log(cP/c0) + DeltaMu_plus*phi + e0*psi;
    muM = (1/beta)*log(cM/c0) + DeltaMu_minus*phi - e0*psi;

    % ----- ion fluxes (conservative, face-based) -----
    % Onsager: J = -M c d(mu)/dx, with M = D/(kBT)
    MP = Dplus*beta;
    MM = Dminus*beta;

    JP = flux_onsager_upwind(cP, muP, MP, dx);
    JM = flux_onsager_upwind(cM, muM, MM, dx);

    % Update concentrations: dc/dt = -dJ/dx (Neumann no-flux => J=0 at boundaries)
    cP = cP - dt*div_face_flux_1d(JP, dx);
    cM = cM - dt*div_face_flux_1d(JM, dx);

    % positivity floor
    cP = max(cP, cFloor);
    cM = max(cM, cFloor);

    % ----- simple diagnostics -----
    if mod(n, 200) == 0 || n==1
        Q = trapz(x, rho);  % net charge (C/m^2 in 1D per unit area)
        fprintf('step %d/%d, t=%.3e s, min/max phi=%.3f/%.3f, min c=%.3e, netQ=%.3e\n',...
            n, nStep, t, min(phi), max(phi), min([cP;cM]), Q);

        % quick plot
        clf;
        subplot(3,1,1); plot(x, phi); ylabel('\phi'); ylim([-0.05 1.05]);
        subplot(3,1,2); plot(x, cP/NA, x, cM/NA); ylabel('c (mol/m^3)'); legend('c_+','c_-');
        subplot(3,1,3); plot(x, psi); ylabel('\psi (V)'); xlabel('x (m)');
        drawnow;
    end
end

end

%% ==================== helper functions ====================

function D2 = laplacian_neumann_1d(N, dx)
% Second derivative with Neumann BC (ghost = mirror)
% u_x=0 at boundaries => u0=u2, u_{N+1}=u_{N-1}
e = ones(N,1);
D2 = spdiags([e -2*e e], [-1 0 1], N, N) / dx^2;

% Modify first and last row for Neumann
% u_x(1)=0 -> u0=u2 => u_xx(1) = (u2 -2u1 + u0)/dx^2 = (2u2 -2u1)/dx^2
D2(1,1) = -2/dx^2; D2(1,2) =  2/dx^2;
% u_x(N)=0 -> u_{N+1}=u_{N-1} => u_xx(N)=(2u_{N-1}-2u_N)/dx^2
D2(N,N) = -2/dx^2; D2(N,N-1) = 2/dx^2;
end

function g = grad_central_1d(u, dx)
% Central gradient with Neumann at boundaries (du/dx=0)
N = numel(u);
g = zeros(N,1);
g(2:N-1) = (u(3:N) - u(1:N-2))/(2*dx);
g(1) = 0;
g(N) = 0;
end

function psi = solve_poisson_var_eps_neumann(eps, rho, dx)
% Solve: -d/dx( eps dpsi/dx ) = rho, with Neumann BC: (eps dpsi/dx)|bdry = 0
% Gauge fixing: set mean(psi)=0 by adding a constraint (replace one equation).

N = numel(eps);

% face eps
eps_f = zeros(N+1,1);
eps_f(2:N) = 0.5*(eps(1:N-1) + eps(2:N));
eps_f(1)   = eps(1);    % boundary faces
eps_f(N+1) = eps(N);

% Build matrix A * psi = b
% Discretization: -(1/dx)[ eps_{i+1/2}(psi_{i+1}-psi_i)/dx - eps_{i-1/2}(psi_i-psi_{i-1})/dx ] = rho_i
main = zeros(N,1);
low  = zeros(N-1,1);
up   = zeros(N-1,1);

for i = 1:N
    em = eps_f(i);     % i-1/2
    ep = eps_f(i+1);   % i+1/2
    main(i) = (em + ep)/dx^2;
    if i>1
        low(i-1) = -em/dx^2;
    end
    if i<N
        up(i) = -ep/dx^2;
    end
end

A = spdiags([ [low;0] main [0;up] ], [-1 0 1], N, N);

b = rho;  % units: C/m^3 ; then psi in V since eps in F/m

% Neumann Poisson is singular: fix gauge by enforcing psi(1)=0 or mean(psi)=0.
% Here: enforce mean(psi)=0 by replacing first row with ones.
A(1,:) = 1;
b(1)   = 0;

psi = A \ b;
end

function J = flux_onsager_upwind(c, mu, M, dx)
% Face flux J_{i+1/2} = - M * c_face * (mu_{i+1}-mu_i)/dx
% Upwind choice for c_face based on sign of -dmu/dx (drift direction):
% If (mu_{i+1}-mu_i) > 0 => dmu/dx positive => flux negative (to left), use right cell? 
% We'll use sign of velocity v = -M*(mu_{i+1}-mu_i)/dx; if v>0 (to right), take left c; else take right c.

N = numel(c);
dmu = diff(mu);                 % size N-1
v   = -M * (dmu/dx);            % "drift" velocity proportional

cL = c(1:N-1);
cR = c(2:N);

c_face = cL;
c_face(v<0) = cR(v<0);          % upwind

Jint = -M .* c_face .* (dmu/dx); % size N-1

% include boundary faces (no-flux): J(1)=J(N+1)=0
J = zeros(N+1,1);
J(2:N) = Jint;
end

function divJ = div_face_flux_1d(Jface, dx)
% divergence at cell centers from face fluxes:
% div = (J_{i+1/2}-J_{i-1/2})/dx
N = numel(Jface)-1;
divJ = (Jface(2:N+1) - Jface(1:N))/dx;
end
