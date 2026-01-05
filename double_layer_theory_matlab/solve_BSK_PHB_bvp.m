function sol = solve_BSK_PHB_bvp()
% Solve the dimensionless 1D BSK + hydration (PHB) model between two plates
% using bvp4c.
%
% Unknowns (dimensionless):
%   y(1)=phi, y(2)=phi', y(3)=phi'', y(4)=phi''', y(5)=psi_h, y(6)=psi_h'
%
% Governing equations:
%   delta_c^2 * phi'''' - phi'' = c_plus - c_minus
%   c_plus  = exp(-z*phi + psi_h)
%   c_minus = exp(phi)
%
%   psi_h'' - kappa2 * psi_h = -alphaHyd * kappa2 * (c_plus - 1)
%   => psi_h'' = kappa2*psi_h - alphaHyd*kappa2*(c_plus-1)
%
% Boundary conditions at xL, xR (dimensionless):
%   phi' - delta_c^2 * phi''' = -S   at xL
%   phi''' = -phi''                 at xL
%   phi' - delta_c^2 * phi''' = +S   at xR
%   phi''' = +phi''                 at xR
%   psi_h' = -B                      at xL
%   psi_h' = +B                      at xR

  % -------------------------
  % Parameters (EDIT HERE)
  % -------------------------
  p.z        = 2;      % counterion valence for z:1 electrolyte (e.g. z=2 for 2:1)
  p.delta_c  = 0.5;    % delta_c = l_c / lambda_D
  p.kappa2   = (2.0)^2; % (kappa_h * lambda_D)^2  (so kappa2 >= 0)
  p.alphaHyd = 1.4;    % hydration coupling prefactor (set from paper's Eq.(14))
  p.S        = 5.0;    % S = (q_s e lambda_D)/(eps kBT)  (dimensionless surface charge param)
  p.B        = 0.2;    % B = RHS magnitude for psi_h' BC (dimensionless, from paper's BCs)

  % Domain endpoints (dimensionless)
  % xL = -(d/(2*lambda_D) + l_h/lambda_D), xR = +(d/(2*lambda_D) + l_h/lambda_D)
  p.d_over_lambdaD  = 10.0;   % d / lambda_D
  p.lh_over_lambdaD = 0.5;   % l_h / lambda_D
  p.xL = -0.5*p.d_over_lambdaD - p.lh_over_lambdaD;
  p.xR =  0.5*p.d_over_lambdaD + p.lh_over_lambdaD;

  % -------------------------
  % Initial mesh + guess
  % -------------------------
  xmesh = linspace(p.xL, p.xR, 200);

  % A simple initial guess (weak potential)
  solinit = bvpinit(xmesh, @(x) guessFcn(x, p));

  % -------------------------
  % Solve BVP
  % -------------------------
  opts = bvpset('RelTol',1e-6,'AbsTol',1e-8,'NMax',50000);
  sol  = bvp4c(@(x,y) odefun(x,y,p), @(ya,yb) bcfun(ya,yb,p), solinit, opts);

  % -------------------------
  % Plot results
  % -------------------------
  x   = sol.x;
  Y   = sol.y;
  phi = Y(1,:);
  psi = Y(5,:);

  % Concentrations
  cplus  = exp(-p.z*phi + psi);
  cminus = exp(phi);

  figure; plot(x, phi, 'LineWidth', 1.5); grid on;
  xlabel('\tilde{x}'); ylabel('\tilde{\phi}');
  title('Dimensionless potential');

  figure; plot(x, psi, 'LineWidth', 1.5); grid on;
  xlabel('\tilde{x}'); ylabel('\psi_h');
  title('Hydration potential');

  figure; plot(x, cplus, 'LineWidth', 1.5); hold on;
  plot(x, cminus, 'LineWidth', 1.5); grid on;
  xlabel('\tilde{x}'); ylabel('c / c_0 (dimensionless)');
  legend('c_+','c_-');
  title('Ion profiles');

end

% =============== ODE system ===============
function dydx = odefun(x, y, p) %#ok<INUSD>
  phi    = y(1);
  dphi   = y(2);
  ddphi  = y(3);
  dddphi = y(4);

  psi    = y(5);
  dpsi   = y(6);

  % Boltzmann distributions (dimensionless)
  cplus  = exp(-p.z*phi + psi);
  cminus = exp(phi);

  % From: delta_c^2*phi'''' - phi'' = cplus - cminus
  % => phi'''' = (phi'' + (cplus - cminus)) / delta_c^2
  phi4   = (ddphi + (cplus - cminus)) / (p.delta_c^2);

  % Hydration Helmholtz:
  % psi'' - kappa2*psi = -alphaHyd*kappa2*(cplus-1)
  % => psi'' = kappa2*psi - alphaHyd*kappa2*(cplus-1)
  psi2   = p.kappa2*psi - p.alphaHyd*p.kappa2*(cplus - 1);

  dydx = zeros(6,1);
  dydx(1) = dphi;
  dydx(2) = ddphi;
  dydx(3) = dddphi;
  dydx(4) = phi4;
  dydx(5) = dpsi;
  dydx(6) = psi2;
end

% =============== Boundary conditions ===============
function res = bcfun(ya, yb, p)
  % ya = y(xL), yb = y(xR)
  % y = [phi, phi', phi'', phi''', psi, psi']'

  res = zeros(6,1);

  % Left boundary xL
  res(1) = ya(2) - p.delta_c^2 * ya(4) + p.S;     % phi' - delta_c^2*phi''' = -S  -> bring to 0
  res(2) = ya(4) + ya(3);                         % phi''' = -phi''  -> phi'''+phi''=0
  res(3) = ya(6) + p.B;                            % psi' = -B -> psi' + B = 0

  % Right boundary xR
  res(4) = yb(2) - p.delta_c^2 * yb(4) - p.S;     % phi' - delta_c^2*phi''' = +S
  res(5) = yb(4) - yb(3);                         % phi''' = +phi'' -> phi'''-phi''=0
  res(6) = yb(6) - p.B;                            % psi' = +B

end

% =============== Initial guess function ===============
function y0 = guessFcn(x, p) %#ok<INUSD>
  % A mild guess: small antisymmetric slope in phi, psi ~ 0
  % You can replace this with something smarter if convergence is hard.
  y0 = zeros(6,1);
  y0(1) = 0.0;   % phi
  y0(2) = 0.0;   % phi'
  y0(3) = 0.0;   % phi''
  y0(4) = 0.0;   % phi'''
  y0(5) = 0.0;   % psi
  y0(6) = 0.0;   % psi'
end
