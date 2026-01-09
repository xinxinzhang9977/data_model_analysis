function sol = bsk_goodwin_unified_demo()
  % -----------------------
  % Parameters (edit here)
  % -----------------------
  p.delta_c = 15;     % δc = lc/lambdaD
  p.gamma   = 0.9;     % bulk volume fraction (0<gamma<1)
  p.ap      = 0.1;     % a_+
  p.am      = 0.1;     % a_-
  p.b       = -0.9;     % b  (cross short-range)
  p.U       = 2.0;     % surface potential u(0)=zeV/kT (dimensionless)
  p.L       = 20;      % domain length in x~ units (large enough -> bulk)
  
  % initial mesh + guess
  xmesh = linspace(0, p.L, 200);
  yguess = @(x) [p.U*(1-x/p.L); -p.U/p.L; 0; 0];  % simple linear guess
  
  solinit = bvpinit(xmesh, yguess);
  opts = bvpset('RelTol',1e-6,'AbsTol',1e-8);
  sol = bvp4c(@(x,y) odefun(x,y,p), @(ya,yb) bcfun(ya,yb,p), solinit, opts);

  % quick plot
  figure; plot(sol.x, sol.y(1,:), 'LineWidth', 1.5);
  xlabel('x~'); ylabel('u = ze\phi/kT'); grid on;
  title('Unified BSK + Goodwin model');
end

function dydx = odefun(x, y, p)
  u  = y(1);
  up = y(2);
  u2 = y(3);
  u3 = y(4);

  % solve local composition (phi+, phi-) from algebraic relations
  [phip, phim] = local_phi_from_u(u, p);

  charge = (phip - phim);  % = phi+ - phi-

  % 4th-order equation: u'' - δc^2 u'''' = -charge
  % -> u'''' = (u'' + charge)/δc^2
  u4 = (u2 + charge) / (p.delta_c^2);

  dydx = [up; u2; u3; u4];
end

function res = bcfun(ya, yb, p)
  % ya: x=0, yb: x=L
  % BC: u(0)=U, u'''(0)=0, u(L)=0, u'(L)=0
  res = [
    ya(1) - p.U
    ya(4)       % u'''(0)=0
    yb(1)       % u(L)=0
    yb(2)       % u'(L)=0
  ];
end

function [phip, phim] = local_phi_from_u(u, p)
  % Solve (R+)(R-) for phip, phim given u
  % Unknowns: phip, phim in (0,1), with phi0=1-phip-phim>0

  phib = p.gamma/2;
  phi0b = 1 - p.gamma;
  if phi0b <= 0
    error('gamma must be < 1.');
  end

  Cplus  = log(phib/phi0b) + p.ap*phib - p.b*phib;
  Cminus = log(phib/phi0b) + p.am*phib - p.b*phib;

  % Newton initial guess: start from bulk, bias by u
  phip = min(max(phib*exp(-0.5*u), 1e-12), 1-1e-12);
  phim = min(max(phib*exp( 0.5*u), 1e-12), 1-1e-12);

  for it = 1:50
    phi0 = 1 - phip - phim;
    if phi0 <= 1e-12
      % push back into feasible region
      s = (1-1e-6)/(phip+phim);
      phip = phip*s; phim = phim*s;
      phi0 = 1 - phip - phim;
    end

    F1 = log(phip/phi0) + p.ap*phip - p.b*phim + u - Cplus;
    F2 = log(phim/phi0) + p.am*phim - p.b*phip - u - Cminus;

    % Jacobian
    dF1_dphip = 1/phip + 1/phi0 + p.ap;
    dF1_dphim = 1/phi0 - p.b;

    dF2_dphip = 1/phi0 - p.b;
    dF2_dphim = 1/phim + 1/phi0 + p.am;

    J = [dF1_dphip, dF1_dphim;
         dF2_dphip, dF2_dphim];

    d = -J\[F1; F2];

    phip_new = phip + d(1);
    phim_new = phim + d(2);

    % damp + clamp
    alpha = 1.0;
    for ls = 1:10
      pp = phip + alpha*d(1);
      pm = phim + alpha*d(2);
      if pp>1e-12 && pm>1e-12 && (1-pp-pm)>1e-12
        phip_new = pp; phim_new = pm;
        break;
      end
      alpha = alpha/2;
    end

    if max(abs([phip_new-phip, phim_new-phim])) < 1e-12
      phip = phip_new; phim = phim_new;
      return;
    end
    phip = phip_new; phim = phim_new;
  end

  % If not converged, still return last iterate (bvp4c may still converge globally)
end
