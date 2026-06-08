function [xnodes, uh] = solve_dtr(mu, b, sigma, n, fhandle, method, tau_param)

%-------------------------------------------------------------------
%  Solve -mu u'' + b u' + sigma u = f  on (0,1)
%  with  u(0) = 0,  u(1) = 1
%  using P1 (linear) finite elements.
%
%  Inputs:
%    mu, b, sigma - positive coefficients
%    n            - number of sub-intervals
%    fhandle      - function handle for f(x), e.g. @(x) 1 or @(x) x
%    method       - 'galerkin' | 'upwind_full' | 'upwind_opt' | 'supg' | 'gals'
%    tau_param    - (optional) vector of stabilization parameters per element
%                   if empty/missing, computed automatically per method
%
%  Outputs:
%    xnodes - nodal coordinates (n+1 x 1)
%    uh     - finite element solution at nodes (n+1 x 1)
%-------------------------------------------------------------------

if nargin < 7, tau_param = []; end

nnel  = 2;   % nodes per element
ndof  = 1;   % dofs per node
sdof  = n+1; % total system dofs
h     = 1/n;

% nodal coordinates
xnodes = linspace(0,1,n+1)';

% initialize system
kk = zeros(sdof, sdof);
ff = zeros(sdof, 1);

for iel = 1:n
  xl    = xnodes(iel);
  xr    = xnodes(iel+1);
  eleng = h;

  % DOF connectivity
  index = feeldof1(iel, nnel, ndof);

  % --- Base Galerkin element matrix ---
  k_el = feode2l(-mu, b, sigma, eleng);

  % --- Load vector ---
  f_el = fef_load(fhandle, xl, xr);

  % --- Stabilization ---
  Pe = abs(b)*eleng / (2*mu + eps);   % element Peclet number

  switch lower(method)
    case 'galerkin'
      % no stabilization

    case 'galerkin_lumped'
      % Reaction stabilization by mass lumping.
      % Lumped mass = diag row-sum of the consistent reaction mass matrix.
      % Equivalent to adding (sigma*h/6)*[1 -1; -1 1] to the Galerkin matrix:
      %   (sigma h/6)[2 1;1 2] + (sigma h/6)[1 -1;-1 1] = (sigma h/2) I
      k_el = k_el + (sigma*eleng/6)*[1 -1; -1 1];

    case 'upwind_full'
      % add full artificial diffusion mu_art = |b|*h/2
      mu_art = abs(b)*eleng/2;
      k_el = k_el + (mu_art/eleng)*[1 -1; -1 1];

    case 'upwind_opt'
      % optimal artificial diffusion (Scharfetter-Gummel / exact upwind)
      if Pe > 1e-6
        xi = coth(Pe) - 1/Pe;
      else
        xi = Pe/3;  % Taylor expansion for small Pe
      end
      mu_art = xi * abs(b)*eleng/2;
      k_el = k_el + (mu_art/eleng)*[1 -1; -1 1];

    case 'supg'
      % SUPG: perturbed test w_i = N_i + tau * b * N_i'
      if isempty(tau_param)
        if Pe > 1e-6
          xi = coth(Pe) - 1/Pe;
        else
          xi = Pe/3;
        end
        tau = xi * eleng / (2*abs(b) + eps);
      else
        tau = tau_param(iel);
      end
      % dN/dx = [-1/L, 1/L]
      dN = [-1/eleng, 1/eleng];
      % SUPG element stiffness addition: tau * b * dNi * (b*dNj + sigma*Nj)
      % = tau*b^2 * (dNi)(dNj) + tau*b*sigma * (dNi)(Nj integrated)
      % For the dNi*dNj term:
      k_supg = tau * b^2 * (1/eleng) * [1 -1; -1 1];
      % For the dNi*Nj term (convection-reaction cross):
      % \int dNi * Nj dx = dNi * \int Nj dx = dNi * [L/2]
      % Actually dNi(x) is constant, Nj(x) is linear:
      % \int N1 dx = L/2, \int N2 dx = L/2
      % dN1 = -1/L, dN2 = 1/L
      k_supg = k_supg + tau * b * sigma * [-1/2, -1/2; 1/2, 1/2];
      k_el = k_el + k_supg;
      % SUPG RHS addition: tau * b * dNi * f  (integrated)
      % \int tau * b * dNi * f(x) dx
      f_supg = [0; 0];
      gp = [-1/sqrt(3), 1/sqrt(3)];
      gw = [1, 1];
      for g = 1:2
        xi_g = gp(g);
        x_g  = 0.5*(1-xi_g)*xl + 0.5*(1+xi_g)*xr;
        jac  = eleng/2;
        fx   = fhandle(x_g);
        for i = 1:2
          f_supg(i) = f_supg(i) + gw(g) * tau * b * dN(i) * fx * jac;
        end
      end
      f_el = f_el + f_supg;

    case 'gals'
      % Galerkin Least Squares: perturbed test w_i = N_i + tau * L(N_i)
      % For P1: L(N_i) = b*N_i' + sigma*N_i  (since N_i'' = 0)
      if isempty(tau_param)
        % Codina-type tau valid for advection, reaction or diffusion dominated:
        %   tau = [ 4*mu/h^2 + 2*|b|/h + sigma ]^{-1}
        tau = 1 / (4*mu/eleng^2 + 2*abs(b)/eleng + sigma);
      else
        tau = tau_param(iel);
      end
      dN = [-1/eleng, 1/eleng];
      % GALS stiffness: tau * \int (b*dNi + sigma*Ni)(b*dNj + sigma*Nj) dx
      % = tau [b^2 (dNi dNj)*L + b*sigma*(dNi Nj + Ni dNj) integrals + sigma^2 mass]
      % b^2 dNidNj term
      k_gals = tau * b^2 * (1/eleng) * [1 -1; -1 1];
      % b*sigma cross terms: tau*b*sigma*[\int dNi*Nj + \int Ni*dNj]
      % \int dNi*Nj: row i col j = dNi * \int Nj = dNi * L/2
      G1 = [-1/2 -1/2; 1/2 1/2];
      G2 = G1';  % \int Ni*dNj
      k_gals = k_gals + tau*b*sigma*(G1 + G2);
      % sigma^2 mass term
      k_gals = k_gals + tau*sigma^2*(eleng/6)*[2 1; 1 2];
      k_el = k_el + k_gals;
      % GALS RHS: tau * \int (b*dNi + sigma*Ni) * f dx
      f_gals = [0; 0];
      gp = [-1/sqrt(3), 1/sqrt(3)];
      gw = [1, 1];
      for g = 1:2
        xi_g = gp(g);
        N1   = 0.5*(1-xi_g);
        N2   = 0.5*(1+xi_g);
        x_g  = N1*xl + N2*xr;
        jac  = eleng/2;
        fx   = fhandle(x_g);
        Nvals = [N1; N2];
        for i = 1:2
          f_gals(i) = f_gals(i) + gw(g)*(tau*b*dN(i) + tau*sigma*Nvals(i))*fx*jac;
        end
      end
      f_el = f_el + f_gals;

    otherwise
      error('Unknown method: %s', method);
  end

  % assemble
  [kk, ff] = feasmbl2(kk, ff, k_el, f_el, index);
end

% boundary conditions: u(0) = 0, u(1) = 1
bcdof = [1, sdof];
bcval = [0, 1];
[kk, ff] = feaplyc2(kk, ff, bcdof, bcval);

% solve
uh = kk \ ff;
