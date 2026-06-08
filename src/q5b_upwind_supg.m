%% Q5b - Transport dominant : Galerkin vs Upwind (diffusion artificielle)
%  -mu u'' + b u' = 0, u(0)=0, u(1)=1   avec  mu=0.01, b=1
%  Mise en application DIRECTE des routines fournies (sans solve_dtr).
%
%  Upwind "plein" : on ajoute a la matrice de raideur la diffusion
%  artificielle rho_h = |b|*h/2.  La contribution elementaire ajoutee est
%      (rho_h/h) * [ 1 -1 ; -1 1 ]
%  (avec b=1 cela vaut 0.5*[1 -1;-1 1], independant de h).

clear; close all;
if exist('OCTAVE_VERSION','builtin'), graphics_toolkit('fltk'); warning('off','all'); end

%--- donnees physiques ---
mu = 0.01; b_coef = 1; sigma = 0;
acoef = -mu;  bcoef = b_coef;  ccoef = sigma;
nnel = 2; ndof = 1;

%--- solution exacte (analytique) en reference ---
x_exact = linspace(0,1,500)';
alpha = b_coef/mu;
u_exact = (exp(alpha*x_exact) - 1) / (exp(alpha) - 1);

n_values = [10, 40];
outfiles = {'../figures/q5b_upwind_supg.png', '../figures/q5b_upwind_supg_n40.png'};

for ic = 1:length(n_values)
  nel = n_values(ic);
  nnode = nel + 1;  sdof = nnode*ndof;
  hh = 1/nel;
  rho_h = abs(b_coef)*hh/2;        % diffusion artificielle (upwind plein)

  %--- coordonnees et connectivite ---
  gcoord = zeros(nnode,1);
  for k1 = 1:nnode, gcoord(k1) = (k1-1)*hh; end
  nodes = zeros(nel,2);
  for iel = 1:nel, nodes(iel,1)=iel; nodes(iel,2)=iel+1; end

  bcdof = [1, nnode];  bcval = [0, 1];

  %--- assemblage Galerkin et Upwind (routines fournies) ---
  kk_gal = zeros(sdof);  ff_gal = zeros(sdof,1);
  kk_upw = zeros(sdof);  ff_upw = zeros(sdof,1);
  Dart = [1 -1; -1 1];             % matrice de diffusion unitaire
  for iel = 1:nel
    xl = gcoord(nodes(iel,1));  xr = gcoord(nodes(iel,2));
    eleng = xr - xl;
    index = feeldof1(iel,nnel,ndof);
    ke = feode2l(acoef,bcoef,ccoef,eleng);     % Galerkin (fourni)
    fe = [0; 0];                               % f = 0
    [kk_gal,ff_gal] = feasmbl2(kk_gal,ff_gal,ke,fe,index);

    % upwind : on AJOUTE rho_h a la matrice de raideur
    ke_up = ke + (rho_h/eleng)*Dart;
    [kk_upw,ff_upw] = feasmbl2(kk_upw,ff_upw,ke_up,fe,index);
  end

  [kk_gal,ff_gal] = feaplyc2(kk_gal,ff_gal,bcdof,bcval);
  [kk_upw,ff_upw] = feaplyc2(kk_upw,ff_upw,bcdof,bcval);
  uh_gal = kk_gal \ ff_gal;
  uh_upw = kk_upw \ ff_upw;

  figure();
  plot(x_exact, u_exact, 'k-', 'LineWidth', 2); hold on;
  plot(gcoord, uh_gal, 'rs--', 'MarkerSize', 5, 'LineWidth', 1);
  plot(gcoord, uh_upw, 'bd-',  'MarkerSize', 5, 'LineWidth', 1);
  xlabel('x'); ylabel('u(x)');
  title(sprintf('Q5b: Galerkin vs Upwind (\\mu=%.2f, b=%g, n=%d, \\rho_h=%.4g)', ...
        mu, b_coef, nel, rho_h));
  legend('Exacte', 'Galerkin', 'Upwind (\rho_h=|b|h/2)', 'Location', 'NorthWest');
  grid on;
  print('-dpng', '-r150', outfiles{ic});

  fprintf('n=%d, h=%.4f, Pe_h=%.2f, rho_h=%.4f | Galerkin min=%.4f | Upwind min=%.4f\n', ...
          nel, hh, b_coef*hh/(2*mu), rho_h, min(uh_gal), min(uh_upw));
end

fprintf('\nQ5b done.\n');
fprintf('- Galerkin standard : oscillations spurieuses pour Pe_h > 1.\n');
fprintf('- Upwind (ajout de rho_h=|b|h/2 a la raideur) : oscillations supprimees (min>=0),\n');
fprintf('  mais schema diffusif (couche limite lissee).\n');
