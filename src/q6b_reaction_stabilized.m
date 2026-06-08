%% Q6b - Reaction dominant: methodes stabilisees
%  -mu u'' + sigma u = 0, b=0, u(0)=0, u(1)=1
%  mu=1/6000, sigma=1
%  Compare: Galerkin, Galerkin+mass-lumping, GALS, SUPG
%  Remarque: pour b=0, SUPG est inactif (terme convectif nul) -> SUPG == Galerkin

clear; close all;
if exist('OCTAVE_VERSION','builtin'), graphics_toolkit('fltk'); warning('off','all'); end

mu = 1/6000; sigma = 1; b_coef = 0;
fhandle = @(x) 0;
alpha = sqrt(sigma/mu);

x_exact = linspace(0,1,1000)';
u_exact = sinh(alpha*x_exact)/sinh(alpha);

for n = [10, 20]
  [xn, ug] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'galerkin');
  [~,  ul] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'galerkin_lumped');
  [~,  ua] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'gals');
  [~,  us] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'supg');

  figure();
  plot(x_exact, u_exact, 'k-', 'LineWidth', 2); hold on;
  plot(xn, ug, 'rs--', 'MarkerSize', 5, 'LineWidth', 1);
  plot(xn, ul, 'g^-',  'MarkerSize', 5, 'LineWidth', 1.2);
  plot(xn, ua, 'bd--', 'MarkerSize', 5, 'LineWidth', 1);
  plot(xn, us, 'mv:',  'MarkerSize', 6, 'LineWidth', 1);
  xlabel('x'); ylabel('u(x)');
  title(sprintf('Q6b: Reaction stabilisee (\\mu=1/%d, \\sigma=%g, n=%d)', round(1/mu), sigma, n));
  legend('Exacte', 'Galerkin', 'Galerkin + mass-lumping', 'GALS', 'SUPG (=Galerkin)', ...
         'Location', 'NorthWest');
  grid on;
  print('-dpng', '-r150', sprintf('../figures/q6b_reaction_stabilized_n%d.png', n));

  fprintf('n=%2d | min: Galerkin=%+.4f  Lumped=%+.4f  GALS=%+.4f  SUPG=%+.4f\n', ...
          n, min(ug), min(ul), min(ua), min(us));
end

fprintf('\nQ6b done.\n');
fprintf('Commentaire:\n');
fprintf('- SUPG est INACTIF pour b=0 (perturbation tau*b*v'' nulle) => identique a Galerkin.\n');
fprintf('- Le mass-lumping supprime totalement les oscillations (min(uh) >= 0) :\n');
fprintf('  il revient a ajouter (sigma*h/6)[1 -1;-1 1], une diffusion artificielle sigma*h^2/6.\n');
fprintf('- GALS avec tau de Codina N''ENLEVE PAS les oscillations en regime reaction-dominant\n');
fprintf('  (il ajoute une masse tau*sigma^2 positive). Le remede adapte ici est le mass-lumping.\n');
