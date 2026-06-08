%% Q6c - Etude de convergence (cas reaction dominant, b=0)
%  -mu u'' + sigma u = 0, u(0)=0, u(1)=1
%  Exact: u(x) = sinh(alpha x)/sinh(alpha), alpha = sqrt(sigma/mu)
%  Erreurs L2 et max nodale en fonction de h, pour
%  Galerkin (masse consistante) et Galerkin + mass-lumping.

clear; close all;
if exist('OCTAVE_VERSION','builtin'), graphics_toolkit('fltk'); warning('off','all'); end

mu = 1/6000; sigma = 1; b_coef = 0;
fhandle = @(x) 0;
alpha = sqrt(sigma/mu);
uexact  = @(x) sinh(alpha*x)/sinh(alpha);
duexact = @(x) alpha*cosh(alpha*x)/sinh(alpha);

n_list = [20, 40, 80, 160, 320, 640];
methods = {'galerkin', 'galerkin_lumped'};
mnames  = {'Galerkin (consistant)', 'Galerkin + mass-lumping'};

fprintf('\n=== Q6c: Convergence (reaction dominant, mu=1/%d, sigma=%g, alpha=%.1f) ===\n', ...
        round(1/mu), sigma, alpha);

h_data = 1 ./ n_list;
for m = 1:length(methods)
  fprintf('\n--- %s ---\n', mnames{m});
  fprintf('%6s %10s %10s %8s %12s %8s\n', 'n', 'h', 'eL2', 'ordre', 'emax', 'ordre');
  eL2_prev = NaN; emax_prev = NaN; h_prev = NaN;
  eL2_all = zeros(size(n_list));
  for i = 1:length(n_list)
    n = n_list(i);
    h = 1/n;
    [xn, uh] = solve_dtr(mu, b_coef, sigma, n, fhandle, methods{m});
    [eL2, eH1, emax] = compute_error(xn, uh, uexact, duexact);
    eL2_all(i) = eL2;
    if i == 1
      fprintf('%6d %10.4e %10.4e %8s %12.4e %8s\n', n, h, eL2, '-', emax, '-');
    else
      pL2 = log(eL2_prev/eL2)/log(h_prev/h);
      pmx = log(emax_prev/emax)/log(h_prev/h);
      fprintf('%6d %10.4e %10.4e %8.2f %12.4e %8.2f\n', n, h, eL2, pL2, emax, pmx);
    end
    eL2_prev = eL2; emax_prev = emax; h_prev = h;
  end
  loglog_data{m} = eL2_all;
end

figure();
loglog(h_data, loglog_data{1}, 'rs-', 'LineWidth', 1.2); hold on;
loglog(h_data, loglog_data{2}, 'g^-', 'LineWidth', 1.2);
loglog(h_data, h_data.^2, 'k--', 'LineWidth', 1);   % reference O(h^2)
xlabel('h'); ylabel('erreur L2');
title('Q6c: Convergence L2 (reaction dominant)');
legend('Galerkin', 'Galerkin + mass-lumping', 'pente 2 (ref)', 'Location', 'NorthWest');
grid on;
print('-dpng', '-r150', '../figures/q6c_convergence.png');

fprintf('\nQ6c done.\n');
fprintf('Commentaire:\n');
fprintf('- Les deux methodes convergent en O(h^2) en L2 une fois la couche limite resolue.\n');
fprintf('- Le mass-lumping garantit une solution sans oscillations (monotone) des les maillages\n');
fprintf('  grossiers, au prix d''une constante d''erreur legerement plus grande.\n');
