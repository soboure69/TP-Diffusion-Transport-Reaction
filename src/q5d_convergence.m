%% Q5d - Etude de convergence (cas transport dominant, sigma=0)
%  -mu u'' + b u' = 0, u(0)=0, u(1)=1
%  Exact: u(x) = (exp(b/mu x)-1)/(exp(b/mu)-1)
%  Erreurs L2, H1-semi et max nodale en fonction de h, pour
%  Galerkin, upwind optimal et SUPG optimal.

clear; close all;
if exist('OCTAVE_VERSION','builtin'), graphics_toolkit('fltk'); warning('off','all'); end

mu = 0.01; b_coef = 1; sigma = 0;
fhandle = @(x) 0;
alpha = b_coef/mu;
uexact  = @(x) (exp(alpha*x)-1)/(exp(alpha)-1);
duexact = @(x) (alpha*exp(alpha*x))/(exp(alpha)-1);

n_list = [10, 20, 40, 80, 160, 320];
methods = {'galerkin', 'upwind_opt', 'supg'};
mnames  = {'Galerkin', 'Upwind optimal', 'SUPG optimal'};

fprintf('\n=== Q5d: Convergence (transport dominant, mu=%.2f, b=%g) ===\n', mu, b_coef);

for m = 1:length(methods)
  fprintf('\n--- %s ---\n', mnames{m});
  fprintf('%6s %10s %10s %8s %12s %8s\n', 'n', 'h', 'eL2', 'ordre', 'emax', 'ordre');
  eL2_prev = NaN; emax_prev = NaN; h_prev = NaN;
  eL2_all = zeros(size(n_list));
  h_all   = zeros(size(n_list));
  for i = 1:length(n_list)
    n = n_list(i);
    h = 1/n;
    [xn, uh] = solve_dtr(mu, b_coef, sigma, n, fhandle, methods{m});
    [eL2, eH1, emax] = compute_error(xn, uh, uexact, duexact);
    eL2_all(i) = eL2; h_all(i) = h;
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
  h_data = h_all;
end

% log-log plot
figure();
loglog(h_data, loglog_data{1}, 'rs-', 'LineWidth', 1.2); hold on;
loglog(h_data, loglog_data{2}, 'g^-', 'LineWidth', 1.2);
loglog(h_data, loglog_data{3}, 'bd-', 'LineWidth', 1.2);
loglog(h_data, 0.5*h_data.^2, 'k--', 'LineWidth', 1);   % reference O(h^2)
xlabel('h'); ylabel('erreur L2');
title('Q5d: Convergence L2 (transport dominant)');
legend('Galerkin', 'Upwind optimal', 'SUPG optimal', 'pente 2 (ref)', ...
       'Location', 'NorthWest');
grid on;
print('-dpng', '-r150', '../figures/q5d_convergence.png');

fprintf('\nQ5d done.\n');
fprintf('Commentaire:\n');
fprintf('- Galerkin: pre-asymptotique tant que Pe_h>1 (erreur grande), puis O(h^2) en L2.\n');
fprintf('- Upwind optimal et SUPG: nodalement exacts -> erreur L2 = erreur d''interpolation O(h^2),\n');
fprintf('  precis des les maillages grossiers (pas d''oscillations).\n');
