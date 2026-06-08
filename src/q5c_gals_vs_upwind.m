%% Q5c - Compare GALS (tau = h/(2b)) with upwind
%  -mu u'' + b u' = 0, u(0)=0, u(1)=1
%  b=1, mu=0.01

clear; close all;

mu = 0.01; b_coef = 1; sigma = 0;
fhandle = @(x) 0;

for n = [10, 40, 60]
  h = 1/n;
  % GALS with tau_j = h/(2b)
  tau_vec = ones(n,1) * h / (2*b_coef);
  [xn, uh_gals] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'gals', tau_vec);

  % Upwind full
  [~, uh_upf] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'upwind_full');

  % Upwind optimal
  [~, uh_upo] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'upwind_opt');

  % exact
  x_exact = linspace(0,1,500)';
  alpha   = b_coef/mu;
  u_exact = (exp(alpha*x_exact)-1)/(exp(alpha)-1);

  figure();
  plot(x_exact, u_exact, 'k-', 'LineWidth', 2); hold on;
  plot(xn, uh_gals, 'rs-', 'MarkerSize', 5, 'LineWidth', 1.2);
  plot(xn, uh_upf,  'bd-', 'MarkerSize', 5, 'LineWidth', 1.2);
  plot(xn, uh_upo,  'g^-', 'MarkerSize', 5, 'LineWidth', 1.2);
  xlabel('x'); ylabel('u(x)');
  title(sprintf('Q5c: GALS vs Upwind (\\mu=%.2f, b=%g, n=%d)', mu, b_coef, n));
  legend('Exacte', sprintf('GALS (\\tau=h/(2b)=%.4f)',h/(2*b_coef)), ...
         'Upwind (plein)', 'Upwind (optimal)', 'Location', 'NorthWest');
  grid on;
  print('-dpng', '-r150', sprintf('../figures/q5c_gals_upwind_n%d.png', n));
end

fprintf('Q5c done.\n');
fprintf('Commentaire:\n');
fprintf('- Pour sigma=0, GALS avec tau=h/(2b) est equivalent a upwind plein\n');
fprintf('  car le terme GALS tau*b^2*int(dNi*dNj) = (h/(2b))*b^2/h = b/2 [1 -1;-1 1]\n');
fprintf('  ce qui ajoute une viscosite artificielle mu_art = b*h/2 (identique a upwind plein).\n');
fprintf('- Upwind optimal (Scharfetter-Gummel) est plus precis car tau adaptatif.\n');
