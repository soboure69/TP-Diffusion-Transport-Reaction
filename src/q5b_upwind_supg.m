%% Q5b - Transport dominant: Galerkin vs upwind vs SUPG
%  -mu u'' + b u' = 0, u(0)=0, u(1)=1
%  b=1, mu=0.01
%  Compare exact, Galerkin, upwind (full), upwind (optimal), SUPG (optimal)

clear; close all;

mu = 0.01; b_coef = 1; sigma = 0;
fhandle = @(x) 0;
n = 10;  % low n to show oscillation differences clearly

% exact solution
x_exact = linspace(0,1,500)';
alpha = b_coef/mu;
u_exact = (exp(alpha*x_exact) - 1) / (exp(alpha) - 1);

% Galerkin
[xn, uh_gal] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'galerkin');

% Upwind (full artificial diffusion)
[~, uh_upf] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'upwind_full');

% Upwind (optimal / Scharfetter-Gummel)
[~, uh_upo] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'upwind_opt');

% SUPG (optimal tau)
[~, uh_supg] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'supg');

figure();
plot(x_exact, u_exact, 'k-', 'LineWidth', 2); hold on;
plot(xn, uh_gal,  'rs--', 'MarkerSize', 5, 'LineWidth', 1);
plot(xn, uh_upf,  'bd--', 'MarkerSize', 5, 'LineWidth', 1);
plot(xn, uh_upo,  'g^--', 'MarkerSize', 5, 'LineWidth', 1);
plot(xn, uh_supg, 'mv--', 'MarkerSize', 5, 'LineWidth', 1);
xlabel('x'); ylabel('u(x)');
title(sprintf('Q5b: Methodes stabilisees (\\mu=%.2f, b=%g, n=%d)', mu, b_coef, n));
legend('Exacte', 'Galerkin', 'Upwind (plein)', 'Upwind (optimal)', 'SUPG (optimal)', ...
       'Location', 'NorthWest');
grid on;
print('-dpng', '-r150', '../figures/q5b_upwind_supg.png');

% same with n=40
n = 40;
[xn, uh_gal]  = solve_dtr(mu, b_coef, sigma, n, fhandle, 'galerkin');
[~, uh_upf]   = solve_dtr(mu, b_coef, sigma, n, fhandle, 'upwind_full');
[~, uh_upo]   = solve_dtr(mu, b_coef, sigma, n, fhandle, 'upwind_opt');
[~, uh_supg]  = solve_dtr(mu, b_coef, sigma, n, fhandle, 'supg');

figure();
plot(x_exact, u_exact, 'k-', 'LineWidth', 2); hold on;
plot(xn, uh_gal,  'rs--', 'MarkerSize', 4, 'LineWidth', 1);
plot(xn, uh_upf,  'bd--', 'MarkerSize', 4, 'LineWidth', 1);
plot(xn, uh_upo,  'g^--', 'MarkerSize', 4, 'LineWidth', 1);
plot(xn, uh_supg, 'mv--', 'MarkerSize', 4, 'LineWidth', 1);
xlabel('x'); ylabel('u(x)');
title(sprintf('Q5b: Methodes stabilisees (\\mu=%.2f, b=%g, n=%d)', mu, b_coef, n));
legend('Exacte', 'Galerkin', 'Upwind (plein)', 'Upwind (optimal)', 'SUPG (optimal)', ...
       'Location', 'NorthWest');
grid on;
print('-dpng', '-r150', '../figures/q5b_upwind_supg_n40.png');

fprintf('Q5b done.\n');
fprintf('Commentaire:\n');
fprintf('- Galerkin standard: oscillations spurieuses pour Pe_h > 1.\n');
fprintf('- Upwind plein: stable mais excessivement diffusif (smearing de la couche limite).\n');
fprintf('- Upwind optimal: stable et nodalement exact pour 1D constant coefficients.\n');
fprintf('- SUPG optimal: equivalent a upwind optimal en 1D, donne la solution nodale exacte.\n');
