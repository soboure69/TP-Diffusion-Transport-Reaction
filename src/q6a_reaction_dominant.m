%% Q6a - Diffusion-Reaction dominant (b=0, sigma>>mu)
%  -mu u'' + sigma u = 0, u(0)=0, u(1)=1
%  Exact: u(x) = (e^{alpha x} - e^{-alpha x}) / (e^alpha - e^{-alpha})
%         = sinh(alpha x) / sinh(alpha),  alpha = sqrt(sigma/mu)
%  sigma=1, mu=1/6000, n=10,20,30,40

clear; close all;

mu = 1/6000; sigma_coef = 1; b_coef = 0;
fhandle = @(x) 0;
alpha = sqrt(sigma_coef/mu);

% exact solution
x_exact = linspace(0,1,1000)';
u_exact = sinh(alpha*x_exact) / sinh(alpha);

n_values = [10, 20, 30, 40];
colors = {'r', 'b', 'g', 'm'};
markers = {'o', 's', 'd', '^'};

figure();
plot(x_exact, u_exact, 'k-', 'LineWidth', 2); hold on;

for i = 1:length(n_values)
  n = n_values(i);
  [xn, uh] = solve_dtr(mu, b_coef, sigma_coef, n, fhandle, 'galerkin');
  plot(xn, uh, [colors{i} markers{i} '-'], 'MarkerSize', 4, 'LineWidth', 1);
  h = 1/n;
  Pe_react = sigma_coef * h^2 / (6*mu);  % reaction Peclet analog
  fprintf('n=%d, h=%.4f, sigma*h^2/(6*mu)=%.2f, max(uh)=%.4f, min(uh)=%.4f\n', ...
          n, h, Pe_react, max(uh), min(uh));
end

xlabel('x'); ylabel('u(x)');
title(sprintf('Q6a: Diffusion-Reaction (\\mu=1/%d, \\sigma=%g)', round(1/mu), sigma_coef));
legend('Exacte', sprintf('n=%d',n_values(1)), sprintf('n=%d',n_values(2)), ...
       sprintf('n=%d',n_values(3)), sprintf('n=%d',n_values(4)), 'Location', 'NorthWest');
grid on;
print('-dpng', '-r150', '../figures/q6a_reaction_dominant.png');

fprintf('\nQ6a done.\n');
fprintf('alpha = sqrt(sigma/mu) = %.2f\n', alpha);
fprintf('Commentaire:\n');
fprintf('- La solution exacte presente une couche limite en x=0 et un plateau a 0 avant\n');
fprintf('  de monter rapidement pres de x=1 vers u=1.\n');
fprintf('- Pour des maillages grossiers, la solution EF peut presenter des oscillations\n');
fprintf('  spurieuses lorsque sigma*h^2/(6*mu) > 1 (regime reaction-dominant).\n');
fprintf('- Avec n suffisamment grand (h petit), la solution EF converge vers la solution exacte.\n');
