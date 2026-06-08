%% Q5a - Transport dominant (sigma=0, b>>mu): standard Galerkin
%  -mu u'' + b u' = 0, u(0)=0, u(1)=1
%  Exact: u(x) = (exp(b/mu * x) - 1)/(exp(b/mu) - 1)
%  b=1, mu=0.01, n=10,40,60

clear; close all;

mu = 0.01; b_coef = 1; sigma = 0;
fhandle = @(x) 0;
n_values = [10, 40, 60];

% exact solution
x_exact = linspace(0,1,500)';
alpha = b_coef/mu;
u_exact = (exp(alpha*x_exact) - 1) / (exp(alpha) - 1);

figure();
plot(x_exact, u_exact, 'k-', 'LineWidth', 2); hold on;
colors = {'r', 'b', 'g'};
markers = {'o', 's', 'd'};

for i = 1:length(n_values)
  n = n_values(i);
  [xnodes, uh] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'galerkin');
  plot(xnodes, uh, [colors{i} markers{i} '-'], 'MarkerSize', 4, 'LineWidth', 1);
  h = 1/n;
  Pe = b_coef*h/(2*mu);
  fprintf('n=%d, h=%.4f, Pe_h=%.2f, max(uh)=%.4f, min(uh)=%.4f\n', ...
          n, h, Pe, max(uh), min(uh));
end

xlabel('x'); ylabel('u(x)');
title(sprintf('Q5a: Transport dominant - Galerkin (\\mu=%.2f, b=%g)', mu, b_coef));
legend('Exacte', sprintf('n=%d',n_values(1)), sprintf('n=%d',n_values(2)), ...
       sprintf('n=%d',n_values(3)), 'Location', 'NorthWest');
grid on;
print('-dpng', '-r150', '../figures/q5a_transport_galerkin.png');
fprintf('\nQ5a done.\n');
fprintf('Oscillations appear when the element Peclet number Pe_h = b*h/(2*mu) > 1.\n');
fprintf('For b=%g, mu=%g: h_crit = 2*mu/b = %g (n_crit = %g).\n', ...
        b_coef, mu, 2*mu/b_coef, ceil(b_coef/(2*mu)));
