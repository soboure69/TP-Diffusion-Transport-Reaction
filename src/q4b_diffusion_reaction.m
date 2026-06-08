%% Q4b - Diffusion-Reaction: mu=sigma=1, b=0, f(x)=x
%  -u'' + u = x on (0,1), u(0)=0, u(1)=1
%  Exact solution: u(x) = x

clear; close all;

mu = 1; b_coef = 0; sigma = 1;
fhandle = @(x) x;
n = 20;

[xnodes, uh] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'galerkin');

% exact solution
x_exact = linspace(0,1,200)';
u_exact = x_exact;   % u(x)=x

figure();
plot(x_exact, u_exact, 'b-', 'LineWidth', 1.5); hold on;
plot(xnodes, uh, 'ro-', 'MarkerSize', 5);
xlabel('x'); ylabel('u(x)');
title(sprintf('Q4b: Diffusion-Reaction (\\mu=%g, \\sigma=%g, b=%g, f(x)=x, n=%d)', ...
      mu, sigma, b_coef, n));
legend('Solution exacte u(x)=x', 'Solution EF u_h', 'Location', 'NorthWest');
grid on;
print('-dpng', '-r150', '../figures/q4b_diffusion_reaction.png');
fprintf('Q4b done. Max error = %e\n', max(abs(uh - xnodes)));
