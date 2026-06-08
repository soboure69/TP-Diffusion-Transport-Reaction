%% Q4a - Diffusion-Transport: mu=b=f=1, sigma=0
%  -u'' + u' = 1 on (0,1), u(0)=0, u(1)=1
%  Exact solution: u(x) = x

clear; close all;

mu = 1; b_coef = 1; sigma = 0;
fhandle = @(x) 1;
n = 20;

[xnodes, uh] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'galerkin');

% exact solution
x_exact = linspace(0,1,200)';
u_exact = x_exact;   % u(x)=x

figure();
plot(x_exact, u_exact, 'b-', 'LineWidth', 1.5); hold on;
plot(xnodes, uh, 'ro-', 'MarkerSize', 5);
xlabel('x'); ylabel('u(x)');
title(sprintf('Q4a: Diffusion-Transport (\\mu=%g, b=%g, f=%g, \\sigma=%g, n=%d)', ...
      mu, b_coef, 1, sigma, n));
legend('Solution exacte u(x)=x', 'Solution EF u_h', 'Location', 'NorthWest');
grid on;
print('-dpng', '-r150', '../figures/q4a_diffusion_transport.png');
fprintf('Q4a done. Max error = %e\n', max(abs(uh - xnodes)));
