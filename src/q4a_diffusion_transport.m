%% Q4a - Diffusion-Transport: mu=b=f=1, sigma=0
%  -u'' + u' = 1 on (0,1), u(0)=0, u(1)=1
%  Objectif: representer graphiquement la solution approchee u_h (EF P1).

clear; close all;
if exist('OCTAVE_VERSION','builtin'), graphics_toolkit('fltk'); warning('off','all'); end

mu = 1; b_coef = 1; sigma = 0;
fhandle = @(x) 1;
n = 20;

[xnodes, uh] = solve_dtr(mu, b_coef, sigma, n, fhandle, 'galerkin');

figure();
plot(xnodes, uh, 'ro-', 'MarkerSize', 5, 'LineWidth', 1.2);
xlabel('x'); ylabel('u_h(x)');
title(sprintf('Q4a: Solution approchee u_h (\\mu=%g, b=%g, f=%g, \\sigma=%g, n=%d)', ...
      mu, b_coef, 1, sigma, n));
legend('Solution EF u_h (P1)', 'Location', 'NorthWest');
grid on;
print('-dpng', '-r150', '../figures/q4a_diffusion_transport.png');

fprintf('Q4a done. n=%d, u_h(0)=%.4f, u_h(1)=%.4f, min=%.4f, max=%.4f\n', ...
        n, uh(1), uh(end), min(uh), max(uh));
